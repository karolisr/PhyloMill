# -*- coding: utf-8 -*-

from __future__ import print_function


def search_genbank(ncbi_db, query_term_str, ncbi_tax_ids, max_seq_length, email):

    from krpy import krncbi

    tax_ncbi_query_strings = list()
    for t in ncbi_tax_ids:
        tnqs = 'txid' + str(t) + '[Organism]'
        tax_ncbi_query_strings.append(tnqs)

    taxa_query_str = ' OR '.join(tax_ncbi_query_strings)
    taxa_query_str = '(' + taxa_query_str + ')'
    query_term_str = '(' + query_term_str + ')'
    seq_length_str = '0:' + str(max_seq_length) + '[Sequence Length]'

    query_str = taxa_query_str + ' AND ' + query_term_str + ' AND ' + \
                seq_length_str

    result_uids = krncbi.esearch(esearch_terms=query_str, db=ncbi_db,
        email=krncbi)

    uid_list = list()

    for uid in result_uids:
        gi = int(uid)
        uid_list.append(gi)

    return uid_list


def rename_organisms_with_record_taxon_mappings(
    kr_seq_db_object, record_taxon_mappings_dict, taxonomy_cache, log_file_path,
    email):

    from krpy import krbionames
    from krpy import krncbi
    from krpy.krother import write_log

    for gi_mapped in record_taxon_mappings_dict.keys():

        in_db = kr_seq_db_object.in_db(
            record_reference=int(gi_mapped),
            record_reference_type='gi')

        if in_db:

            rec = kr_seq_db_object.get_record(
                record_reference=int(gi_mapped), record_reference_type='gi')
            org_flat = rec.annotations['organism']
            acc_name_flat = record_taxon_mappings_dict[gi_mapped]

            acc_name = krbionames.parse_organism_name(
                name=acc_name_flat, sep=' ', ncbi_authority=False)

            acc_name['status'] = 'record_taxon_mapping'

            tax_id = krncbi.get_ncbi_tax_id_for_tax_term(
                email=email, tax_term=acc_name_flat)
            ncbi_tax_id_list = list()
            if tax_id:
                ncbi_tax_id_list = [tax_id]

            taxonomy_list = None
            if acc_name['genus'] in taxonomy_cache.keys():
                taxonomy_list = taxonomy_cache[acc_name['genus']]
            else:
                taxonomy_list = krncbi.get_lineage(
                    email=email, tax_term=acc_name['genus'])
                taxonomy_cache[acc_name['genus']] = taxonomy_list

            org_id_new = kr_seq_db_object.add_organism(
                organism_dict=acc_name,
                taxonomy_list=taxonomy_list,
                ncbi_tax_id_list=ncbi_tax_id_list)[0]

            kr_seq_db_object.update_records(
                values_dict={'org_id': org_id_new},
                where_dict={'ncbi_gi': int(gi_mapped)})

            if org_flat != acc_name_flat:

                msg = 'renaming: ' + org_flat + ' -> ' + acc_name_flat + \
                ' tax_id:' + str(tax_id) + ' gi:' + str(gi_mapped) + \
                ' note:' + acc_name['status']
                write_log(msg, log_file_path)

    kr_seq_db_object.delete_orphaned_organisms()
    kr_seq_db_object.delete_orphaned_taxonomies()
    kr_seq_db_object.save()


def resolve_org_name_using_synonymy(
    tax_id,
    organism,
    ncbi_names_table,
    synonymy_table,
    auth_file
    ):

    from krpy import krbionames

    acc_name = None
    source_name = krbionames.parse_organism_name(
        organism,
        sep=' ',
        ncbi_authority=True)

    resolved = krbionames.resolve_taxid(
        tax_id,
        ncbi_names_table,
        synonymy_table,
        auth_file,
        sorting='authority')

    resolved_list = list()
    if resolved:
        resolved_list = resolved[2]

    genus_set = set()
    species_set = set()

    for n in resolved_list:
        genus_set.add(n[1]['genus'])
        species_set.add(n[1]['species'])

    if len(resolved_list) == 0:
        acc_name = source_name
    elif len(genus_set) > 1 or len(species_set) > 1:
        acc_name = source_name
        acc_name['status'] = 'collision'
    else:
        acc_name = resolved[0]
        source_name = resolved[1]

    if acc_name['hybrid'] != '':
        acc_name['status'] = 'hybrid'

    if not acc_name['genus'] or acc_name['genus'] == '':
        acc_name['status'] = 'unknown'

    return (acc_name, source_name)


def rename_organisms_using_taxids(
    kr_seq_db_object,
    taxid_blacklist_set,
    taxid_taxon_mappings_dict,
    taxonomy_cache,
    log_file_path,
    email,
    rem_hybrids,
    rem_nonspecific,
    skip_prev_syn,
    tax_groups_to_syn,
    authority_alternates_file,
    ncbi_names_table,
    synonymy_table):

    from krpy import krbionames
    from krpy import krncbi
    from krpy import krcl
    from krpy import krother
    from krpy.krother import write_log

    organisms = kr_seq_db_object.get_organisms()
    organism_count = len(organisms)

    for i, org_dict in enumerate(organisms):

        resolved = False
        deleted = False
        not_in_synonymy = False

        delete_note = ''

        acc_name = None
        acc_name_flat = None

        org_id = org_dict['id']
        org_synonymy_check_done = org_dict['synonymy_check_done']

        org_flat = krbionames.flatten_organism_name(
            parsed_name=org_dict, sep=' ')

        krcl.print_progress(
            current=i+1, total=organism_count, length=0,
            prefix=krother.timestamp() + ' - ',
            postfix=' - ' + org_flat,
            show_bar=False)

        tax_id = None
        tax_id_list = org_dict['ncbi_tax_ids']

        if tax_id_list:

            for tax_id in tax_id_list:
                if str(tax_id) in taxid_blacklist_set:

                    delete_note = 'blacklisted_taxid'

                    resolved = False
                    deleted = True

                    break

            tax_id = tax_id_list[0]

        if (not deleted) and (str(tax_id) in taxid_taxon_mappings_dict.keys()):
            acc_name_flat = taxid_taxon_mappings_dict[str(tax_id)]
            acc_name = krbionames.parse_organism_name(
                name=acc_name_flat, sep=' ', ncbi_authority=False)
            acc_name['status'] = 'taxid_taxon_mapping'
            resolved = True
            deleted = False

        if (not (resolved or deleted)
            and rem_hybrids
            and ((org_dict['hybrid'] is not None)
                and (org_dict['hybrid'] != ''))):

            delete_note = 'hybrid'

            resolved = False
            deleted = True

        if (not (resolved or deleted)
            and rem_nonspecific
            and ((org_dict['species'] is None) or (org_dict['species'] == ''))):

            delete_note = 'unknown_species'

            resolved = False
            deleted = True

        check_syn = False
        taxonomy = org_dict['taxonomy']
        if taxonomy:
            taxonomy = taxonomy.split(',')
            taxonomy_lower = [x.lower() for x in taxonomy]
            for tax_term in tax_groups_to_syn:
                if tax_term in taxonomy_lower:
                    check_syn = True
                    break

        skip = False
        if skip_prev_syn and org_synonymy_check_done:
            skip = True

        synonymy_check_done = False

        if (not skip) and (not (resolved or deleted)) and check_syn:

            acc_name = resolve_org_name_using_synonymy(
                tax_id=tax_id,
                organism=org_flat,
                ncbi_names_table=ncbi_names_table,
                synonymy_table=synonymy_table,
                auth_file=authority_alternates_file)

            acc_name = acc_name[0]
            acc_name_flat = krbionames.flatten_organism_name(
                parsed_name=acc_name, sep=' ')

            if (acc_name['genus'] is None) or (acc_name['genus'] == ''):
                resolved = False
                deleted = False
                not_in_synonymy = True
                print('\n\n\nnot_in_synonymy', org_flat, '->', acc_name_flat, '\n\n\n')
            else:
                acc_name['status'] = 'synonymy'
                resolved = True
                deleted = False
                synonymy_check_done = True

        if resolved or synonymy_check_done:

            if org_flat != acc_name_flat:

                msg = 'renaming: ' + org_flat + ' -> ' + acc_name_flat + \
                ' tax_id:' + str(tax_id) + ' note:' + acc_name['status']
                write_log(msg, log_file_path)

            values_dict = {
                'genus': acc_name['genus'],
                'species': acc_name['species'],
                'subspecies': acc_name['subspecies'],
                'variety': acc_name['variety'],
                'hybrid': acc_name['hybrid'],
                'other': acc_name['other'],
                'authority': acc_name['authority']
                }

            taxonomy_list = None
            if acc_name['genus'] in taxonomy_cache.keys():
                taxonomy_list = taxonomy_cache[acc_name['genus']]
            else:
                taxonomy_list = krncbi.get_lineage(
                    email=email, tax_term=acc_name['genus'])
                taxonomy_cache[acc_name['genus']] = taxonomy_list

            ncbi_tax_id_list = None
            if tax_id:
                ncbi_tax_id_list=[tax_id]

            org_id_new = kr_seq_db_object.add_organism(
                organism_dict=values_dict,
                taxonomy_list=taxonomy_list,
                ncbi_tax_id_list=ncbi_tax_id_list,
                synonymy_check_done=synonymy_check_done)[0]

            kr_seq_db_object.update_records(
                values_dict={'org_id': org_id_new},
                where_dict={'org_id': org_id})

            kr_seq_db_object.delete_orphaned_organisms()

        elif deleted:

            msg = 'tax_id:' + str(tax_id) + ', status:deleting : ' + org_flat

            msg = 'deleting: ' + org_flat + \
            ' tax_id:' + str(tax_id) + ' note:' + delete_note

            write_log(msg, log_file_path)

            where_dict = {'org_id': org_id}
            blacklist_notes = delete_note + ' ' + org_flat
            kr_seq_db_object.delete_records(
                where_dict=where_dict, blacklist=True,
                blacklist_notes=blacklist_notes)

            where_dict = {'id': org_id}
            kr_seq_db_object.delete_organisms(where_dict=where_dict)

            kr_seq_db_object.delete_orphaned_taxonomies()

        kr_seq_db_object.save()

    kr_seq_db_object.delete_orphaned_organisms()
    kr_seq_db_object.delete_orphaned_taxonomies()
    kr_seq_db_object.save()


# Hacks! -----------------------------------------------------------------------


def clean_tgrc_voucher(raw_voucher):

    import re

    voucher = str(raw_voucher)

    voucher = re.findall('LA\d+', voucher)
    if len(voucher) == 0:
        return(None)
    voucher = voucher[0]
    v_split_2 = re.findall('(\d+|[a-zA-Z]+)', voucher)
    if len(v_split_2) > 1:
        voucher = v_split_2[0] + "%04d" % (int(v_split_2[1]),)

    return voucher


def resolve_org_name_using_tgrc(voucher):

    '''
        Some Lycopersicon / Solanum members have voucher numbers which resolve
        most current species names. This function will take a voucher number
        and will return a list:

            [OLD_NAME, NEW_NAME]
    '''

    import urllib2
    from bs4 import BeautifulSoup
    import krbionames
    import time
    import sys

    # Let's not overwhelm the server.
    time.sleep(0.5)

    url = 'http://tgrc.ucdavis.edu/Data/Acc/AccDetail.aspx?AccessionNum='

    response = None

    i = 0
    while True:
        try:
            response = urllib2.urlopen(url + voucher)
        except:
            i = i + 1
            if i == 10:
                break
            time.sleep(2 * i)
            continue
        break

    html = ''
    if response:
        html = response.read()
    else:
        # HTTP Error
        sys.exit(1)
    doc = BeautifulSoup(html)

    old_line = doc.find(id='lTaxon')
    if old_line.font is None:
        return(None)
    old = str(old_line.font.text).strip()
    new_line = doc.find(id='lTaxon2')
    new = str(new_line.font.text).strip()

    o_parsed = krbionames.parse_organism_name(
        name=old,
        sep=' ',
        ncbi_authority=False)

    n_parsed = krbionames.parse_organism_name(
        name=new,
        sep=' ',
        ncbi_authority=False)

    if o_parsed['genus'] == 'L.':
        o_parsed['genus'] = 'Lycopersicon'
    if o_parsed['genus'] == 'S.':
        o_parsed['genus'] = 'Solanum'
    if n_parsed['genus'] == 'L.':
        n_parsed['genus'] = 'Lycopersicon'
    if n_parsed['genus'] == 'S.':
        n_parsed['genus'] = 'Solanum'

    return [o_parsed, n_parsed]


def get_tgrc_gis(kr_seq_db_object):

    kr_seq_db_conn=kr_seq_db_object.get_cursor()

    sql_string = '''
        SELECT records.ncbi_gi, qualifier
        FROM record_feature_qualifier_types
        LEFT OUTER JOIN record_feature_qualifiers
        ON record_feature_qualifier_types.id IS record_feature_qualifiers.rec_feat_qual_type_id
        LEFT OUTER JOIN record_features
        ON record_features.id IS rec_feat_id
        LEFT OUTER JOIN records
        ON records.id IS rec_id
        WHERE (type IN ('specimen_voucher', 'cultivar', 'strain', 'isolate'))
        AND ((qualifier like 'LA___') OR (qualifier like 'LA____') OR (qualifier like 'LA_____') OR (qualifier like 'LA______'));
        '''

    kr_seq_db_conn.execute(sql_string)
    results = kr_seq_db_conn.fetchall()

    return results


def rename_tgrc_organisms(kr_seq_db_object, taxonomy_cache, log_file_path,
        email):

    from krpy import krbionames
    from krpy import krncbi
    from krpy import krcl
    from krpy import krother
    from krpy.krother import write_log

    tgrc_list = get_tgrc_gis(kr_seq_db_object=kr_seq_db_object)

    tgrc_cache = dict()

    tgrc_count = len(tgrc_list)

    for i, tgrc in enumerate(tgrc_list):

        gi = tgrc['ncbi_gi']

        rec = kr_seq_db_object.get_record(
            record_reference=gi, record_reference_type='gi')
        org_flat = rec.annotations['organism']

        raw_voucher = tgrc['qualifier']
        voucher = clean_tgrc_voucher(raw_voucher)

        krcl.print_progress(
            current=i+1, total=tgrc_count, length=0,
            prefix=krother.timestamp() + ' - ',
            postfix= ' - ' + voucher + ' - ' + org_flat,
            show_bar=False)

        tgrc_resolved = None
        if voucher in tgrc_cache.keys():
            tgrc_resolved = tgrc_cache[voucher]
        else:
            tgrc_resolved = resolve_org_name_using_tgrc(voucher)
            tgrc_cache[voucher] = tgrc_resolved

        old_name = tgrc_resolved[0]
        old_name_flat = krbionames.flatten_organism_name(
            parsed_name=old_name, sep=' ')
        new_name = tgrc_resolved[1]
        new_name_flat = krbionames.flatten_organism_name(
            parsed_name=new_name, sep=' ')

        tax_id = krncbi.get_ncbi_tax_id_for_tax_term(
            email=email, tax_term=new_name_flat)
        ncbi_tax_id_list = list()
        if tax_id:
            ncbi_tax_id_list = [tax_id]

        taxonomy_list = None
        if new_name['genus'] in taxonomy_cache.keys():
            taxonomy_list = taxonomy_cache[new_name['genus']]
        else:
            taxonomy_list = krncbi.get_lineage(
                email=email, tax_term=new_name['genus'])
            taxonomy_cache[new_name['genus']] = taxonomy_list

        org_id_new = kr_seq_db_object.add_organism(
            organism_dict=new_name,
            taxonomy_list=taxonomy_list,
            ncbi_tax_id_list=ncbi_tax_id_list)[0]

        kr_seq_db_object.update_records(
            values_dict={'org_id': org_id_new},
            where_dict={'ncbi_gi': gi})

        if org_flat != new_name_flat:

            msg = 'renaming: ' + org_flat + ' -> ' + new_name_flat + \
            ' tax_id:' + str(tax_id) + ' gi:' + str(gi) + ' voucher:' + \
            str(voucher) + ' note:tgrc'
            write_log(msg, log_file_path)

    kr_seq_db_object.delete_orphaned_organisms()
    kr_seq_db_object.delete_orphaned_taxonomies()
    kr_seq_db_object.save()
