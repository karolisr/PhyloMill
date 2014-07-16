# -*- coding: utf-8 -*-

from __future__ import print_function

LOG_LINE_SEP = '################################################################################'

def blacklist_gis(gis, kr_seq_db_object, log_file_path):

    from krpy.krother import write_log

    DB = kr_seq_db_object
    LFP = log_file_path

    msg = 'Inactivating blacklisted gis.'
    write_log(msg, LFP, newlines_before=0, newlines_after=0)

    for gi in gis:

        GI = int(gi)

        in_blacklist = kr_seq_db_object.in_blacklist(
            record_reference=GI,
            record_reference_type='gi')

        if not in_blacklist:

            blacklist_notes = 'user_deleted'

            where_dict = {'ncbi_gi': GI}

            row_ids = DB.db_get_row_ids(table_name='records', where_dict=where_dict)

            if row_ids:

                msg = 'inactivating:' + \
                ' gi:' + str(GI) + ' note:' + blacklist_notes
                write_log(msg, LFP)

                record = DB.get_record(
                    record_reference=GI,
                    record_reference_type='gi'  # gi version internal raw
                    )

                del_rec_id = int(record.annotations['kr_seq_db_id'])

                DB.set_inactive(
                    table_name='records',
                    where_dict=where_dict)

                DB.add_record_to_blacklist(
                    ncbi_gi=GI,
                    ncbi_version=record.id,
                    internal_reference=record.annotations['internal_reference'],
                    notes=blacklist_notes)

                where_dict = {'parent_rec_id': del_rec_id}

                DB.db_delete(
                    table_name='record_ancestry',
                    where_dict=where_dict)

            else:

                DB.add_record_to_blacklist(
                    ncbi_gi=GI,
                    ncbi_version=None,
                    internal_reference=None,
                    notes=blacklist_notes)

            DB.save()


def download_new_records(locus_name, ncbi_db, gis, dnld_file_path, kr_seq_db_object, email, log_file_path):

    from krpy import krncbi
    # from krpy import krother
    from krpy import krbioio
    from krpy.krother import write_log

    # gis_clean = list()
    # for gi in gis:
    #     if gi not in LOCI[locus_name]['bad_gis']:
    #         gis_clean.append(gi)

    gis_clean = gis
    LFP = log_file_path
    DB = kr_seq_db_object
    # DNLD_DIR_PATH = dnld_dir_path
    EMAIL = email

    msg = 'Found ' + str(len(gis_clean)) + ' records.'
    write_log(msg, LFP)

    gis_in_blacklist = list()
    gis_good = list()

    for gi in gis_clean:

        in_blacklist = DB.in_blacklist(
            record_reference = gi,
            record_reference_type='gi')

        if in_blacklist:
            gis_in_blacklist.append(gi)
        else:
            gis_good.append(gi)

    msg = 'There are ' + str(len(gis_in_blacklist)) + \
          ' blacklisted records.'
    write_log(msg, LFP)

    # gis_in_db = list()
    # gis_new = list()

    all_gis_in_db = DB.get_all_record_ids(
        record_reference_type='gi')

    # print('all_gis_in_db', len(all_gis_in_db))
    # print('gis_good', len(gis_good))

    # gis_new = list()
    # for gi in gis_good:
    #     if gi not in all_gis_in_db:
    #         gis_new.append(gi)

    gis_new = list(set(gis_good) - set(all_gis_in_db))

    # for gi in gis_good:

    #     in_db = gi in all_gis_in_db

    #     if in_db:
    #         gis_in_db.append(gi)
    #     else:
    #         gis_new.append(gi)

    msg = 'There are ' + str(len(gis_new)) + \
          ' new records.'
    write_log(msg, LFP)

    records_new = list()
    if len(gis_new) > 0:

        msg = 'Downloading new records.'
        write_log(msg, LFP)

        # print('')

        krncbi.download_sequence_records(
            file_path=dnld_file_path,
            uids=gis_new,
            db=ncbi_db,
            entrez_email=EMAIL,
            verbose=True)

        # print('')

        records_new_temp = krbioio.read_sequence_file(
            file_path=dnld_file_path,
            file_format='gb',
            ret_type='dict',
            key='gi')

        downloaded_gis = records_new_temp.keys()
        downloaded_gis = [int(x) for x in downloaded_gis]

        gis_new_after_dnld = list(set(downloaded_gis) - set(all_gis_in_db))

        for gi in gis_new_after_dnld:
            records_new.append(records_new_temp[str(gi)])

        # print('gis_new_after_dnld', len(gis_new_after_dnld))

    return records_new


def produce_ncbi_query_string(ncbi_tax_ids, exclude_tax_ids, query_term_str, min_seq_length, max_seq_length):

    tax_ncbi_query_strings = list()
    for t in ncbi_tax_ids:
        tnqs = 'txid' + str(t) + '[Organism]'
        tax_ncbi_query_strings.append(tnqs)

    taxa_query_str = ' OR '.join(tax_ncbi_query_strings)
    taxa_query_str = ' AND (' + taxa_query_str + ')'

    tax_exclude_ncbi_query_strings = list()
    for te in exclude_tax_ids:
        tnqs = 'txid' + str(te) + '[Organism]'
        tax_exclude_ncbi_query_strings.append(tnqs)

    taxa_exclude_query_str = ' OR '.join(tax_exclude_ncbi_query_strings)
    taxa_exclude_query_str = ' NOT (' + taxa_exclude_query_str + ')'

    query_term_str = '(' + query_term_str + ')'
    seq_length_str = ''
    if max_seq_length > 0:
        seq_length_str = ' AND (' + str(min_seq_length) + ':' + str(max_seq_length) + '[Sequence Length])'

    query_str = query_term_str + \
                seq_length_str + \
                taxa_query_str + \
                taxa_exclude_query_str

    return query_str


def search_genbank(ncbi_db, query_term_str, ncbi_tax_ids, exclude_tax_ids, min_seq_length, max_seq_length, email, log_file_path):

    from krpy import krncbi

    from krpy.krother import write_log

    query_str = produce_ncbi_query_string(
        ncbi_tax_ids=ncbi_tax_ids,
        exclude_tax_ids=exclude_tax_ids,
        query_term_str=query_term_str,
        min_seq_length=min_seq_length,
        max_seq_length=max_seq_length
    )

    msg = 'ENTREZ query: ' + query_str
    write_log(msg, log_file_path, newlines_before=0, newlines_after=0,
        to_file=True, to_screen=False)

    result_uids = krncbi.esearch(esearch_terms=query_str, db=ncbi_db,
        email=email)

    uid_list = list()

    for uid in result_uids:
        gi = int(uid)
        uid_list.append(gi)

    return uid_list


def regular_search(kr_seq_db_object, log_file_path, email, loci, locus_name, ncbi_tax_ids, exclude_tax_ids, max_seq_length, dnld_dir_path):

    from krpy import krother
    # from krpy import krncbi
    # from krpy import krbioio
    from krpy import krcl
    from krpy.krother import write_log

    LOCI = loci
    LFP = log_file_path
    DB = kr_seq_db_object
    TAX_IDS = ncbi_tax_ids
    EXCLUDE_TAX_IDS = exclude_tax_ids
    MAX_SEQ_LENGTH = max_seq_length
    EMAIL = email
    DNLD_DIR_PATH = dnld_dir_path

    ncbi_db = LOCI[locus_name]['database']
    query_term_str = LOCI[locus_name]['query']

    min_seq_length = None
    strategies = LOCI[locus_name]['strategies']
    for s in strategies:
        l_msl = s['min_length']
        if not min_seq_length:
            min_seq_length = l_msl
        else:
            min_seq_length = min(l_msl, min_seq_length)

    write_log(LOG_LINE_SEP, LFP, newlines_before=0, newlines_after=0)

    msg = 'Searching NCBI ' + ncbi_db + ' database for ' + \
          locus_name + '.'
    write_log(msg, LFP, newlines_before=0, newlines_after=0)

    gis = search_genbank(
        ncbi_db=ncbi_db,
        query_term_str=query_term_str,
        ncbi_tax_ids=TAX_IDS,
        exclude_tax_ids=EXCLUDE_TAX_IDS,
        min_seq_length=min_seq_length,
        max_seq_length=MAX_SEQ_LENGTH,
        email=EMAIL,
        log_file_path=log_file_path)

    # gis_clean = gis

    timestamp = krother.timestamp()
    timestamp = timestamp.replace('-', '_')
    timestamp = timestamp.replace(':', '_')
    timestamp = timestamp.replace(' ', '_')

    gb_file_name = locus_name + '_' + timestamp + '.gb'
    gb_file_path = DNLD_DIR_PATH + gb_file_name

    records_to_add = download_new_records(
        locus_name=locus_name,
        ncbi_db=ncbi_db,
        gis=gis,
        dnld_file_path=gb_file_path,
        kr_seq_db_object=DB,
        email=EMAIL,
        log_file_path=LFP)

    if len(records_to_add) > 0:

        msg = 'Adding downloaded records to database. This may take a bit.'
        write_log(msg, LFP)

        record_count = len(records_to_add)

        # DB._DB_CURSOR.execute('BEGIN TRANSACTION;')

        for i, record in enumerate(records_to_add):

            krcl.print_progress(
                current=i+1, total=record_count, length=0,
                prefix=krother.timestamp() + ' ',
                postfix=' - ' + record.annotations['gi'],
                show_bar=False)

            DB.add_genbank_record(
                record=record,
                action_str='Genbank search result.')

            gi = int(record.annotations['gi'])

            DB.add_record_annotation(
                record_reference=gi,
                type_str='locus',
                annotation_str=locus_name,
                record_reference_type='gi')

            DB.add_record_annotation(
                record_reference=gi,
                type_str='locus_not_extracted',
                annotation_str=locus_name,
                record_reference_type='gi')

            DB.add_record_annotation(
                record_reference=gi,
                type_str='source_file',
                annotation_str=gb_file_name,
                record_reference_type='gi')

        # DB._DB_CURSOR.execute('END TRANSACTION;')

        print()

    DB.save()


def update_lineage_info(kr_seq_db_object, log_file_path, email):

    from krpy import krbionames
    from krpy import krncbi
    from krpy import krcl
    from krpy import krother
    from krpy.krother import write_log

    ########################################################################

    # Get common names and full lineage information

    LFP = log_file_path
    DB = kr_seq_db_object
    EMAIL = email

    write_log(LOG_LINE_SEP, LFP, newlines_before=0, newlines_after=0)

    msg = 'Updating lineage information, common names, and NCBI taxids.'
    write_log(msg, LFP, newlines_before=0, newlines_after=0)

    organisms = DB.get_organisms(where_dict={'active': 1, 'taxonomy_check_done': 0})

    if not organisms:
        return

    organism_count = len(organisms)
    tax_id_list = list()

    for org_dict in organisms:
        org_tax_id_list = org_dict['ncbi_tax_ids']
        tax_id_list = tax_id_list + org_tax_id_list

    # print('tax_id_list BEFORE', len(tax_id_list))

    tax_id_list = list(set(tax_id_list))

    # print('tax_id_list AFTER', len(tax_id_list))

    msg = 'Downloading lineage information.'
    write_log(msg, LFP, newlines_before=0, newlines_after=0)
    lineages = krncbi.get_lineages(email=EMAIL, tax_terms=None, tax_ids=tax_id_list)

    msg = 'Downloading common names.'
    write_log(msg, LFP, newlines_before=0, newlines_after=0)
    common_names = krncbi.get_common_names(email=EMAIL, tax_ids=tax_id_list)

    msg = 'Downloading NCBI taxids.'
    write_log(msg, LFP, newlines_before=0, newlines_after=0)
    all_org_tax_id_dict = krncbi.get_taxids(email=EMAIL, tax_terms=None, tax_ids=tax_id_list)

    # print('all_org_tax_id_dict', len(all_org_tax_id_dict.keys()))

    for i, org_dict in enumerate(organisms):

        org_id = org_dict['id']

        org_flat = krbionames.flatten_organism_name(
            parsed_name=org_dict, sep=' ')

        org_tax_id_list_temp = org_dict['ncbi_tax_ids']
        org_tax_id_list = list()

        for org_tax_id in org_tax_id_list_temp:
            org_tax_id_list = org_tax_id_list + all_org_tax_id_dict[str(org_tax_id)]

        org_tax_id_list = list(set(org_tax_id_list))

        org_tax_id_list_str = list()
        org_common_names = list()
        org_lineage = None
        for ti in org_tax_id_list:
            org_tax_id_list_str.append(str(ti))
            if str(ti) in common_names.keys():
                org_common_names.append(common_names[str(ti)])
            if str(ti) in lineages.keys():
                org_lineage = lineages[str(ti)]

        common_names_str = ''
        if org_common_names:
            common_names_str = ', '.join(org_common_names)
        else:
            common_names_str = org_flat

        org_tax_id_str = ', '.join(org_tax_id_list_str)

        print_str = org_flat + ' -> ' + common_names_str + ' (' + org_tax_id_str + ')'

        krcl.print_progress(
            current=i+1, total=organism_count, length=0,
            prefix=krother.timestamp() + ' ',
            postfix=' - ' + print_str + '\n',
            show_bar=False)

        msg = print_str
        write_log(msg, LFP, newlines_before=0, newlines_after=0,
        to_file=True, to_screen=False)

        DB.db_update(
            table_name='organisms',
            values_dict={'common_name': common_names_str},
            where_dict={'id': org_id})

        values_dict = {
            'genus': org_dict['genus'],
            'species': org_dict['species'],
            'subspecies': org_dict['subspecies'],
            'variety': org_dict['variety'],
            'hybrid': org_dict['hybrid'],
            'other': org_dict['other'],
            'authority': org_dict['authority'],
            'common_name': common_names_str
            }

        org_id_new = DB.add_organism(
            organism_dict=values_dict,
            taxonomy_list=org_lineage,
            ncbi_tax_id_list=org_tax_id_list,
            synonymy_check_done=org_dict['synonymy_check_done'],
            taxonomy_check_done=True)[0]

        DB.update_records(
            values_dict={'org_id': org_id_new},
            where_dict={'org_id': org_id})

    # DB._DB_CURSOR.execute('END TRANSACTION;')
    DB.save()
    DB.delete_orphaned_organisms()
    DB.save()
    DB.delete_orphaned_taxonomies()
    DB.save()


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
                # taxonomy_list = krncbi.get_lineage(
                #     email=email, tax_term=acc_name['genus'])
                taxonomy_list = krncbi.get_lineages(
                    email=email, tax_terms=[acc_name['genus']], tax_ids=None).values()[0]
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

                pass

    kr_seq_db_object.delete_orphaned_organisms()
    kr_seq_db_object.delete_orphaned_taxonomies()
    # kr_seq_db_object.save()


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
    merge_subspecies,
    merge_varieties,
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

    organisms = kr_seq_db_object.get_organisms(where_dict={'active': 1})
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
            prefix=krother.timestamp() + ' ',
            postfix=' - ' + org_flat,
            show_bar=False)

        tax_id = None
        tax_id_list = org_dict['ncbi_tax_ids']

        if tax_id_list:

            for tax_id in tax_id_list:
                if str(tax_id) in taxid_blacklist_set:

                    delete_note = 'blacklisted_taxid_' + str(tax_id)

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
        taxonomy_list_orig = None
        if taxonomy:
            taxonomy = taxonomy.split('$$$')
            taxonomy_list_orig = taxonomy
            taxonomy = krncbi.parse_lineage_string_list(taxonomy)[1]
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
                # not_in_synonymy = True
                # print('\n\n\nnot_in_synonymy', org_flat, '->', acc_name_flat, '\n\n\n')
            else:
                acc_name['status'] = 'synonymy'
                resolved = True
                deleted = False
                synonymy_check_done = True

        if (not (resolved or deleted)
            and merge_varieties
            and ((org_dict['variety'] is not None)
                and (org_dict['variety'] != ''))):

            org_dict['status'] = 'merge_varieties'
            org_dict['variety'] = None
            acc_name = org_dict
            acc_name_flat = krbionames.flatten_organism_name(
                parsed_name=acc_name, sep=' ')

            resolved = True
            deleted = False

        if (not (resolved or deleted)
            and merge_subspecies
            and ((org_dict['subspecies'] is not None)
                and (org_dict['subspecies'] != ''))):

            org_dict['status'] = 'merge_subspecies'
            org_dict['subspecies'] = None
            acc_name = org_dict
            acc_name_flat = krbionames.flatten_organism_name(
                parsed_name=acc_name, sep=' ')

            resolved = True
            deleted = False

        if resolved or synonymy_check_done:

            if org_flat != acc_name_flat:

                msg = 'renaming: ' + org_flat + ' -> ' + acc_name_flat + \
                ' tax_id:' + str(tax_id) + ' note:' + acc_name['status']
                write_log(msg, log_file_path)

                pass

            values_dict = {
                'genus': acc_name['genus'],
                'species': acc_name['species'],
                'subspecies': acc_name['subspecies'],
                'variety': acc_name['variety'],
                'hybrid': acc_name['hybrid'],
                'other': acc_name['other'],
                'authority': acc_name['authority'],
                'common_name': acc_name['common_name']
                }

            taxonomy_list = None
            if acc_name['genus'] in taxonomy_cache.keys():
                taxonomy_list = taxonomy_cache[acc_name['genus']]
            else:
                taxonomy_list = krncbi.get_lineages(
                    email=email, tax_terms=[acc_name['genus']], tax_ids=None).values()
                if taxonomy_list:
                    taxonomy_list = taxonomy_list[0]
                elif taxonomy_list_orig:
                    taxonomy_list = taxonomy_list_orig
                if taxonomy_list:
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

            msg = 'inactivating: ' + org_flat + \
            ' tax_id:' + str(tax_id) + ' note:' + delete_note
            write_log(msg, log_file_path)

            where_dict = {'org_id': org_id}

            rec_ids = kr_seq_db_object.db_get_row_ids(
                'records', where_dict=where_dict)

            blacklist_notes = delete_note + ' ' + org_flat

            kr_seq_db_object.set_inactive(
                table_name='records',
                where_dict=where_dict)

            if rec_ids:

                for rec_id in rec_ids:

                    rec = kr_seq_db_object.get_record(
                        record_reference=rec_id,
                        record_reference_type='raw'
                        )

                    # Have to check if the record is 'flattened', then it will
                    # not have an associated GI

                    rec_gi = rec.annotations['gi']
                    if rec_gi and rec_gi != 'None':
                        rec_gi = int(rec_gi)

                    kr_seq_db_object.add_record_to_blacklist(
                        ncbi_gi=rec_gi,
                        ncbi_version=rec.id,
                        internal_reference=rec.annotations['internal_reference'],
                        notes=blacklist_notes)

            where_dict = {'id': org_id}

            kr_seq_db_object.set_inactive(
                table_name='organisms',
                where_dict=where_dict)

        # kr_seq_db_object.save()

    kr_seq_db_object.delete_orphaned_organisms()
    kr_seq_db_object.delete_orphaned_taxonomies()
    # kr_seq_db_object.save()


def accept_records_by_similarity(records, seeds, identity_threshold=0.85, cpu=1):

    from krpy import kralign
    # from krpy import krusearch

    accept = list()
    reject = list()

    clusters = dict()

    if len(records) < 10:
        for r in records:
            accept.append(['+', r.annotations['gi'], 1.0])
    else:

        # if len(records) <= 500:

        aln_options = ''
        rec_lengths = list()
        for r in records:
            rec_lengths.append(len(r.seq))
        mean_seq_length = float(sum(rec_lengths)) / float(len(rec_lengths))
        if mean_seq_length >= 10000:
            aln_options = '--auto --nuc --reorder --adjustdirection --thread ' + str(cpu)
        else:
            aln_options = '--genafpair --jtt 10 --maxiterate 1000 --nuc --reorder --adjustdirection --thread ' + str(cpu)

        clusters = kralign.cluster(
            records=records,
            threshold=identity_threshold,
            unknown='N',
            key='gi',
            aln_program='mafft',
            aln_executable='mafft',
            aln_options=aln_options,
            seeds=seeds,
            seed_coverage=0.1,
            query_coverage=0.5)

        # else:

        #     clusters = krusearch.cluster_records(
        #         records=records,
        #         identity_threshold=identity_threshold,
        #         temp_dir=temp_dir,
        #         sorted_input=False,
        #         algorithm='smallmem',  # fast smallmem
        #         strand='both',  # plus both aa
        #         threads=1,
        #         quiet=True,
        #         program='usearch',
        #         heuristics=True,
        #         query_coverage=0.1,
        #         target_coverage=0.1,
        #         sizein=False,
        #         sizeout=False,
        #         usersort=False,
        #         seq_id='gi',
        #         cluster_key='centroid')

        min_clust_size = 5

        for key in clusters.keys():

            clust_size = len(clusters[key])

            if (key == 'unclustered') or (clust_size < min_clust_size):
                reject = reject + clusters[key]
            else:
                accept = accept + clusters[key]

            # if clust_size >= min_clust_size:
            #     accept = accept + clusters[key]
            # else:
            #     reject = reject + clusters[key]

            # print()
            # print()
            # print(key, clust_size)
            # print('--- --- --- --- ---')
            # print(clusters[key])
            # print('=== === === === ===')

    return {'accept': accept, 'reject': reject}


def feature_for_locus(record, feature_type, qualifier_label, locus_name_list,
    match_stringency='strict'):

    from krpy import krseq

    from Bio.SeqFeature import SeqFeature

    feature_indexes = list()
    loose_matching = False
    if match_stringency.lower() == 'loose':
        loose_matching = True
    for l in locus_name_list:
        fi = krseq.get_features_with_qualifier(
            record=record, qualifier_label=qualifier_label,
            qualifier=l.strip(), feature_type=feature_type,
            loose=loose_matching)
        feature_indexes = feature_indexes + fi

    feature_indexes = list(set(feature_indexes))
    feature_indexes.sort()
    feature = None
    log_message = ''

    if len(feature_indexes) == 0:
        log_message = 'No locus annotation.'
    elif len(feature_indexes) > 1:

        # print('MORE THAN ONE', record.annotations['gi'])

        loc_list_temp = list()
        # strand_list_temp = list()

        for fi in feature_indexes:
            # print("\n\t\tIndex: "+str(fi))
            # print("\n\t\t"+str(record.id))
            f = record.features[fi]
            if f.type == feature_type:
                for fq in f.qualifiers.keys():
                    if fq == qualifier_label:
                        # print("\t\t"+fq + ': ' + str(f.qualifiers[fq]) + ' ' + str(f.location))
                        loc_list_temp.append(int(f.location.start))
                        loc_list_temp.append(int(f.location.end))
                        # strand_list_temp.append(f.location.strand)

        # print(loc_list_temp)
        # print(strand_list_temp)

        location_string = '[' + str(min(loc_list_temp)) + ':' + str(max(loc_list_temp)) + '](' + '+' + ')'

        # print(location_string)

        location = krseq.location_from_string(
            location_string=location_string)

        feature = SeqFeature(
            location=location,
            type='fet_type',
            qualifiers=None)

        # ### THIS BREAKS IN A MULTITHREADED VERSION. While one thread waits for input
        # ### others print log messages and the pipeline crashes.

        # # Let user pick which of the indexes to use
        # print("\n\n\tFound more than one annotation for "+str(locus_name_list)+" please pick the correct index:")

        # for fi in feature_indexes:

        #     # print("\n\t\tIndex: "+str(fi))
        #     print("\n\t\t"+str(record.id))
        #     print()

        #     rng = range(1, 6)
        #     rng.reverse()
        #     for n in rng:
        #         if fi-n >= 0:
        #             f = record.features[fi-n]
        #             if f.type == feature_type:
        #                 for fq in f.qualifiers.keys():
        #                     if fq == qualifier_label:
        #                         print("\t\t"+fq + ': ' + str(f.qualifiers[fq]) + ' ' + str(f.location))

        #     print('\t\t---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----')

        #     f = record.features[fi]
        #     for fq in f.qualifiers.keys():
        #         if fq == qualifier_label:
        #             print("\t\t" + str(fi) + " > " + fq + ': ' + str(f.qualifiers[fq]) + ' ' + str(f.location))

        #     print('\t\t---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----')

        #     rng.reverse()
        #     for n in rng:
        #         if fi+n < len(record.features):
        #             f = record.features[fi+n]
        #             if f.type == feature_type:
        #                 for fq in f.qualifiers.keys():
        #                     if fq == qualifier_label:
        #                         print("\t\t"+fq + ': ' + str(f.qualifiers[fq]) + ' ' + str(f.location))

        #     print('\n\t\t==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====')

        # class BadChoiceException(Exception):
        #     pass

        # while True:
        #     picked_fi = raw_input("\n\tPick index or type 'exclude' to not use this sequence: ")
        #     # picked_fi = 1
        #     try:
        #         if str(picked_fi).lower().startswith('exclude'):
        #             print()
        #             log_message = 'More than one locus annotation. User excluded.'
        #             break
        #         else:
        #             try:
        #                 int(picked_fi)
        #             except ValueError:
        #                 print("\n\tBad choice.")
        #                 continue
        #             if int(picked_fi) in feature_indexes:
        #                 feature = record.features[int(picked_fi)]
        #                 print()
        #                 break
        #             else:
        #                 raise(BadChoiceException)
        #     except BadChoiceException:
        #         print("\n\tBad choice.")
    else:
        feature = record.features[feature_indexes[0]]

    return (feature, log_message)


def extract_loci(locus_dict, records, log_file_path, kr_seq_db_object, temp_dir, seed_recs_file_path, cpu=1):

    import os
    import copy

    from krpy import krcl
    from krpy import krother
    from krpy import kralign
    from krpy import krbioio
    from krpy.krother import write_log

    strategies = locus_dict['strategies']
    locus_name = locus_dict['name']
    # locus_short_name = locus_dict['short_name']

    acc_rej_gi_dict = dict()

    trimmed_records = list()

    no_feature_record_gi_list = list()
    short_feature_record_gi_list = list()

    records_count = len(records)

    prelim_record_feat_loc_list = list()

    msg = 'Extracting sequences for locus ' + locus_name + ' from downloaded records.'
    write_log(msg, log_file_path, newlines_before=0, newlines_after=0, to_file=True, to_screen=True)

    # if records_count < 3:
    #     return {'no_feature': list(), 'short_feature': list(), 'accept': list(), 'reject': list()}

    msg = 'There are ' + str(records_count) + ' records for locus ' + locus_name + '.'
    write_log(msg, log_file_path, newlines_before=0, newlines_after=0, to_file=True, to_screen=True)

    for i, record in enumerate(records):

        krcl.print_progress(
            current=i+1, total=records_count, length=0,
            prefix=krother.timestamp() + ' ',
            postfix='',
            show_bar=False)

        rec_id = int(record.annotations['kr_seq_db_id'])
        gi = int(record.annotations['gi'])

        location_list = list()

        feature_found = False

        for strategy in strategies:

            locus_relative_position = strategy['locus_relative_position']
            feature_type = strategy['feature_type']
            qualifier_label = strategy['qualifier_label']
            qualifier_value = strategy['qualifier_value']
            # l_re = strategy['regex']
            strict_value_match = strategy['strict_value_match']
            l_ml = strategy['min_length']
            # l_el = strategy['extra_length']

            match_stringency = 'loose'
            if strict_value_match:
                match_stringency = 'strict'

            start = None
            end = None
            strand = None

            if '->' in qualifier_value:

                locus_name_list = qualifier_value.split('->')
                locus_A = locus_name_list[0]
                locus_B = locus_name_list[1]

                # Locus may have different names
                locus_name_list_A = locus_A.split('|')
                locus_name_list_B = locus_B.split('|')

                feature_tuple_A = feature_for_locus(
                    record, feature_type, qualifier_label, locus_name_list_A,
                    match_stringency)
                feature_A = feature_tuple_A[0]

                feature_tuple_B = feature_for_locus(
                    record, feature_type, qualifier_label, locus_name_list_B,
                    match_stringency)
                feature_B = feature_tuple_B[0]

                if feature_A and feature_B:

                    feature_found = True

                    start_A = int(feature_A.location.start)
                    end_A = int(feature_A.location.end)
                    # strand_A = int(feature_A.location.strand)

                    start_B = int(feature_B.location.start)
                    end_B = int(feature_B.location.end)
                    # strand_B = int(feature_B.location.strand)

                    d1 = abs(start_A - start_B)
                    d2 = abs(end_A - start_B)
                    d3 = abs(start_A - end_B)
                    d4 = abs(end_A - end_B)

                    distances = (d1, d2, d3, d4)

                    min_distance = min(distances)
                    min_index = distances.index(min_distance)

                    if min_index == 0:
                        if start_B < start_A:
                            strand = -1
                        start = min(start_A, start_B)
                        end = max(start_A, start_B)
                    elif min_index == 1:
                        if start_B < end_A:
                            strand = -1
                        start = min(end_A, start_B)
                        end = max(end_A, start_B)
                    elif min_index == 2:
                        if end_B < start_A:
                            strand = -1
                        start = min(start_A, end_B)
                        end = max(start_A, end_B)
                    elif min_index == 3:
                        if end_B < end_A:
                            strand = -1
                        start = min(end_A, end_B)
                        end = max(end_A, end_B)

                    location_list.append((start, end, locus_relative_position, strand))

            else:
                # Locus may have different names
                locus_name_list = qualifier_value.split('|')

                feature_tuple = feature_for_locus(
                    record=record,
                    feature_type=feature_type,
                    qualifier_label=qualifier_label,
                    locus_name_list=locus_name_list,
                    match_stringency=match_stringency)

                feature = feature_tuple[0]

                if feature:

                    feature_found = True

                    # There should be only one matching index.
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    strand = int(feature.location.strand)

                    location_list.append((start, end, locus_relative_position, strand))


        if not feature_found:
            no_feature_record_gi_list.append(gi)

        else:

            location_list = sorted(location_list, key=lambda x: (x[2], x[0]), reverse=False)

            location_list_deduplicated = list()

            prev_lrp = None
            prev_loc_start = None
            prev_loc_end = None

            for loc in location_list:

                lrp = loc[2]
                loc_start = loc[0]
                loc_end = loc[1]
                strand = loc[3]

                if not prev_lrp:
                    prev_lrp = lrp
                    prev_loc_start = loc_start
                    prev_loc_end = loc_end
                    prev_strand = strand
                    continue

                if prev_lrp == lrp:
                    prev_loc_end = loc_end
                else:
                    location_list_deduplicated.append((prev_loc_start, prev_loc_end, prev_lrp, prev_strand))
                    prev_lrp = lrp
                    prev_loc_start = loc_start
                    prev_loc_end = loc_end
                    prev_strand = strand

            location_list_deduplicated.append((prev_loc_start, prev_loc_end, prev_lrp, prev_strand))

            location_list_deduplicated.append((location_list_deduplicated[0][0], location_list_deduplicated[-1][1], 'x', location_list_deduplicated[0][3]))
            location_list_deduplicated = sorted(location_list_deduplicated, key=lambda x: (x[2], x[0]), reverse=False)

            feat_dict = {'rec_id': rec_id, 'gi': gi, 'loc_list': location_list_deduplicated}
            prelim_record_feat_loc_list.append(feat_dict)

            loc = location_list_deduplicated[-1]
            trimmed_rec = record[loc[0]:loc[1]]
            if len(trimmed_rec.seq) < l_ml:
                short_feature_record_gi_list.append(gi)
            else:
                if loc[3] and loc[3] < 0:
                    trimmed_rec = trimmed_rec.reverse_complement()
                trimmed_rec.annotations['organism'] = record.annotations['organism']
                trimmed_rec.annotations['gi'] = record.annotations['gi']
                # trimmed_rec.id = record.id
                # trimmed_rec.name = record.name
                trimmed_records.append(trimmed_rec)

    msg = 'Filtering extracted records for locus ' + locus_name + '.'
    write_log(msg, log_file_path, newlines_before=0, newlines_after=0, to_file=True, to_screen=False)

    # Produce seeds. Longest records from each genus.

    seed_records = list()

    if os.path.exists(seed_recs_file_path):

        msg = '\tReading previously produced seed sequences for locus ' + locus_name + ' .'
        write_log(msg, log_file_path, newlines_before=0, newlines_after=0, to_file=True, to_screen=False)

        seed_records = krbioio.read_sequence_file(
            file_path=seed_recs_file_path,
            file_format='fasta',
            ret_type='list')

        for seed_rec in seed_records:
            seed_rec.annotations['gi'] = seed_rec.id

    elif len(trimmed_records) > 10:

        msg = '\tPreparing seed sequences for locus ' + locus_name + ' .'
        write_log(msg, log_file_path, newlines_before=0, newlines_after=0, to_file=True, to_screen=False)

        seed_records_dict = dict()
        for r in trimmed_records:
            genus = r.annotations['organism'].split(' ')[0]
            if genus not in seed_records_dict.keys():
                seed_records_dict[genus] = list()
            seed_records_dict[genus].append(r)

        seed_records_prelim = list()
        for key in seed_records_dict.keys():
            recs = seed_records_dict[key]
            recs = sorted(recs, key=lambda x: len(x.seq), reverse=True)
            seed_records_prelim = seed_records_prelim + recs[0:min(2, len(recs))]

        aln_options = ''
        rec_lengths = list()
        for r in seed_records_prelim:
            rec_lengths.append(len(r.seq))
        mean_seq_length = 0
        if rec_lengths:
            mean_seq_length = float(sum(rec_lengths)) / float(len(rec_lengths))
        if mean_seq_length >= 10000:
            aln_options = '--auto --nuc --reorder --adjustdirection --thread ' + str(cpu)
        else:
            aln_options = '--genafpair --jtt 10 --maxiterate 1000 --nuc --reorder --adjustdirection --thread ' + str(cpu)

        clusters = kralign.cluster(
            records=seed_records_prelim,
            threshold=0.90,
            unknown='N',
            key='gi',
            aln_program='mafft',
            aln_executable='mafft',
            aln_options=aln_options,
            seeds=None,
            seed_coverage=0.1,
            query_coverage=0.5)

        for key in clusters.keys():
            clust_size = len(clusters[key])
            if clust_size > 1:
                for sr in seed_records_prelim:
                    gi = sr.annotations['gi']
                    if gi == key:
                        seed_records.append(sr)

        seed_records = copy.copy(seed_records)

        for sr in seed_records:
            sr.description = ''
            sr.name = ''
            sr.id = sr.annotations['gi']

        krbioio.write_sequence_file(
            records=seed_records,
            file_path=seed_recs_file_path,
            file_format='fasta')

    msg = '\tFiltering records for locus ' + locus_name + '.'
    write_log(msg, log_file_path, newlines_before=0, newlines_after=0, to_file=True, to_screen=False)

    acc_rej_gi_dict = accept_records_by_similarity(
        records=trimmed_records,
        seeds=seed_records,
        identity_threshold=0.80,
        cpu=cpu)

    msg = 'Annotating sequence features in database for locus ' + locus_name + '.'
    write_log(msg, log_file_path, newlines_before=0, newlines_after=0, to_file=True, to_screen=False)

    acc_rej_gi_dict['no_feature'] = no_feature_record_gi_list
    acc_rej_gi_dict['short_feature'] = short_feature_record_gi_list

    # Fix sequence direction based on majority of good sequences
    acc_gi_list = acc_rej_gi_dict['accept']
    rev_comp_gi_list = list()
    for acc in acc_gi_list:

        kr_seq_db_object.remove_record_annotation(
            record_reference=int(acc[1]),
            type_str='locus_not_extracted',
            annotation_str=locus_name,
            record_reference_type='gi'  # gi version internal raw
            )

        if acc[0] == '-':
            rev_comp_gi_list.append(int(acc[1]))

    prfl_count = len(prelim_record_feat_loc_list)
    for i, r in enumerate(prelim_record_feat_loc_list):

        location_list_deduplicated = r['loc_list']
        rec_id = r['rec_id']
        gi = r['gi']

        krcl.print_progress(
            current=i+1, total=prfl_count, length=0,
            prefix=krother.timestamp() + ' ',
            postfix='',
            show_bar=False)

        for loc in location_list_deduplicated:

            strand_int = 1
            if gi in rev_comp_gi_list:
                strand_int = -1

            strand = ''
            if loc[3] and loc[3] * strand_int == 1:
                strand = '(+)'
            elif loc[3] and loc[3] * strand_int == -1:
                strand = '(-)'

            location_str = '[' + str(loc[0]) + ':' + str(loc[1]) + ']' + strand

            feat_type_str = 'pwf' + str(loc[2])

            feat_id = kr_seq_db_object.add_record_feature(
                rec_id=rec_id, type_str=feat_type_str,
                location_str=location_str)[0]

            qual_type_str = 'note'
            qual_str = locus_name + '|' + str(loc[2])

            kr_seq_db_object.add_record_feature_qualifier(
                rec_feat_id=feat_id,
                type_str=qual_type_str,
                qualifier_str=qual_str)

    kr_seq_db_object.save()

    return acc_rej_gi_dict


def trim_record_to_locus(record, locus_name):

    feat = feature_for_locus(
        record=record,
        feature_type='pwfx',
        qualifier_label='note',
        locus_name_list=[locus_name + '|x'],
        match_stringency='strict')[0]

    rec_trimmed = None

    if feat:

        start = feat.location.start
        end = feat.location.end
        strand = feat.location.strand

        rec_trimmed = record[start:end]

        if strand and (int(strand) < 0):
            rec_trimmed = rec_trimmed.reverse_complement()

        # rec_trimmed.id = record.id
        # rec_trimmed.name = record.name
        # rec_trimmed.description = record.description
        rec_trimmed.id = record.annotations[b'gi']
        rec_trimmed.name = ''
        rec_trimmed.description = ''

        rec_trimmed.annotations[b'gi'] = record.annotations[b'gi']
        # rec_trimmed.annotations[b'organism'] = record.annotations[b'organism']
        # rec_trimmed.annotations[b'kr_seq_db_org_id'] = record.annotations[b'kr_seq_db_org_id']
        # rec_trimmed.annotations[b'kr_seq_db_id'] = record.annotations[b'kr_seq_db_id']

    return rec_trimmed


def improve_alignment_using_reference_records(
    org_name,
    records,
    reference_records,
    locus_name,
    log_file_path,
    current_aln,
    aln_program='mafft',
    aln_program_executable='mafft',
    aln_options='--auto',
    min_locus_sequence_identity_range=[0.90, 0.95]):

    # import numpy

    from Bio.Align import MultipleSeqAlignment

    from Bio.SeqRecord import SeqRecord

    from krpy import kralign
    from krpy import krother
    from krpy.krother import write_log

    new_aln = None

    ident_step = 0.005

    bottom_ident = min_locus_sequence_identity_range[0]
    max_ident = min_locus_sequence_identity_range[1]
    current_ident = max_ident

    ref_alignments = list()
    ref_alignments.append(current_aln)

    for i, ref_rec in enumerate(reference_records):
        ref_rec_original_gi = ref_rec.annotations['gi']
        # ref_rec_trimmed = trim_record_to_locus(
        #     record=ref_rec, locus_name=locus_name)
        ref_rec_trimmed = SeqRecord(ref_rec.seq)
        ref_rec_trimmed.id = str(krother.random_id(20))
        ref_loc_records = records + [ref_rec_trimmed]

        ref_aln = kralign.align(
            records=ref_loc_records,
            program=aln_program,
            options=aln_options,
            program_executable=aln_program_executable)

        # ident = kralign.identity(
        #     alignment=ref_aln,
        #     unknown_letters=set(['N']),
        #     unknown_id=0.0,
        #     free_unknowns=True,
        #     gap_id=0.0,
        #     free_gaps=True,
        #     end_gap_id=0.0,
        #     free_end_gaps=True)

        ident_with_ref = kralign.identity(
            alignment=ref_aln,
            unknown_letters=set(['N']),
            free_unknowns=True,
            free_gaps=True,
            free_end_gaps=True,
            return_all_pairwise_identities=False)

        ref_aln_record_list = list()
        for r in ref_aln:
            if r.id != ref_rec_trimmed.id:
                ref_aln_record_list.append(r)
        ref_aln = MultipleSeqAlignment(ref_aln_record_list)

        ident_without_ref = kralign.identity(
            alignment=ref_aln,
            unknown_letters=set(['N']),
            free_unknowns=True,
            free_gaps=True,
            free_end_gaps=True,
            return_all_pairwise_identities=False)

        # print(ident_with_ref, ident_without_ref)

        ident = max(ident_with_ref, ident_without_ref)

        # print(ident)

        # ident = ident[1]
        # if ident:
        #     median_ident = numpy.median(ident)
        #     ident = median_ident
        #     # mean_ident = numpy.mean(ident)
        #     # std_ident = numpy.std(ident)
        # else:
        #     ident = ident[0]

        # print()
        msg = '\tLocus alignment identity: ' + str(ident) + ' (threshold=' + str(current_ident) + ') ref: ' + ref_rec_original_gi
        write_log(msg, log_file_path, newlines_before=0, newlines_after=0,
            to_file=True, to_screen=False)

        if ident < current_ident:

            ref_alignments.append([ref_aln, ident])
            ref_alignments = sorted(ref_alignments, key=lambda x: x[1], reverse=True)
            temp_ident = ref_alignments[0][1]
            if temp_ident >= current_ident:
                new_aln = ref_alignments[0]
                break
            current_ident = max(current_ident-ident_step, bottom_ident)

        else:

            new_aln = [ref_aln, ident]

            # Move the successful reference sequence to the top of the list
            reference_records.insert(0, reference_records.pop(i))

            break

    if not new_aln:
        new_aln = ref_alignments[0]

    # new_aln is a list: [aln, identity]
    return new_aln


def flatten_locus(
    org_name,
    records,
    reference_records,
    locus_dict,
    log_file_path,
    already_trimmed=False,
    aln_program='mafft',
    aln_program_executable='mafft',
    aln_options='--auto',
    min_locus_sequence_identity_range=[0.90, 0.98]):

    from krpy.krother import write_log
    from krpy import kralign
    # import numpy

    locus_name = locus_dict['name']

    records_trimmed = list()

    if already_trimmed:
        records_trimmed = records

    else:
        for record in records:

            rec_trimmed = trim_record_to_locus(record, locus_name)

            if rec_trimmed:
                records_trimmed.append(rec_trimmed)

    aln = None

    msg = org_name + ' - ' + str(len(records_trimmed)) + ' records.'
    write_log(msg, log_file_path, newlines_before=0, newlines_after=0,
        to_file=True, to_screen=False)

    if len(records_trimmed) > 1:

        aln = kralign.align(
            records=records_trimmed,
            program=aln_program,
            options=aln_options,
            program_executable=aln_program_executable
            )

        # ident = kralign.identity(
        #     alignment=aln,
        #     unknown_letters=set(['N']),
        #     unknown_id=0.0,
        #     free_unknowns=True,
        #     gap_id=0.0,
        #     free_gaps=True,
        #     end_gap_id=0.0,
        #     free_end_gaps=True)

        ident = kralign.identity(
            alignment=aln,
            unknown_letters=set(['N']),
            free_unknowns=True,
            free_gaps=True,
            free_end_gaps=True,
            return_all_pairwise_identities=False)

        # print(ident)

        # ident = ident[1]
        # if ident:
        #     median_ident = numpy.median(ident)
        #     ident = median_ident
        #     # mean_ident = numpy.mean(ident)
        #     # std_ident = numpy.std(ident)
        # else:
        #     ident = ident[0]

        # print()
        msg = '\tLocus alignment identity: ' + str(ident) + ' (threshold=' + str(min_locus_sequence_identity_range[1]) + ')'
        write_log(msg, log_file_path, newlines_before=0, newlines_after=0,
            to_file=True, to_screen=False)

        aln = [aln, ident]

        # if (ident < min_locus_sequence_identity_range[1]) and (len(records_trimmed) <= 40):
        if (ident < min_locus_sequence_identity_range[1]):

            new_aln = improve_alignment_using_reference_records(
                org_name=org_name,
                records=records_trimmed,
                reference_records=reference_records,
                locus_name=locus_name,
                log_file_path=log_file_path,
                current_aln=aln,
                aln_program=aln_program,
                aln_program_executable=aln_program_executable,
                aln_options=aln_options,
                min_locus_sequence_identity_range=min_locus_sequence_identity_range)

            aln = new_aln


    # aln is a list: [aln, identity]
    return aln


def update_record_alignment(rec_id, new_aln, aln_name, kr_seq_db_object):

    kr_seq_db_object.db_delete(
        table_name='alignments',
        where_dict={'rec_id': rec_id})

    seq_rep_id_list = list()

    for ar in new_aln:

        seq = kr_seq_db_object.get_sequence_for_record(
            record_reference=int(ar.id),
            record_reference_type='gi'
            )

        seq_rep_list = kr_seq_db_object.produce_seq_edits(
            s1=str(seq).upper(),
            s2=str(ar.seq).upper())

        ar_id = kr_seq_db_object.db_get_row_ids(
            table_name='records',
            where_dict={'ncbi_gi': int(ar.id)})[0]

        seq_id = kr_seq_db_object.db_get_row_ids(
            table_name='sequences',
            where_dict={'rec_id': ar_id})[0]

        seq_rep_id = kr_seq_db_object.add_sequence_representation(
            seq_id=seq_id,
            repr_list=seq_rep_list,
            rec_id=None,
            aln_id=None)[0]

        seq_rep_id_list.append(seq_rep_id)

    kr_seq_db_object.add_alignment(
        name=aln_name,
        seq_rep_id_list=seq_rep_id_list,
        description=None,
        rec_id=rec_id)


def produce_reference_sequences(locus_name, records, ref_recs_file_path, log_file_path, cpu=1):

    import os
    import sys

    import numpy

    from krpy import kralign
    from krpy import krbioio
    from krpy.krother import write_log

    LFP = log_file_path

    reference_records = list()

    if os.path.exists(ref_recs_file_path):

        msg = 'Reading previously produced reference sequences for locus ' + locus_name + '.'
        write_log(msg, LFP, newlines_before=0, newlines_after=0)

        reference_records = krbioio.read_sequence_file(
            file_path=ref_recs_file_path,
            file_format='fasta',
            ret_type='list')

        for ref_rec in reference_records:
            ref_rec.annotations['gi'] = ref_rec.id
    else:

        msg = 'Producing dereplicated set of reference sequences for locus ' + locus_name + ', this may take a bit.'
        write_log(msg, LFP, newlines_before=0, newlines_after=0)

        reference_records_temp = list()

        ref_rec_lengths = list()
        for record in records:
            rec_trimmed = trim_record_to_locus(
                record=record,
                locus_name=locus_name)
            reference_records_temp.append(rec_trimmed)
            ref_rec_lengths.append(len(rec_trimmed.seq))

        if len(reference_records_temp) == 0:
            msg = 'There are no annotated records for locus ' + locus_name + ". Was 'extract_loci' command run?"
            write_log(msg, LFP, newlines_before=0, newlines_after=0)
            sys.exit(1)

        ref_median_length = numpy.median(ref_rec_lengths)
        ref_mean_length = numpy.mean(ref_rec_lengths)
        # ref_std_length = numpy.std(ref_rec_lengths)
        # ref_cutoff_length = ref_mean_length - (1 * ref_std_length)
        # ref_cutoff_length = min(ref_median_length, ref_mean_length) - (0.5 * ref_std_length)
        ref_cutoff_length = min(ref_median_length, ref_mean_length)

        # print(len(ref_rec_lengths), ref_median_length, ref_std_length, ref_cutoff_length)

        for ref_rec in reference_records_temp:
            ref_rec_len = len(ref_rec.seq)
            if ref_rec_len >= ref_cutoff_length:
                unk_in_ref_rec = str(ref_rec.seq).upper().count('N')
                prop_unk = float(unk_in_ref_rec) / float(ref_rec_len)
                # print(prop_unk, ref_rec.annotations['gi'])
                if prop_unk < 0.1:
                    reference_records.append(ref_rec)

        if len(reference_records) > 200:
            reference_records = reference_records[0:200]

        # print(len(reference_records))

        # reference_records = sorted(reference_records, key=lambda x: len(x.seq), reverse=True)
        # reference_records = reference_records[0:min(100, len(reference_records))]

        aln_options = ''
        if ref_mean_length >= 10000:
            aln_options = '--auto --nuc --reorder --adjustdirection --thread ' + str(cpu)
        else:
            aln_options = '--genafpair --jtt 10 --maxiterate 1000 --nuc --reorder --adjustdirection --thread ' + str(cpu)

        reference_records = kralign.dereplicate(
            records=reference_records,
            threshold=0.90,
            unknown='N',
            key='gi',
            aln_program='mafft',
            aln_executable='mafft',
            aln_options=aln_options,
            seed_coverage=0.30,
            query_coverage=0.80)

        krbioio.write_sequence_file(
            records=reference_records,
            file_path=ref_recs_file_path,
            file_format='fasta')

    return reference_records


def get_outgroup_names(alignment):

    outgroup_taxa = list()
    for a in alignment:
        last_id_term = a.id.split('||')[-1]
        if last_id_term == 'OUT':
            outgroup_taxa.append(a.id)

    return outgroup_taxa


def produce_raxml_input_files(name, aln_file_path, out_dir_path, outgroup_taxa=None, partitions=None, threads=1, locus_name_list=None):

    import random
    from krpy import krio

    part_file_path = out_dir_path + 'locus-partitions-RAxML' + '.txt'
    if partitions and locus_name_list:
        f_part_raxml = open(part_file_path, 'wb')
        for i, part in enumerate(partitions):
            raxml_part_line = 'DNA, ' + locus_name_list[i] + ' = ' + str(part[0]) + '-' + str(part[1]) + '\n'
            f_part_raxml.write(raxml_part_line)
        f_part_raxml.close()

    rand_seed = str(random.randrange(0, 1000000000))
    raxml_dir = out_dir_path + 'RAxML_' + name
    raxml_commands_file = out_dir_path + 'RAxML_commands.txt'
    f_raxml = open(raxml_commands_file, 'wb')

    raxml_line_1 = 'raxml \\\n'
    f_raxml.write(raxml_line_1)

    raxml_line_2 = '-s ' + aln_file_path + ' \\\n'
    f_raxml.write(raxml_line_2)

    if partitions and locus_name_list:
        raxml_line_3 = '-q ' + part_file_path + ' \\\n'
        f_raxml.write(raxml_line_3)

    if outgroup_taxa:
        outgroup_taxa_raxml = ','.join(outgroup_taxa)
        raxml_line_4 = '-o "' + outgroup_taxa_raxml + '" \\\n'
        f_raxml.write(raxml_line_4)

    krio.prepare_directory(raxml_dir)

    raxml_line_5 = '-w ' + raxml_dir + ' \\\n'
    f_raxml.write(raxml_line_5)

    raxml_line_6 = '-m GTRCAT \\\n-j \\\n-T ' + str(threads) + ' \\\n-N 1 \\\n'
    f_raxml.write(raxml_line_6)

    raxml_line_7 = '-p ' + rand_seed + ' \\\n'
    f_raxml.write(raxml_line_7)

    raxml_line_8 = '-n ' + '' + name + '\n'
    f_raxml.write(raxml_line_8)

    f_raxml.close()


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
    from krpy import krbionames
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
        FROM record_feature_qualifiers
        LEFT OUTER JOIN record_features
        ON record_features.id IS rec_feat_id
        LEFT OUTER JOIN records
        ON records.id IS rec_id
        WHERE (record_feature_qualifiers.type IN ('specimen_voucher', 'cultivar', 'strain', 'isolate'))
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

        if not voucher:
            continue

        krcl.print_progress(
            current=i+1, total=tgrc_count, length=0,
            prefix=krother.timestamp() + ' ',
            postfix= ' - ' + voucher + ' - ' + org_flat,
            show_bar=False)

        tgrc_resolved = None
        if voucher in tgrc_cache.keys():
            tgrc_resolved = tgrc_cache[voucher]
        else:
            tgrc_resolved = resolve_org_name_using_tgrc(voucher)
            tgrc_cache[voucher] = tgrc_resolved

        # old_name = tgrc_resolved[0]
        # old_name_flat = krbionames.flatten_organism_name(
            # parsed_name=old_name, sep=' ')
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
            # taxonomy_list = krncbi.get_lineage(
            #     email=email, tax_term=new_name['genus'])
            taxonomy_list = krncbi.get_lineages(
                    email=email, tax_terms=[new_name['genus']], tax_ids=None).values()[0]
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

            pass

    kr_seq_db_object.delete_orphaned_organisms()
    kr_seq_db_object.delete_orphaned_taxonomies()
    # kr_seq_db_object.save()
