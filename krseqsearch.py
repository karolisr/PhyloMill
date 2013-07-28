#!/usr/bin/env python

from __future__ import print_function
#from __future__ import unicode_literals


# Hacks! ----------------------------------------------------------------------

def hack_tgrc(voucher):

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
    import re

    # Let's not overwhelm the server.
    time.sleep(1)

    url = 'http://tgrc.ucdavis.edu/Data/Acc/AccDetail.aspx?AccessionNum='

    voucher = str(voucher)

    # More hacking, based on experience ---------------------------------------
    voucher = re.findall('LA\d+', voucher)
    if len(voucher) == 0:
        return(None)
    voucher = voucher[0]
    v_split_2 = re.findall('(\d+|[a-zA-Z]+)', voucher)
    if len(v_split_2) > 1:
        voucher = v_split_2[0] + "%04d" % (int(v_split_2[1]),)
    # -------------------------------------------------------------------------

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

    return([o_parsed, n_parsed, voucher])


# Pipeline functions ----------------------------------------------------------

def search_and_download(queries, output_dir, file_name_sep, email, input_dir, log_dir):

    '''
    This will search NCBI and download sequences.
    '''

    import os
    import krio
    import krncbi
    import krbioio

    ps = os.path.sep

    print('\nSearching NCBI.')
    print('\n\tPreparing output directory "', output_dir, '"', sep='')

    krio.prepare_directory(output_dir)

    all_uids = set()

    for query_dict in queries:

        name1 = query_dict['name1']
        name2 = query_dict['name2']
        db = query_dict['database']
        query = query_dict['query']

        '''
        File name is generated based on the search query (periods will be
        replaced with file_name_sep):

            name1.name2.locus.minlen.feature_type.qualifier_label
            .database.[other].extension

        Short explanation:

            Full name of a query consists of a name1 and a name2. This
            is to allow for multiple query strings for the same locus.
            An example is our good friend trnK/matK. trnK contains matK,
            so we want to be sure we get all trnK and matK loci because it
            may happen that matK will not be annotated within trnK. While
            it is possible to introduce "OR" statements into NCBI queries
            and collapse the search query into one string, it does not help
            us later if we wish to treat the results independently.
        '''

        file_path = (output_dir.rstrip(ps) + ps +
                     name1 + file_name_sep +
                     name2 + file_name_sep +
                     '.gb')

        if query.lower().startswith('donotdownload'):
            all_records = dict()
            file_list = krio.parse_directory(output_dir, file_name_sep)
            for f in file_list:
                if not f['ext'].startswith('gb'):
                    continue
                if f['split'][0] == name1:
                    records = krbioio.read_sequence_file(
                        file_path=f['path'],
                        file_format='genbank',
                        ret_type='dict')
                    for key in records.keys():
                        if records[key] not in all_records.keys():
                            all_records[key] = records[key]

            all_records = all_records.values()
            print('\n\tWill not download:', name1, name2, '\n', sep=' ')
            krbioio.write_sequence_file(records, file_path, 'genbank')

        else:

            print('\n\tSearching for:', name1, name2, 'in', db, 'database.\n',
                  sep=' ')

            # Search NCBI.
            result_uids = krncbi.esearch(query, db, email)
            all_uids |= set(result_uids)
            # Download records.
            krncbi.download_sequence_records(file_path, result_uids, db, email)

    downloaded_gis_file = input_dir.rstrip(ps) + ps + 'downloaded-gis.csv'

    # First find previously downloaded GIs and see how many of them are the same
    previous_uids = krio.read_table_file(
        path=downloaded_gis_file,
        has_headers=False,
        headers=None,
        delimiter=',',
        quotechar=None,
        stripchar='"',
        rettype='set')

    novel_gis_file = log_dir.rstrip(ps) + ps + '01-novel-gis.csv'
    missing_gis_file = log_dir.rstrip(ps) + ps + '01-missing-gis.csv'

    novel_uids = all_uids - previous_uids
    missing_uids = previous_uids - all_uids

    novel_uids = list(novel_uids)
    handle = open(novel_gis_file, 'w')
    for uid in novel_uids:
        handle.write(uid+'\n')
    handle.close()

    missing_uids = list(missing_uids)
    handle = open(missing_gis_file, 'w')
    for uid in missing_uids:
        handle.write(uid+'\n')
    handle.close()

    all_uids = list(all_uids)
    handle = open(downloaded_gis_file, 'w')
    for uid in all_uids:
        handle.write(uid+'\n')
    handle.close()


def filter_records(search_results_dir, output_dir, cutlist_records_file,
                   cutlist_taxonomy_file):

    import os
    import krio
    import krbioio
    import krcl
    import krseq
    import krncbi

    ps = os.path.sep

    print('\nFiltering records.')

    print('\n\tPreparing output directory "', output_dir, '"', sep='')
    krio.prepare_directory(output_dir)

    cutlist_records = None
    cutlist_taxonomy = None

    if cutlist_records_file:

        cutlist_records = krio.read_table_file(
            path=cutlist_records_file,
            has_headers=False,
            headers=None,
            delimiter=',',
            quotechar='"',
            rettype='set')

    if cutlist_taxonomy_file:

        cutlist_taxonomy = krio.read_table_file(
            path=cutlist_taxonomy_file,
            has_headers=False,
            headers=None,
            delimiter=',',
            quotechar='"',
            rettype='set')

    file_list = krio.parse_directory(search_results_dir, ' ')

    # Iterate over search results
    for f in file_list:

        if not f['ext'].startswith('gb'):
            continue

        output_file = output_dir + ps + f['full']

        # Read search results
        records = krbioio.read_sequence_file(f['path'], 'genbank')

        records_filtered_id = list()
        records_filtered_tax = list()

        print('\n\tProcessing: ', f['full'], sep='')

        krcl.hide_cursor()

        print('\n\tFiltering records for ids.')
        if cutlist_records:
            records_count = len(records)
            for i, record in enumerate(records):
                krcl.print_progress(i + 1, records_count, 50, '\t')
                if cutlist_records:
                    if ((record.id not in cutlist_records) and
                        (krseq.get_annotation(record, 'gi') not in
                            cutlist_records)):
                        records_filtered_id.append(record)
        else:
            records_filtered_id = records
        print('\n\tAccepted', len(records_filtered_id), 'records.')

        print('\n\tFiltering records for taxonomy.')
        if cutlist_taxonomy:
            records_count = len(records_filtered_id)
            for i, record in enumerate(records_filtered_id):
                krcl.print_progress(i + 1, records_count, 50, '\t')
                if cutlist_taxonomy:
                    if ((krncbi.get_ncbi_tax_id(record) not in
                        cutlist_taxonomy) and
                        krseq.get_annotation(record, 'organism') not in
                            cutlist_taxonomy):
                        records_filtered_tax.append(record)
        else:
            records_filtered_tax = records_filtered_id
        print('\n\tAccepted', len(records_filtered_tax), 'records.')

        # Write results
        if cutlist_taxonomy:
            krbioio.write_sequence_file(records_filtered_tax, output_file,
                                        'gb')
        else:
            krbioio.write_sequence_file(records_filtered_id, output_file, 'gb')

        krcl.show_cursor()


def produce_seq_id(locus, reference_id, taxid, old_name, new_name, separator_0='|'):
    sequence_record_id = (
        locus + separator_0 +
        # NCBI accession
        reference_id + separator_0 +
        # Organism name as it appears in search results
        taxid + separator_0 +
        # NCBI taxid
        old_name + separator_0 +
        # Accepted name
        new_name)
    return(sequence_record_id)


separator_lrps = '%'
separator_lrp = '_'
separator_unique = '@'
separator_accession = '$'


def parse_seq_id(seq_id, separator_0='|'):
    parsed = dict()
    seq_id_split_0 = seq_id.split(separator_0)
    parsed['locus'] = seq_id_split_0[0]
    parsed['reference_id'] = seq_id_split_0[1]
    parsed['taxid'] = seq_id_split_0[2]
    parsed['old_name'] = seq_id_split_0[3]
    parsed['new_name'] = seq_id_split_0[4]
    reference_id_split_lrps = parsed['reference_id'].split(separator_lrps)
    reference_id_split_unique = parsed['reference_id'].split(separator_unique)
    lrps = False
    if len(reference_id_split_lrps) > 1:
        lrps = reference_id_split_lrps[0].split(separator_lrp)
    parsed['lrps'] = lrps
    unique = False
    if len(reference_id_split_unique) > 1:
        unique = reference_id_split_unique[1]
    parsed['unique'] = unique
    accessions = False
    if lrps and unique:
        accessions = parsed['reference_id'].split(separator_lrps)[1].split(separator_unique)[0].split(separator_accession)
    elif lrps:
        accessions = parsed['reference_id'].split(separator_lrps)[1].split(separator_accession)
    elif unique:
        accessions = parsed['reference_id'].split(separator_unique)[0].split(separator_accession)
    else:
        accessions = parsed['reference_id'].split(separator_accession)
    parsed['accessions'] = accessions
    return(parsed)


def extract_loci(
    search_results_dir,
    output_dir,
    queries,
    ncbi_names_table,
    temp_dir,
    file_name_sep,
    synonymy_table=None,
    auth_file=None,
    hacks=None,
    hacks_data_location=None,
    unresolvable_taxonomy_list=None,
    keeplist_taxonomy_list=None,
    force_reverse_complement_list=None,
    log_dir=None
):

    '''
    Extract relevant loci from the search results, do some filtering by length
    and similarity.
    '''

    print('\nExtracting relevant loci.')

    import os
    import shutil

    from Bio import SeqRecord

    import krio
    import krbioio
    import krseq
    import krncbi
    import krcl
    import krbionames

    ps = os.path.sep

    # HACKS ###################################################################

    hack_sol_genus = None
    hack_sol_species = None

    found_previously = dict()
    do_not_repeat = list()

    if hacks and 'tgrc' in hacks:
        hack_sol_species_set = krio.read_table_file(
            path=hacks_data_location['tgrc'],
            has_headers=False,
            headers=None,
            delimiter=',',
            quotechar='"',
            stripchar='',
            commentchar="#",
            rettype='set'  # dict, list, set
        )

        hack_sol_genus = list()
        hack_sol_species = list()

        for r in hack_sol_species_set:
            name = krbionames.parse_organism_name(
                r,
                sep=' ',
                ncbi_authority=False)

            hack_sol_genus.append(name['genus'])
            hack_sol_species.append(name['species'])

    ###########################################################################

    print('\n\tPreparing output directory "', output_dir, '"', sep='')
    krio.prepare_directory(output_dir)
    print('\tPreparing temporary directory "', temp_dir, '"', sep='')
    krio.prepare_directory(temp_dir)

    file_list = krio.parse_directory(search_results_dir, file_name_sep)

    unresolvable_taxonomy_dict = dict()
    for ud in unresolvable_taxonomy_list:
        unresolvable_taxonomy_dict[ud['Name']] = ud['Danger']

    # Logging -----------------------------------------------------------------

    def tax_log(accession, taxid, input_name_dict, output_name_dict, handle):
        i = input_name_dict
        o = output_name_dict

        org_name = i['genus'] + ' ' + i['species']
        danger = 0
        if org_name in unresolvable_taxonomy_dict.keys():
            danger = unresolvable_taxonomy_dict[org_name]

        handle.write(
            str(accession) + ',' +
            str(taxid) + ',' +
            str(danger) + ',' +
            i['genus'] + ',' +
            i['species'] + ',' +
            i['authority'] + ',' +
            i['subspecies'] + ',' +
            i['variety'] + ',' +
            i['cross'] + ',' +
            i['form'] + ',' +
            o['genus'] + ',' +
            o['species'] + ',' +
            o['authority'] + ',' +
            o['subspecies'] + ',' +
            o['variety'] + ',' +
            o['status'] + ',' +
            o['id'] + ',' +
            'http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=' + str(taxid) + ',' +
            'http://www.ncbi.nlm.nih.gov/nuccore/' + str(accession) + '\n')

    def tax_log_html(accession, taxid, input_name_dict, output_name_dict, handle):
        i = input_name_dict
        o = output_name_dict

        org_name = i['genus'] + ' ' + i['species']
        danger = 0
        if org_name in unresolvable_taxonomy_dict.keys():
            danger = int(unresolvable_taxonomy_dict[org_name])
        danger_td = '<td>'

        if danger == 1:
            danger_td = '<td bgcolor="#FFFF66">'
        elif danger == 2:
            danger_td = '<td bgcolor="#FFCC33">'
        elif danger == 3:
            danger_td = '<td bgcolor="#FF6600">'

        oid = '?'
        if o['status'] == 'TGRC':
            oid = o['id'].split('tgrc-')[1]
            oid = '<a href="http://tgrc.ucdavis.edu/Data/Acc/AccDetail.aspx?AccessionNum=' + oid + '">' + oid + '</a>'
        else:
            oid = o['id']

        tr = '<tr onMouseOver="this.bgColor=\'gold\';" onMouseOut="this.bgColor=\'#FFFFFF\';">'
        if o['genus'] == '':
            tr = '<tr bgcolor="#FFFF33">'

        handle.write(
            tr +
            '<td>' + '<a href="http://www.ncbi.nlm.nih.gov/nuccore/' + str(accession) + '">' + str(accession) + '</a></td>' +
            '<td>' + '<a href="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=' + str(taxid) + '">' + str(taxid) + '</a></td>' +
            danger_td + str(danger) + '</td>' +
            '<td>' + i['genus'] + '</td>' +
            '<td>' + i['species'] + '</td>' +
            '<td>' + i['authority'] + '</td>' +
            '<td>' + i['subspecies'] + '</td>' +
            '<td>' + i['variety'] + '</td>' +
            '<td>' + i['cross'] + '</td>' +
            '<td>' + i['form'] + '</td>' +
            '<td>' + o['genus'] + '</td>' +
            '<td>' + o['species'] + '</td>' +
            '<td>' + o['authority'] + '</td>' +
            '<td>' + o['subspecies'] + '</td>' +
            '<td>' + o['variety'] + '</td>' +
            '<td>' + o['status'] + '</td>' +
            '<td>' + oid + '</td>' +
            '</tr>\n')

    tax_log_file = log_dir + ps + '03-extracted-0-' + 'synonymy' + '.csv'
    tax_log_handle = open(tax_log_file, 'w')
    tax_log_handle.write('ncbiaccession,ncbitaxid,danger,igenus,ispecies,iauthority,isubspecies,ivariety,icross,iform,ogenus,ospecies,oauthority,osubspecies,ovariety,ostatus,oid,ncbitaxlink,ncbiacclink\n')
    tax_log_html_file = log_dir + ps + '03-extracted-0-' + 'synonymy' + '.html'
    tax_log_html_handle = open(tax_log_html_file, 'w')
    tax_log_html_handle.write('<!DOCTYPE html>\n<html>\n<style type="text/css">body{font-family:monospace;} table { border-collapse:collapse; } table td, table tr { border:1px solid #666;padding:3px; } </style>\n<body>\n<table>\n')
    tax_log_html_handle.write(
        '<tr bgcolor="#99FF66"><td>ncbiaccession</td><td>ncbitaxid</td><td>danger</td><td>igenus</td><td>ispecies</td><td>iauthority</td><td>isubspecies</td><td>ivariety</td><td>icross</td><td>iform</td><td>ogenus</td><td>ospecies</td><td>oauthority</td><td>osubspecies</td><td>ovariety</td><td>ostatus</td><td>oid</td></tr>\n')

    # End Logging -------------------------------------------------------------

    # Iterate over search results
    for f in file_list:

        if not f['ext'].startswith('gb'):
            continue

        # Stores all loci that have passed quality control
        sequences = list()

        # Read search results
        records = krbioio.read_sequence_file(f['path'], 'genbank')

        # Get information from the search result file name
        file_name = f['name']
        name1 = f['split'][0]
        name2 = f['split'][1]

        locus = ''
        minlen = 0
        feature_type = ''
        qualifier_label = ''
        match_stringency = ''

        for query_dict in queries:
            if name1 == query_dict['name1'] and name2 == query_dict['name2']:
                force_rev_comp = query_dict['force_rev_comp']
                locus = query_dict['locus']
                minlen = int(query_dict['minlen'])
                feature_type = query_dict['ncbi_feature_type']
                qualifier_label = query_dict['ncbi_qualifier_label']
                match_stringency = query_dict['match_stringency']
                break

        log_file = log_dir + ps + '03-extracted-' + file_name + '.tsv'
        log_handle = open(log_file, 'w')

        output_file = output_dir + ps + file_name + '.fasta'

        records_count = len(records)

        print('\n\tProcessing: ', name1, ' ', name2, ' / ', locus, sep='')

        krcl.hide_cursor()

        # Locus may have different names
        locus = locus.split(',')

        for i, record in enumerate(records):
            krcl.print_progress(i + 1, records_count, 50, '\t')

            # genbank records contain "features" which contain annotation
            #   information. Here we look for the feature that contains our
            #   target locus and note its index

            feature_indexes = list()

            loose_matching = False

            if match_stringency.lower() == 'loose':
                loose_matching = True

            for l in locus:
                fi = krseq.get_features_with_qualifier(
                    record=record, qualifier_label=qualifier_label,
                    qualifier=l.strip(), feature_type=feature_type,
                    loose=loose_matching)
                feature_indexes = feature_indexes + fi
            feature_indexes = list(set(feature_indexes))

            feature = None

            if len(feature_indexes) == 0:
                log_handle.write(
                    name1 + '_' + name2 + '\t' +
                    record.id + '\t' +
                    '\t' +
                    '\t' +
                    'No locus annotation.\n')
                continue

            elif len(feature_indexes) > 1:
                # log_handle.write(
                #     name1 + '_' + name2 + '\t' +
                #     record.id + '\t' +
                #     '\t' +
                #     '\t' +
                #     'More than one locus annotation.\n')
                # Let user pick which of the indexes to use
                print("\n\n\tFound more than one annotation for "+str(locus)+" please pick the correct index:")
                for fi in feature_indexes:
                    print("\n\t\tIndex: "+str(fi))
                    print("\t\t"+str(record.id))
                    f = record.features[fi]
                    print("\t\tStrand: "+str(f.strand))
                    print("\t\t"+str(f.location))
                    for fq in f.qualifiers.keys():
                        print("\t\t\t"+fq + ': ' + str(f.qualifiers[fq]))
                    print('\t\t---- ---- ---- ---- ---- ---- ---- ---- ---- ----')

                # picked_fi = raw_input("\n\tPick index: ")
                # if int(picked_fi) in feature_indexes:
                #     feature = record.features[int(picked_fi)]
                # else:
                #     print("\n\t Bad choice.")

                picked_index = False

                class BadChoiceException(Exception):
                    pass

                while True:
                    picked_fi = raw_input("\n\tPick index or type 'exclude' to not use this sequence: ")
                    try:
                        if str(picked_fi).lower().startswith('exclude'):
                            picked_index = False
                            print()
                            break
                        else:
                            try:
                                int(picked_fi)
                            except ValueError:
                                print("\n\tBad choice.")
                                continue
                            if int(picked_fi) in feature_indexes:
                                feature = record.features[int(picked_fi)]
                                picked_index = True
                                print()
                                break
                            else:
                                raise(BadChoiceException)
                    except BadChoiceException:
                        print("\n\tBad choice.")

                if not picked_index:
                    continue

            else:
                feature = record.features[feature_indexes[0]]

            # There should be only one matching index.
            start = int(feature.location.start)
            end = int(feature.location.end)
            strand = int(feature.location.strand)
            # Extract relevant region
            seq = record.seq[start:end]
            # If the feature is in reverse orientation, reverse-complement.
            if strand == -1 or force_rev_comp == 'yes':
                seq = seq.reverse_complement()
            # If the record is in the force_reverse_complement_list, we reverse
            # complement it, no matter what happend to it before
            if record.id in force_reverse_complement_list:
                seq = seq.reverse_complement()

            # Deal with the organism name
            tax_id = krncbi.get_ncbi_tax_id(record)

            organism = krseq.get_annotation(record, 'organism')
            organism = organism.replace(' ', '_')

            acc_name = None
            acc_name_flat = None
            taxid_name_list = None
            tried_name = None

            # HACKS ###########################################################
            hack_species_found = False
            if hacks and 'tgrc' in hacks:

                tried_name = krbionames.parse_organism_name(
                    organism,
                    sep='_',
                    ncbi_authority=True)

                if (tried_name['genus'] in hack_sol_genus and
                        tried_name['species'] in hack_sol_species):

                    hack_fi_specimen_voucher = krseq.get_features_with_qualifier_label(
                        record=record,
                        qualifier_label='specimen_voucher',
                        feature_type='source')

                    hack_fi_cultivar = krseq.get_features_with_qualifier_label(
                        record=record,
                        qualifier_label='cultivar',
                        feature_type='source')

                    hack_fi = None
                    if hack_fi_specimen_voucher:
                        hack_fi = hack_fi_specimen_voucher
                    elif hack_fi_cultivar:
                        hack_fi = hack_fi_cultivar

                    if(hack_fi):
                        hack_f = record.features[hack_fi[0]]
                        voucher = None
                        if hack_fi_specimen_voucher:
                            voucher = hack_f.qualifiers['specimen_voucher'][0]
                        elif hack_fi_cultivar:
                            voucher = hack_f.qualifiers['cultivar'][0]
                        if voucher not in do_not_repeat:
                            if voucher in found_previously.keys():
                                hack_species_found = found_previously[voucher]
                            else:
                                hack_species_found = hack_tgrc(
                                    voucher)
                            if hack_species_found:
                                found_previously[voucher] = hack_species_found
                                acc_name = hack_species_found[1]
                                acc_name['id'] = (
                                    'tgrc-' + hack_species_found[2])
                                acc_name['status'] = 'TGRC'
                                tax_log(record.id, tax_id, tried_name, acc_name, tax_log_handle)
                                tax_log_html(record.id, tax_id, tried_name, acc_name, tax_log_html_handle)
                            else:
                                do_not_repeat.append(voucher)
            ###################################################################

            if (tax_id not in keeplist_taxonomy_list) and (not hack_species_found and synonymy_table and auth_file):
                # Resolve name
                resolved = krbionames.resolve_taxid(
                    tax_id,
                    ncbi_names_table,
                    synonymy_table,
                    auth_file,
                    sorting='authority')

                acc_name = resolved[0]
                tried_name = resolved[1]

                tax_log(record.id, tax_id, tried_name, acc_name, tax_log_handle)
                tax_log_html(record.id, tax_id, tried_name, acc_name, tax_log_html_handle)

            elif (tax_id in keeplist_taxonomy_list) or (not hack_species_found):

                taxid_name_list = krbionames.names_for_ncbi_taxid(
                    tax_id, ncbi_names_table, sorting='class')

                status = 'NA'
                if tax_id in keeplist_taxonomy_list:
                    status = 'keeplist'

                acc_name = {
                    'status': status,
                    'variety': taxid_name_list[0]['variety'],
                    'authority': taxid_name_list[0]['authority'],
                    'id': 'ncbi-' + str(tax_id),
                    'subspecies': taxid_name_list[0]['subspecies'],
                    'genus': taxid_name_list[0]['genus'],
                    'species': taxid_name_list[0]['species'],
                    'form': taxid_name_list[0]['form'],
                    'cross': taxid_name_list[0]['cross']}

                tried_name = acc_name

                tax_log(record.id, tax_id, tried_name, acc_name, tax_log_handle)
                tax_log_html(record.id, tax_id, tried_name, acc_name, tax_log_html_handle)

            tried_name_flat = krbionames.flatten_organism_name(tried_name, '_')
            acc_name_flat = krbionames.flatten_organism_name(acc_name, '_')

            if acc_name['status'] == '':
                log_handle.write(
                    name1 + '_' + name2 + '\t' +
                    record.id + '\t' +
                    tried_name_flat.replace('_', ' ') + '\t' +
                    tax_id + '\t' +
                    'NO_MATCH\n')
            elif len(seq) <= minlen:
                log_handle.write(
                    name1 + '_' + name2 + '\t' +
                    record.id + '\t' +
                    tried_name_flat.replace('_', ' ') + '\t' +
                    tax_id + '\t' +
                    'SHORT' + '\t' +
                    acc_name_flat.replace('_', ' ') + '\n')
            else:
                if acc_name['status'] == 'TGRC':
                    log_handle.write(
                        name1 + '_' + name2 + '\t' +
                        record.id + '\t' +
                        tried_name_flat.replace('_', ' ') + '\t' +
                        tax_id + '\t' +
                        'TGRC' + '\t' +
                        acc_name_flat.replace('_', ' ') + '\n')
                else:
                    log_handle.write(
                        name1 + '_' + name2 + '\t' +
                        record.id + '\t' +
                        tried_name_flat.replace('_', ' ') + '\t' +
                        tax_id + '\t' +
                        'MATCH' + '\t' +
                        acc_name_flat.replace('_', ' ') + '\n')

                # Record id for the fasta output
                sequence_record_id = produce_seq_id(
                    name1, record.id, tax_id, tried_name_flat, acc_name_flat)

                # Produce a biopython sequence record object
                sequence_record = SeqRecord.SeqRecord(
                    seq=seq, id=sequence_record_id, name='', description='')

                sequences.append(sequence_record)

        log_handle.close()

        # Write results
        krbioio.write_sequence_file(sequences, output_file, 'fasta')

        print('\n\tAccepted', len(sequences), 'sequences.')

        krcl.show_cursor()

    tax_log_handle.close()
    tax_log_html_handle.write('</table>\n</body>\n</html>')
    tax_log_html_handle.close()

    all_logs = set()
    file_list = krio.parse_directory(log_dir, file_name_sep)
    for f in file_list:
        if f['ext'].startswith('tsv') and f['name'].startswith('03-extracted-'):
            handle = open(f['path'])
            for l in handle:
                all_logs.add(l)
            handle.close()
    all_logs = list(all_logs)
    all_logs_file = log_dir + ps + '03-extracted-1-all.tsv'
    handle = open(all_logs_file, 'w')
    handle.writelines(all_logs)
    handle.close()

    shutil.rmtree(temp_dir)


def one_locus_per_organism(
    extracted_results_dir,
    output_dir,
    queries,
    cutlist_records,
    cutlist_records_auto,
    # keeplist_records,
    temp_dir,
    file_name_sep,
    usearch_exe,
    aln_program_exe,
    aln_program,
    aln_options,
    good_sequences_dir,
    log_dir,
    additional_sequences
):

    '''
    Produce single sequence per locus, per organism.
    '''

    print('\nProduce single sequence per locus, per organism.')

    import os
    import shutil

    import numpy

    from Bio import SeqRecord
    from Bio import AlignIO
    from Bio import SeqIO
    from Bio import Seq
    from Bio.Align import MultipleSeqAlignment

    import krio
    import krbioio
    import kralign

    ps = os.path.sep
    ls = '\t'
    nl = '\n'

    good_sequences_dir_base = good_sequences_dir.rstrip(ps) + ps

    review_0_dir_base = output_dir + '-review-0' + ps
    review_1_dir_base = output_dir + '-review-1' + ps
    review_0_new_dir_base = output_dir + '-review-0-updated' + ps
    review_1_new_dir_base = output_dir + '-review-1-updated' + ps
    reviewed_dir_base = output_dir + '-reviewed' + ps

    lengths_log_file = log_dir + ps + '04-flat-0-lengths.csv'
    lengths_log_handle = open(lengths_log_file, 'w')
    lengths_log_handle.write('locus,mean,median,stdev\n')

    # Check if the output directory exists. If it does we will treat this as a
    # second run.
    first_run = not os.path.exists(reviewed_dir_base)
    if not first_run:
        cont = raw_input("\n\tThis is a second run, continue? Y/N: ")
        if not cont.lower() == 'y':
            exit(0)

    krio.prepare_directory(review_0_dir_base)
    krio.prepare_directory(review_1_dir_base)
    krio.prepare_directory(review_0_new_dir_base)
    krio.prepare_directory(review_1_new_dir_base)
    krio.prepare_directory(reviewed_dir_base)

    print('\n\tPreparing output directory "', output_dir, '"', sep='')
    krio.prepare_directory(output_dir)

    print('\tPreparing temporary directory "', temp_dir, '"', sep='')
    krio.prepare_directory(temp_dir)

    file_list = krio.parse_directory(extracted_results_dir, file_name_sep)

    all_records = dict()
    all_loci_dict = dict()

    for f in file_list:

        if not f['ext'].startswith('fasta'):
            continue

        name1 = f['split'][0]
        name2 = f['split'][1]

        locus_relative_position = 'x'
        for query_dict in queries:
            if name1 == query_dict['name1'] and name2 == query_dict['name2']:
                locus_relative_position = query_dict['locus_relative_position']
                break

        if not name1 in all_loci_dict:
            all_loci_dict[name1] = list()

        # separator_lrps = '%'
        # separator_lrp = '_'
        # separator_unique = '@'
        # separator_accession = '$'

        records = krbioio.read_sequence_file(f['path'], 'fasta')
        for r in records:
            parsed_id = parse_seq_id(r.id)
            r.id = produce_seq_id(
                locus=name1,
                reference_id=locus_relative_position+separator_lrps+parsed_id['reference_id'],
                taxid=parsed_id['taxid'],
                old_name=parsed_id['old_name'],
                new_name=parsed_id['new_name'])

        for r in additional_sequences:
            parsed_id = parse_seq_id(r.id)
            if parsed_id['locus'] == name1:
                records.append(r)

        all_loci_dict[name1] = all_loci_dict[name1] + records

        for r in records:
            all_records[parse_seq_id(r.id)['accessions'][0]] = r

    all_record_ids = all_records.keys()

    # LOOP ALL NAME1 ----------------------------------------------------------

    for name1 in all_loci_dict.keys():
        print('\n\tProcessing', name1)

        review_0_dir = review_0_dir_base + name1 + '_'
        review_1_dir = review_1_dir_base + name1 + '_'
        review_0_new_dir = review_0_new_dir_base + name1 + '_'
        review_1_new_dir = review_1_new_dir_base + name1 + '_'
        reviewed_dir = reviewed_dir_base + name1 + '_'
        good_sequences_dir = good_sequences_dir_base + name1 + '_'

        log_file = log_dir + ps + '04-flat-' + name1 + '.tsv'
        log_handle = open(log_file, 'w')

        name1_results = list()
        name1_results_file = output_dir + ps + name1 + '.fasta'

        name1_records = all_loci_dict[name1]

        # Keys will be new_name and value will be the list of all sequences
        # for name1 with that new_name. As new_name we will use synonymy-resolved
        # organism name
        name1_records_by_new_name = dict()
        for record in name1_records:
            parsed_id = parse_seq_id(record.id)
            new_name = parsed_id['new_name']
            if not new_name in name1_records_by_new_name:
                name1_records_by_new_name[new_name] = list()
            name1_records_by_new_name[new_name].append(record)

        # LOOP ALL new_name ------------------------------------------------------

        for new_name in name1_records_by_new_name.keys():
            # Get all the records for this new_name
            name1_records_for_new_name = list()
            for record in name1_records_by_new_name[new_name]:
                accession = parse_seq_id(record.id)['accessions'][0]
                # Filter records that appear in cutlists
                if (accession in cutlist_records) or (accession in cutlist_records_auto):
                    if accession in all_record_ids:
                        all_record_ids.remove(accession)
                else:
                    record.name = ''
                    record.description = ''
                    name1_records_for_new_name.append(record)
            name1_records_for_new_name.sort(reverse=True)

            # Deduplicate records, if there is more than one record with the
            # same accession, we keep them only if the sequences differ
            dedupe_records_dict = dict()
            for i, r in enumerate(name1_records_for_new_name):

                id_match = False
                lrp_match = False
                seq_longer = False

                r_accession = parse_seq_id(r.id)['accessions'][0]
                r_lrp = parse_seq_id(r.id)['lrps'][0]

                for k in dedupe_records_dict.keys():

                    dedupe_accession = parse_seq_id(k)['accessions'][0]
                    dedupe_lrp = parse_seq_id(k)['lrps'][0]

                    if r_accession == dedupe_accession:
                        id_match = True
                        if str(r.seq).lower() == str(dedupe_records_dict[k].seq).lower():
                            if dedupe_lrp.lower() == 'x':
                                if r_lrp.lower() != 'x':
                                    dedupe_records_dict[k] = r
                                    lrp_match = True
                        if not lrp_match and r_lrp == dedupe_lrp:
                            lrp_match = True
                            if len(r.seq) > len(dedupe_records_dict[k].seq):
                                seq_longer = True

                if id_match:
                    if lrp_match:
                        if seq_longer:
                            dedupe_records_dict[r.id] = r
                    else:
                        dedupe_records_dict[r.id] = r
                else:
                    dedupe_records_dict[r.id] = r

            # Check if there are previous alignments or sequences for this
            # new_name
            new_records = False

            good_sequences_seq = os.path.exists(good_sequences_dir + new_name + '.fasta')
            good_sequences_aln = os.path.exists(good_sequences_dir + new_name + '.phy')

            if good_sequences_seq and good_sequences_aln:
                print('\n\tWARNING: Both the alignment and sequence file exists for ' + name1 + ' ' + new_name + '.')

            if first_run and (good_sequences_seq or good_sequences_aln):
                all_previous_ids = set()
                if good_sequences_seq:
                    seq = SeqIO.read(good_sequences_dir + new_name + '.fasta', "fasta")
                    all_previous_ids.add(parse_seq_id(seq.id)['accessions'][0])
                elif good_sequences_aln:
                    aln = AlignIO.read(good_sequences_dir + new_name + '.phy', "phylip-relaxed")
                    for s in aln:
                        all_previous_ids.add(parse_seq_id(s.id)['accessions'][0])
                all_ids_in_dedupe = set()
                for k in dedupe_records_dict.keys():
                    all_ids_in_dedupe.add(parse_seq_id(k)['accessions'][0])
                ids_union = all_previous_ids.union(all_ids_in_dedupe)
                if len(ids_union) > len(all_previous_ids):
                    new_records = True

            # Combine sequences by locus relative position if from the same accession
            dedupe_records_lrp = dedupe_records_dict.values()
            if len(dedupe_records_lrp) > 1:
                dedupe_records_by_accession = dict()
                dedupe_records_lrp = list()
                for r in dedupe_records_dict.values():
                    accession = parse_seq_id(r.id)['accessions'][0]
                    if not accession in dedupe_records_by_accession:
                        dedupe_records_by_accession[accession] = list()
                    dedupe_records_by_accession[accession].append(r)

                for acc_list in dedupe_records_by_accession.values():
                    if len(acc_list) > 1:

                        parsed_id = parse_seq_id(acc_list[0].id)
                        accession = parsed_id['accessions'][0]

                        acc_list.sort(key=lambda x: x.id, reverse=False)

                        cat_seq = ''

                        lrps = list()

                        for r in acc_list:
                            lrp = parse_seq_id(r.id)['lrps'][0]
                            if lrp.lower().startswith('x'):
                                dedupe_records_lrp.append(r)
                                continue
                            cat_seq = cat_seq + str(r.seq)
                            lrps.append(lrp)

                        if len(lrps) > 0:
                            cat_id = produce_seq_id(
                                locus=name1,
                                reference_id=separator_lrp.join(lrps) + separator_lrps + accession,
                                taxid=parsed_id['taxid'],
                                old_name=parsed_id['old_name'],
                                new_name=parsed_id['new_name'])

                            sequence_record = SeqRecord.SeqRecord(
                                seq=Seq.Seq(cat_seq), id=cat_id, name='',
                                description='')
                            dedupe_records_lrp.append(sequence_record)
                    else:
                        dedupe_records_lrp.append(acc_list[0])

            dir_0 = None
            dir_1 = None

            if first_run and new_records:
                dir_0 = review_0_new_dir
                dir_1 = review_1_new_dir
            else:
                dir_0 = review_0_dir
                dir_1 = review_1_dir

            if first_run and len(dedupe_records_lrp) > 1:
                if new_records or not good_sequences_aln:
                    aln = kralign.align(
                        records=dedupe_records_lrp,
                        program=aln_program,
                        options=aln_options,
                        program_executable=aln_program_exe)
                    aln.sort()
                    cons = kralign.consensus(
                        aln, threshold=0.1, unknown='N',
                        resolve_ambiguities=False)
                    bases_at_sites = cons[1]
                    tot_bs = 0
                    for bs in bases_at_sites:
                        tot_bs = tot_bs + len(bs)
                    similarity = float(len(bases_at_sites)) / float(tot_bs)
                    if similarity >= 0.95:
                        AlignIO.write(aln, dir_1 + new_name + '.phy', "phylip-relaxed")
                    else:
                        #######################################################
                        n_count = 10
                        new_aln = None
                        all_lrps = set()
                        for s in aln:
                            seq_id = s.id
                            parsed_id = parse_seq_id(seq_id)
                            lrps = parsed_id['lrps']
                            all_lrps |= set(lrps)
                        if len(all_lrps) == 1 and list(all_lrps)[0].lower() == 'x':
                            new_aln = aln
                        else:
                            seq_records = list()
                            last_lrp = None
                            prev_seq_len = 0
                            pos = 0
                            for s in aln:
                                seq = str(s.seq).replace('-', '')
                                seq_len = len(seq)
                                seq_id = s.id
                                parsed_id = parse_seq_id(seq_id)
                                lrps = parsed_id['lrps']
                                if not last_lrp:
                                    last_lrp = lrps[-1]
                                else:
                                    if lrps[-1].lower() == 'x':
                                        pass
                                    elif last_lrp < lrps[-1]:
                                        last_lrp = lrps[-1]
                                        pos = pos + prev_seq_len
                                        seq = pos * '-' + n_count * 'N' + seq
                                        pos = pos + n_count
                                        prev_seq_len = seq_len + n_count
                                    elif last_lrp == lrps[-1]:
                                        seq = pos * '-' + seq
                                prev_seq_len = max(prev_seq_len, seq_len)
                                sequence_record = SeqRecord.SeqRecord(
                                    seq=Seq.Seq(seq), id=seq_id, name='',
                                    description='')
                                seq_records.append(sequence_record)
                            max_len = 0
                            for s in seq_records:
                                max_len = max(max_len, len(s.seq))
                            for s in seq_records:
                                len_diff = max_len - len(s.seq)
                                seq = str(s.seq) + len_diff * '-'
                                s.seq = Seq.Seq(seq)
                            new_aln = MultipleSeqAlignment(seq_records)
                        #######################################################
                        AlignIO.write(new_aln, dir_0 + new_name + '.phy', "phylip-relaxed")
                elif good_sequences_aln:
                    aln = AlignIO.read(good_sequences_dir + new_name + '.phy', "phylip-relaxed")
                    AlignIO.write(aln, reviewed_dir + new_name + '.phy', "phylip-relaxed")
            elif first_run and len(dedupe_records_lrp) == 1:
                if new_records or not good_sequences_seq:
                    SeqIO.write(dedupe_records_lrp[0], dir_1 + new_name + '.fasta', "fasta")
                elif good_sequences_seq:
                    seq = SeqIO.read(good_sequences_dir + new_name + '.fasta', "fasta")
                    SeqIO.write(seq, reviewed_dir + new_name + '.fasta', "fasta")
            elif not first_run and os.path.exists(reviewed_dir + new_name + '.phy'):
                aln = AlignIO.read(reviewed_dir + new_name + '.phy', "phylip-relaxed")
                cons_ids = set()
                for s in aln:
                    parsed_id = parse_seq_id(s.id)
                    accession = parsed_id['accessions'][0]
                    if accession in all_record_ids:
                        all_record_ids.remove(accession)
                    cons_ids.add(accession)

                ###############################################################
                consensus = kralign.consensus(aln, threshold=0.4, unknown='N')
                consensus_seq = consensus[0]
                bases_per_site = consensus[3]

                if len(aln) > 3 or bases_per_site < 1.7:
                    pass
                else:
                    len_of_seq_in_aln = list()
                    for a in aln:
                        seq_str = str(a.seq).lower().replace('n', '').replace('-', '').replace('.', '')
                        len_of_seq_in_aln.append(len(seq_str))
                    max_len_index = len_of_seq_in_aln.index(max(len_of_seq_in_aln))
                    consensus_seq = str(aln[max_len_index].seq).replace('-', '').replace('.', '')
                    consensus_seq = Seq.Seq(consensus_seq)

                sequence_record = SeqRecord.SeqRecord(
                    seq=consensus_seq, id=new_name, name='', description='')
                ###############################################################

                name1_results.append(sequence_record)
                log_message = new_name + ls + ls.join(cons_ids) + nl
                log_handle.write(log_message)
            elif not first_run and os.path.exists(reviewed_dir + new_name + '.fasta'):
                seq = SeqIO.read(reviewed_dir + new_name + '.fasta', "fasta")
                parsed_id = parse_seq_id(seq.id)
                accession = parsed_id['accessions'][0]
                if accession in all_record_ids:
                    all_record_ids.remove(accession)
                log_message = new_name + ls + accession + nl
                log_handle.write(log_message)
                seq.id = new_name
                seq.name = ''
                seq.description = ''
                name1_results.append(seq)
            elif not first_run and len(dedupe_records_lrp) > 0:
                print('\n\tWARNING: Make sure to review all the alignments and sequences.')
                cont = raw_input("\tWARNING: Could not find any sequences for "+name1+" "+new_name+". Continue? Y/N: ")
                if not cont.lower().startswith('y'):
                    shutil.rmtree(temp_dir)
                    exit(0)

        # END LOOP ALL new_name --------------------------------------------------

        if not first_run:

            # Produce sequence length stats
            s_lengths = list()
            for r in name1_results:
                s_lengths.append(len(r.seq))

            s_mean = str(numpy.mean(s_lengths))
            s_median = str(numpy.median(s_lengths))
            s_stdev = str(numpy.std(s_lengths))

            lengths_log_handle.write(name1 + ',' + s_mean + ',' + s_median + ',' + s_stdev + '\n')

            krbioio.write_sequence_file(name1_results, name1_results_file, 'fasta')
            print('\tAccepted sequences from', len(name1_results), 'taxa.')

        log_handle.close()

    # END LOOP ALL NAME1 ------------------------------------------------------

    lengths_log_handle.close()

    if not first_run:
        cutlist_records_auto |= set(all_record_ids)

    shutil.rmtree(temp_dir)

    if first_run:
        exit(0)

    return(cutlist_records_auto)


def align_loci(
    processed_results_dir,
    output_dir,
    aln_program_exe,
    aln_program,
    aln_options,
    temp_dir
):

    '''
    Align individual loci.
    '''

    print('\nAlign individual loci.')

    import os
    # import shutil

    import krio
    import krbioio
    import kralign

    ps = os.path.sep

    print('\n\tPreparing output directory "', output_dir, '"', sep='')
    krio.prepare_directory(output_dir)
    # print('\tPreparing temporary directory "', temp_dir, '"', sep='')
    # krio.prepare_directory(temp_dir)
    print()

    file_list = krio.parse_directory(processed_results_dir, ' ')

    for f in file_list:

        if not f['ext'].startswith('fasta'):
            continue

        file_name = f['name']
        output_file = output_dir + ps + file_name + '.phy'
        records = krbioio.read_sequence_file(f['path'], 'fasta')

        print('\tAligning', file_name)
        aln = kralign.align(
            records=records,
            program=aln_program,
            options=aln_options,
            program_executable=aln_program_exe)
        if aln:
            krbioio.write_alignment_file(aln, output_file, 'phylip-relaxed')

    # shutil.rmtree(temp_dir)


def concatenate(
    aligned_loci_dir,
    output_dir,
    order_of_loci,
    number_of_gaps_between_loci,
    log_dir
):

    '''
    Concatenate the alignments of individual loci.

    number_of_gaps_between_loci - the number of gaps to insert between individual alignments.
    '''

    print('\nConcatenate the alignments of individual loci.')

    import os

    from Bio import AlignIO

    import krio
    import krbioio
    import kralign
    import copy

    ps = os.path.sep

    print('\n\tPreparing output directory "', output_dir, '"', sep='')
    krio.prepare_directory(output_dir)

    file_list = krio.parse_directory(aligned_loci_dir, ' ')

    order_list = [x.strip() for x in order_of_loci]
    alignments = [x.strip() for x in order_of_loci]

    for f in file_list:

        if not f['ext'].startswith('phy'):
            continue

        file_name = f['name']

        aln = AlignIO.read(f['path'], "phylip-relaxed")
        if aln:
            i = alignments.index(file_name)
            alignments[i] = (aln, file_name)

    for aln in alignments:
        if isinstance(aln, basestring):
            alignments.remove(aln)

    print('\n\tProducing concatenated alignment.')
    if alignments:

        # Produce presence/absence matrix
        presence_list = list()
        length_list = list()
        for p in range(0, len(order_list)):
            presence_list.append('0')
        matrix = dict()
        for a in alignments:
            length_list.append(str(a[0].get_alignment_length()))
            for s in a[0]:
                if not s.id in matrix:
                    matrix[s.id] = copy.copy(presence_list)
        for a in alignments:
            for s in a[0]:
                idx = order_list.index(a[1])
                matrix[s.id][idx] = '1'
        matrix_output_file = log_dir + ps + '06-locus-presence' + '.csv'
        f = open(matrix_output_file, 'wb')
        f.write('taxon' + ',' + 'count' + ',' + ','.join(order_list) + '\n')
        f.write('' + ',' + '' + ',' + ','.join(length_list) + '\n')
        for key in matrix.keys():
            f.write(key + ',' + str(matrix[key].count('1')) + ',' +
                    ','.join(matrix[key]) + '\n')
        f.close()

        # Concatenate
        raw_alignments = list()
        for a in alignments:
            raw_alignments.append(a[0])
        concatenated = kralign.concatenate(raw_alignments, int(number_of_gaps_between_loci))
        concatenated_output_file = output_dir + ps + 'concatenated' + '.phy'
        krbioio.write_alignment_file(concatenated, concatenated_output_file,
                                     'phylip-relaxed')

# End pipeline functions ------------------------------------------------------
