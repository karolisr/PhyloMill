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

def search_and_download(queries, output_dir, file_name_sep, email):

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
            # Download records.
            krncbi.download_sequence_records(file_path, result_uids, db, email)


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


def extract_loci(search_results_dir, output_dir, queries, ncbi_names_table,
                 temp_dir, file_name_sep, synonymy_table=None, auth_file=None,
                 hacks=None, hacks_data_location=None, unresolvable_taxonomy_list=None):

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
            # path='..' + ps + 'data' + ps + 'sol_sp_with_vouchers',
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

    tax_log_file = output_dir + ps + 'synonymy' + '.log'
    tax_log_handle = open(tax_log_file, 'w')
    tax_log_handle.write('ncbiaccession,ncbitaxid,danger,igenus,ispecies,iauthority,isubspecies,ivariety,icross,iform,ogenus,ospecies,oauthority,osubspecies,ovariety,ostatus,oid,ncbitaxlink,ncbiacclink\n')

    tax_log_html_file = output_dir + ps + 'synonymy' + '.html'
    tax_log_html_handle = open(tax_log_html_file, 'w')

    tax_log_html_handle.write('<!DOCTYPE html>\n<html>\n<style type="text/css">body{font-family:monospace;} table { border-collapse:collapse; } table td, table tr { border:1px solid #666;padding:3px; } </style>\n<body>\n<table>\n')

    tax_log_html_handle.write(
        '<tr bgcolor="#99FF66"><td>ncbiaccession</td><td>ncbitaxid</td><td>danger</td><td>igenus</td><td>ispecies</td><td>iauthority</td><td>isubspecies</td><td>ivariety</td><td>icross</td><td>iform</td><td>ogenus</td><td>ospecies</td><td>oauthority</td><td>osubspecies</td><td>ovariety</td><td>ostatus</td><td>oid</td></tr>\n')

    # Iterate over search results
    for f in file_list:

        if not f['ext'].startswith('gb'):
            continue

        # Stores all loci that have passed quality control
        loci = list()

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

        log_file = output_dir + ps + file_name + '.log'
        log_handle = open(log_file, 'w')

        output_file = output_dir + ps + file_name + '.fasta'

        records_count = len(records)

        print('\n\tProcessing: ', name1, ' ', name2, ' / ', locus, sep='')

        krcl.hide_cursor()

        # Locus may have different names
        locus = locus.split(',')

        for i, record in enumerate(records):
            krcl.print_progress(i + 1, records_count, 50, '\t')

            # genbank records contain "features" which conatin annotation
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
            # TODO: this should never occur, and occured only once
            # Same gene annotated more than once
            if len(feature_indexes) > 1:
                log_handle.write(
                    name1 + '_' + name2 + '\t' +
                    record.id + '\t' +
                    '\t' +
                    '\t' +
                    'More than one locus annotation.\n')
                continue
            if len(feature_indexes) == 0:
                log_handle.write(
                    name1 + '_' + name2 + '\t' +
                    record.id + '\t' +
                    '\t' +
                    '\t' +
                    'No locus annotation.\n')
                continue
            # There should be only one matching index.
            feature = record.features[feature_indexes[0]]
            start = int(feature.location.start)
            end = int(feature.location.end)
            strand = int(feature.location.strand)
            # Extract relevant region
            seq = record.seq[start:end]
            # If the feature is in reverse orientation, reverse-complement.
            if strand == -1 or force_rev_comp == 'yes':
                seq = seq.reverse_complement()

            # Deal with the organism name
            tax_id = krncbi.get_ncbi_tax_id(record)

            organism = krseq.get_annotation(record, 'organism')
            organism = organism.replace(' ', '_')

            acc_name = None
            acc_name_flat = None
            organism_authority = None

            # HACKS ###########################################################
            hack_species_found = False
            if hacks and 'tgrc' in hacks:

                hack_org = krbionames.parse_organism_name(
                    organism,
                    sep='_',
                    ncbi_authority=True)

                if (hack_org['genus'] in hack_sol_genus and
                        hack_org['species'] in hack_sol_species):

                    hack_fi = krseq.get_features_with_qualifier_label(
                        record=record,
                        qualifier_label='specimen_voucher',
                        feature_type='source')

                    if(hack_fi):
                        hack_f = record.features[hack_fi[0]]
                        voucher = hack_f.qualifiers['specimen_voucher'][0]
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
                                acc_name_flat = (
                                    krbionames.flatten_organism_name(
                                        acc_name, '_'))
                                tax_log(record.id, tax_id, hack_org, acc_name, tax_log_handle)
                                tax_log_html(record.id, tax_id, hack_org, acc_name, tax_log_html_handle)
                            else:
                                do_not_repeat.append(voucher)
            ###################################################################

            loose_mode = False

            if not hack_species_found and synonymy_table and auth_file:
                # A list of organism names based on NCBI taxid. This is a
                #   sorted list with the most complete names at lower indexes.
                organism_authority = krbionames.names_for_ncbi_taxid(
                    tax_id, ncbi_names_table, sorting='authority')

                # Iterate over the list of NCBI names and try to resolve
                #   accepted name.

                # First we look for matches using STRICT mode
                # First look if there is a best possible match "acc"
                found_match = False
                for oa in organism_authority:
                    acc_name = krbionames.accepted_name(
                        name=oa,
                        synonymy_table=synonymy_table,
                        auth_file=auth_file,
                        allow_loose_matching=False)
                    if acc_name and acc_name['status'].lower() == 'acc':
                        found_match = True
                        tax_log(record.id, tax_id, oa, acc_name, tax_log_handle)
                        tax_log_html(record.id, tax_id, oa, acc_name, tax_log_html_handle)
                        break

                # If no, then look for the next best thing "prov"
                if not found_match:
                    for oa in organism_authority:
                        acc_name = krbionames.accepted_name(
                            name=oa,
                            synonymy_table=synonymy_table,
                            auth_file=auth_file,
                            allow_loose_matching=False)
                        if acc_name and acc_name['status'].lower() == 'prov':
                            found_match = True
                            tax_log(record.id, tax_id, oa, acc_name, tax_log_handle)
                            tax_log_html(record.id, tax_id, oa, acc_name, tax_log_html_handle)
                            break

                # Otherwise, let's find something that isn't blank
                if not found_match:
                    for oa in organism_authority:
                        acc_name = krbionames.accepted_name(
                            name=oa,
                            synonymy_table=synonymy_table,
                            auth_file=auth_file,
                            allow_loose_matching=False)
                        if (acc_name and (
                            acc_name['status'].lower() == 'as' or
                            acc_name['status'].lower() == 'nn' or
                            acc_name['status'].lower() == 'unc' or
                                acc_name['status'].lower() == 'unr')):
                            found_match = True
                            tax_log(record.id, tax_id, oa, acc_name, tax_log_handle)
                            tax_log_html(record.id, tax_id, oa, acc_name, tax_log_html_handle)
                            break

                # If we find nothing using STRICT mode we look for matches
                # using LOOSE mode
                # First look if there is a best possible match "acc"
                if not found_match:
                    for oa in organism_authority:
                        acc_name = krbionames.accepted_name(
                            name=oa,
                            synonymy_table=synonymy_table,
                            auth_file=auth_file,
                            allow_loose_matching=True)
                        if acc_name and acc_name['status'].lower().startswith('acc'):
                            loose_mode = True
                            found_match = True
                            tax_log(record.id, tax_id, oa, acc_name, tax_log_handle)
                            tax_log_html(record.id, tax_id, oa, acc_name, tax_log_html_handle)
                            break

                # If no, then look for the next best thing "prov"
                if not found_match:
                    for oa in organism_authority:
                        acc_name = krbionames.accepted_name(
                            name=oa,
                            synonymy_table=synonymy_table,
                            auth_file=auth_file,
                            allow_loose_matching=True)
                        if acc_name and acc_name['status'].lower().startswith('prov'):
                            loose_mode = True
                            found_match = True
                            tax_log(record.id, tax_id, oa, acc_name, tax_log_handle)
                            tax_log_html(record.id, tax_id, oa, acc_name, tax_log_html_handle)
                            break

                # Otherwise, let's find something that isn't blank
                if not found_match:
                    for oa in organism_authority:
                        acc_name = krbionames.accepted_name(
                            name=oa,
                            synonymy_table=synonymy_table,
                            auth_file=auth_file,
                            allow_loose_matching=True)
                        if (acc_name and (
                            acc_name['status'].lower().startswith('as') or
                            acc_name['status'].lower().startswith('nn') or
                            acc_name['status'].lower().startswith('unc') or
                                acc_name['status'].lower().startswith('unr'))):
                            loose_mode = True
                            found_match = True
                            tax_log(record.id, tax_id, oa, acc_name, tax_log_handle)
                            tax_log_html(record.id, tax_id, oa, acc_name, tax_log_html_handle)
                            break

                if not found_match:
                    log_handle.write(
                        name1 + '_' + name2 + '\t' +
                        record.id + '\t' +
                        organism.replace('_', ' ') + '\t' +
                        tax_id + '\t'
                        'No taxonomic match.\n')
                    for oa in organism_authority:
                        tax_log(
                            record.id,
                            tax_id,
                            oa,
                            # {'form': '', 'variety': '', 'cross': '', 'authority': '', 'other': '', 'name_class': '', 'subspecies': '', 'genus': '', 'species': ''},
                            {'status': '', 'variety': '', 'authority': '', 'id': '', 'subspecies': '', 'genus': '', 'species': ''},
                            tax_log_handle)

                        tax_log_html(
                            record.id,
                            tax_id,
                            oa,
                            # {'form': '', 'variety': '', 'cross': '', 'authority': '', 'other': '', 'name_class': '', 'subspecies': '', 'genus': '', 'species': ''},
                            {'status': '', 'variety': '', 'authority': '', 'id': '', 'subspecies': '', 'genus': '', 'species': ''},
                            tax_log_html_handle)
                    continue

                acc_name_flat = krbionames.flatten_organism_name(acc_name, '_')
                if loose_mode:
                    pass

            elif not hack_species_found:
                acc_name = dict()
                acc_name['id'] = tax_id
                acc_name['status'] = 'NA'
                organism_authority = krbionames.names_for_ncbi_taxid(
                    tax_id, ncbi_names_table, sorting='class')
                acc_name_flat = krbionames.flatten_organism_name(
                    organism_authority[0], '_')

            # Record id for the fasta output
            sequence_record_id = (
                # NCBI accession
                record.id + '|' +
                # Organism name as it appears in search results
                organism + '|' +
                # NCBI taxid
                tax_id + '|' +
                # Accepted name
                acc_name_flat + '|' +
                # Accepted name id from the synonymy table
                acc_name['id'])
            # Produce a biopython sequence record object
            sequence_record = SeqRecord.SeqRecord(
                seq=seq, id=sequence_record_id, name='', description='')

            if acc_name['status'] == '':
                log_handle.write(
                    name1 + '_' + name2 + '\t' +
                    record.id + '\t' +
                    organism.replace('_', ' ') + '\t' +
                    tax_id + '\t' +
                    'No taxonomic match.\n')
            elif len(seq) <= minlen:
                log_handle.write(
                    name1 + '_' + name2 + '\t' +
                    record.id + '\t' +
                    organism.replace('_', ' ') + '\t' +
                    tax_id + '\t' +
                    'Sequence is too short.\n')

            else:
                if loose_mode:
                    log_handle.write(
                        name1 + '_' + name2 + '\t' +
                        record.id + '\t' +
                        organism.replace('_', ' ') + '\t' +
                        tax_id + '\t' +
                        'Taxonomic match using LOOSE mode.' + '\t' +
                        acc_name_flat.replace('_', ' ') + '\n')
                elif acc_name['status'] == 'TGRC':
                    log_handle.write(
                        name1 + '_' + name2 + '\t' +
                        record.id + '\t' +
                        organism.replace('_', ' ') + '\t' +
                        tax_id + '\t' +
                        'Match TGRC' + '\t' +
                        acc_name_flat.replace('_', ' ') + '\n')
                else:
                    log_handle.write(
                        name1 + '_' + name2 + '\t' +
                        record.id + '\t' +
                        organism.replace('_', ' ') + '\t' +
                        tax_id + '\t' +
                        'Match' + '\t' +
                        acc_name_flat.replace('_', ' ') + '\n')
                loci.append(sequence_record)

        log_handle.close()

        # Write results
        krbioio.write_sequence_file(loci, output_file, 'fasta')

        print('\n\tAccepted', len(loci), 'sequences.')

        krcl.show_cursor()

    tax_log_handle.close()
    tax_log_html_handle.write('</table>\n</body>\n</html>')
    tax_log_html_handle.close()

    all_logs = set()
    file_list = krio.parse_directory(output_dir, file_name_sep)
    for f in file_list:
        if not f['ext'].startswith('log'):
            continue
        handle = open(f['path'])
        for l in handle:
            all_logs.add(l)
        handle.close()
    all_logs = list(all_logs)
    all_logs_file = output_dir + ps + 'all.log'
    handle = open(all_logs_file, 'w')
    handle.writelines(all_logs)
    handle.close()

    shutil.rmtree(temp_dir)


def one_locus_per_organism(
    extracted_results_dir,
    output_dir,
    queries,
    identity_threshold,
    cutlist_records_file,
    keeplist_records_file,
    temp_dir,
    file_name_sep,
    usearch_exe,
    aln_program_exe,
    aln_program,
    aln_options,
    threads=1,
    proceed_without_review=False
):

    '''
    Produce single sequence per locus, per organism.
    '''

    print('\nProduce single sequence per locus, per organism.')

    import os
    import shutil

    from Bio import Seq
    from Bio import SeqRecord
    from Bio import AlignIO

    import krio
    import krbioio
    import krusearch
    import krcl
    import kralign
    import krseq

    krcl.hide_cursor()

    ps = os.path.sep

    print('\n\tPreparing output directory "', output_dir, '"', sep='')
    krio.prepare_directory(output_dir)

    print('\tPreparing temporary directory "', temp_dir, '"', sep='')
    krio.prepare_directory(temp_dir)

    file_list = krio.parse_directory(extracted_results_dir, file_name_sep)

    # Read in accessions to filter
    cutlist_records = None
    if cutlist_records_file:

        cutlist_records = krio.read_table_file(
            path=cutlist_records_file,
            has_headers=False,
            headers=None,
            delimiter=',',
            quotechar='"',
            rettype='set')

    # Read in white-listed accessions
    keeplist_records = None
    if keeplist_records_file:

        keeplist_records = krio.read_table_file(
            path=keeplist_records_file,
            has_headers=False,
            headers=None,
            delimiter=',',
            quotechar='"',
            rettype='set')

    #
    # Iterate over each file that corresponds to a single query line
    # There are as many files as there are queries. Using these we will build
    # a structure that looks like this:
    #
    #   Nomenclature:
    #       name1: general query name
    #       name2: specific query name
    #
    #   Structure:
    #       all_loci_dict: dictionary with key name1, each value is a list of
    #       dictionaries with keys:
    #           name2:      specific name (string)
    #           records:    list of BioSeqRecord objects (list)
    #           lrp:        locus relative position (int or 'X')
    #

    all_records = dict()
    all_loci_dict = dict()

    for f in file_list:

        if not f['ext'].startswith('fasta'):
            continue

        name1 = f['split'][0]
        name2 = f['split'][1]

        locus_relative_position = 0
        for query_dict in queries:
            if name1 == query_dict['name1'] and name2 == query_dict['name2']:
                locus_relative_position = query_dict['locus_relative_position']
                break

        if not name1 in all_loci_dict:
            all_loci_dict[name1] = list()

        records = krbioio.read_sequence_file(f['path'], 'fasta')
        all_loci_dict[name1].append({
            'name2': name2,
            'records': records,
            'lrp': locus_relative_position})

        for r in records:
            all_records[r.id.split('|')[0]] = r

    # Because for each locus, there could be several records from the same
    # species, we need to pick the best representative sequence. Only one
    # per locus, per organism. Each name1 will have one and only one file
    # associated with it.

    review_needed = False

    # LOOP ALL NAME1 ----------------------------------------------------------
    for name1 in all_loci_dict.keys():
        print('\n\tProcessing', name1)

        results = list()

        log_file = output_dir + ps + name1 + '.log'
        log_handle = open(log_file, 'w')
        output_file = output_dir + ps + name1 + '.fasta'

        # Each row will contain a dictionary with a list of results and a
        # relative locus position:
        name1_results = list()

        # All name2 associated with name1
        list_of_name2 = all_loci_dict[name1]

        # LOOP ALL NAME2 ------------------------------------------------------
        for name2_dict in list_of_name2:

            name2 = name2_dict['name2']
            records = name2_dict['records']
            lrp = name2_dict['lrp']

            if lrp.lower() != 'x':
                lrp = int(lrp)

            # Cluster top X percent of records within a locus within NAME2
            # Get a representative consensus sequence
            records.sort(key=lambda x: len(x), reverse=True)
            records_count = len(records)
            fraction_to_cluster = 0.01
            number_to_cluster = int(records_count * fraction_to_cluster)
            number_to_cluster = max(10, number_to_cluster)
            # TODO: If there are not enough records for a meaningful cluster,
            # they should be automatically sent for manual review.
            if number_to_cluster >= records_count:
                number_to_cluster = max(5, records_count/2)

            records_to_cluster = records[0:number_to_cluster]

            cluster_dict = krusearch.cluster_records(
                records=records_to_cluster,
                identity_threshold=identity_threshold,
                temp_dir=temp_dir,
                sorted_input=False,
                algorithm='smallmem',
                strand='both',
                threads=1,
                quiet=True,
                program=usearch_exe,
                heuristics=False,
                query_coverage=0.1,
                target_coverage=0.1,
                sizein=False,
                sizeout=False,
                usersort=False
            )

            clusters = cluster_dict.values()
            clusters.sort(key=lambda x: len(x), reverse=True)
            largest_cluster = clusters[0]
            largest_cluster_names = list()
            for n in largest_cluster:
                largest_cluster_names.append(n[1])
            records_to_align = list()
            for record in records:
                if record.id in largest_cluster_names:
                    records_to_align.append(record)
            aln = kralign.align(
                records=records_to_align,
                program=aln_program,
                options=aln_options,
                program_executable=aln_program_exe)
            # TODO: Better consensus; This consensus should not include
            # ambiguous characters because uclust does not understand them.
            # Randomly pick a representative character whenever an ambiguity
            # occurs.
            # summary_aln = AlignInfo.SummaryInfo(aln)
            # consensus = summary_aln.dumb_consensus(
            #     threshold=0.001, ambiguous='N')
            consensus = kralign.consensus(aln, threshold=0.4, unknown='N', resolve_ambiguities=True)
            rep_cons_record = SeqRecord.SeqRecord(
                seq=consensus, id='representative', name='', description='')

            # Keys will be taxid and value will be the list of all sequences
            # for name2 with that taxid. As taxid we will use synonymy-resolved
            # organism name
            raw_name2_dict = dict()
            for record in records:
                taxid = record.description.split('|')[3]
                if not taxid in raw_name2_dict:
                    raw_name2_dict[taxid] = list()
                raw_name2_dict[taxid].append(record)

            # Will hold one sequence record per taxid, keys are taxid
            flat_name2_dict = dict()
            taxa_count = len(raw_name2_dict.keys())

            print('\n\t\tProcessing', name1, name2 + '. There are ' +
                  str(records_count) + ' records from ' + str(taxa_count) +
                  ' taxa.')

            for i, taxid in enumerate(raw_name2_dict.keys()):

                krcl.print_progress(i + 1, taxa_count, 50, '\t\t')

                all_records_for_taxid_prefiltered = raw_name2_dict[taxid]
                all_records_for_taxid = list()

                # Filter unwanted accessions
                for tr in all_records_for_taxid_prefiltered:
                    if tr.id.split('|')[0] not in cutlist_records:
                        all_records_for_taxid.append(tr)

                if len(all_records_for_taxid) > 1:
                    # If the records were not white-listed we need to check
                    # if they look like our locus
                    all_records_for_taxid_to_cluster = list()
                    for tr in all_records_for_taxid:
                        if tr.id.split('|')[0] not in keeplist_records:
                            resolved_tr_seq = Seq.Seq(krseq.resolve_ambiguities(str(tr.seq)))
                            resolved_tr_seq_record = SeqRecord.SeqRecord(
                                seq=resolved_tr_seq, id=tr.id, name='',
                                description='')
                            all_records_for_taxid_to_cluster.append(resolved_tr_seq_record)
                    all_records_for_taxid_to_cluster.append(rep_cons_record)
                    # ...we cluster these sequences and hope for only one
                    # cluster
                    cluster_dict = None
                    if len(all_records_for_taxid_to_cluster) > 1:
                        cluster_dict = krusearch.cluster_records(
                            records=all_records_for_taxid_to_cluster,
                            identity_threshold=identity_threshold,
                            temp_dir=temp_dir,
                            sorted_input=False,
                            algorithm='smallmem',
                            strand='plus',
                            threads=1,
                            quiet=True,
                            program=usearch_exe,
                            heuristics=False,
                            query_coverage=0.1,
                            target_coverage=0.1,
                            sizein=False,
                            sizeout=False,
                            usersort=False
                        )

                    # ...if there is only one cluster
                    if not cluster_dict or len(cluster_dict.keys()) == 1:
                        # We align the sequences in a cluster and produce
                        # a consensus
                        aln = kralign.align(
                            records=all_records_for_taxid,
                            program=aln_program,
                            options=aln_options,
                            program_executable=aln_program_exe)
                        # TODO: Better consensus. This consensus should include
                        # correct ambiguities.
                        # summary_aln = AlignInfo.SummaryInfo(aln)
                        # consensus = summary_aln.dumb_consensus(
                        #     threshold=0.001, ambiguous='N')
                        consensus = kralign.consensus(aln, threshold=0.4, unknown='N')

                        # We want the consensus record description to reflect
                        # all the sequence ids that were combined to form it
                        cons_ids = list()
                        for n in all_records_for_taxid:
                            cons_ids.append(n.id.split('|')[0])
                        new_cons_id = '$'.join(cons_ids)

                        old_record_id_split = (
                            all_records_for_taxid[0].id.split('|'))

                        sequence_record_id = (
                            'consensus$' + new_cons_id + '|' +
                            old_record_id_split[1] + '|' +
                            old_record_id_split[2] + '|' +
                            old_record_id_split[3] + '|' +
                            old_record_id_split[4])

                        # Produce a biopython sequence record object
                        sequence_record = SeqRecord.SeqRecord(
                            seq=consensus, id=sequence_record_id, name='',
                            description='')

                        flat_name2_dict[taxid] = sequence_record
                    else:
                        review_needed = True
                        # If there is more than one cluster we produce files
                        # for manual review.
                        review_dir = output_dir + '-review' + ps + name1 + ps
                        krio.prepare_directory(review_dir)
                        # Name the file by the id of the first record in the list
                        review_name = all_records_for_taxid_to_cluster[0].id.split('|')[0]

                        # Write all accessions we are unsure about to a file
                        review_ids_handle = open(review_dir+review_name+'_accessions.txt', 'w')
                        for r in all_records_for_taxid_to_cluster:
                            label = r.id.split('|')[0]
                            if label is not 'representative':
                                review_ids_handle.write(label+'\n')
                        review_ids_handle.close()

                        # Align all accessions we are unsure about, and write to a file for review
                        aln = kralign.align(
                            records=all_records_for_taxid_to_cluster,
                            program=aln_program,
                            options=aln_options,
                            program_executable=aln_program_exe)
                        AlignIO.write(aln, review_dir + review_name + '.phy', "phylip-relaxed")
                elif len(all_records_for_taxid) == 1:
                    # If there is only one record available, we hope for the
                    # best and keep it.
                    # TODO: Confirm this record is actually good.
                    flat_name2_dict[taxid] = all_records_for_taxid[0]
            name1_results.append({'results': flat_name2_dict, 'lrp': lrp})
            print('\n')
        # END LOOP ALL NAME2 --------------------------------------------------

        if not proceed_without_review and review_needed:
            continue

        #
        # This is where the final locus construction happens. Pieces are
        # arranged using lrps and a final locus consensus is produced.
        #
        # Keys will be taxid, values will be a list of dictionaries with keys:
        #   record
        #   lrp
        #

        final_dict = dict()
        name1_results.sort(key=lambda x: x['lrp'], reverse=False)
        for r_dict in name1_results:
            r = r_dict['results']
            lrp = r_dict['lrp']
            for taxid in r.keys():
                if taxid not in final_dict:
                    final_dict[taxid] = list()
                final_dict[taxid].append({'record': r[taxid], 'lrp': lrp})

        for i, taxid in enumerate(final_dict.keys()):
            # Produce a list of sequences that will represent a consensus of
            # the locus. We concatenate numbered lrps in order and stick them
            # at index 0 of the list. Further indexes will contain string named
            # lrps ('x')
            consensus_list = list()
            for r in final_dict[taxid]:
                lrp = r['lrp']
                record = r['record']

                # Based on previous ordering we assume that lrps in the list
                # will occur in order: numbers first, then strings:
                #   1,2,3...X,X,X

                # We have to concatenate records based on their lrp order
                if not isinstance(lrp, basestring):
                    if not len(consensus_list):
                        consensus_list.append(record)
                    else:
                        # for record in consensus list
                        cons_ids_1 = list()
                        desc_split = consensus_list[0].id.split('|')
                        desc_0_split = desc_split[0].split('$')
                        if desc_0_split[0] == 'consensus':
                            desc_0_split.remove('consensus')
                            for acc in desc_0_split:
                                cons_ids_1.append(acc)
                        else:
                            cons_ids_1.append(desc_split[0])

                        # for new record
                        cons_ids_2 = list()
                        desc_split = record.id.split('|')
                        desc_0_split = desc_split[0].split('$')
                        if desc_0_split[0] == 'consensus':
                            desc_0_split.remove('consensus')
                            for acc in desc_0_split:
                                cons_ids_2.append(acc)
                        else:
                            cons_ids_2.append(desc_split[0])

                        cons_ids = list(set(cons_ids_1 + cons_ids_2))
                        consensus_list[0] = consensus_list[0] + record

                        cons_id = '$'.join(cons_ids)
                        old_record_id_split = record.id.split('|')
                        sequence_record_id = (
                            'consensus$' + cons_id + '|' +
                            old_record_id_split[1] + '|' +
                            old_record_id_split[2] + '|' +
                            old_record_id_split[3] + '|' +
                            old_record_id_split[4])

                        sequence_record = SeqRecord.SeqRecord(
                            seq=consensus_list[0].seq,
                            id=sequence_record_id, name='', description='')

                        consensus_list[0] = sequence_record

                # If lrp is 'x' or any string, it means the record contains
                # sequence that covers entire locus that we are trying to
                # reconstruct, so we just put it in the list from which the
                # locus consensus will be determined
                else:
                    consensus_list.append(record)
            # END Produce a list of sequences that will represent a consensus
            # of the locus.

            cons_ids = list()
            for c in consensus_list:
                desc_split = c.id.split('|')
                desc_0_split = desc_split[0].split('$')
                if desc_0_split[0] == 'consensus':
                    desc_0_split.remove('consensus')
                    for acc in desc_0_split:
                        cons_ids.append(acc)
                else:
                    cons_ids.append(desc_split[0])
            cons_ids = list(set(cons_ids))

            taxa_dir = output_dir + '-taxa' + ps + name1 + ps
            krio.prepare_directory(taxa_dir)

            if len(consensus_list) > 1:
                aln = kralign.align(
                    records=consensus_list,
                    program=aln_program,
                    options=aln_options,
                    program_executable=aln_program_exe)
                # TODO: Better consensus. This consensus should include
                # correct ambiguities.
                # summary_aln = AlignInfo.SummaryInfo(aln)
                # consensus = summary_aln.dumb_consensus(threshold=0.001,
                #                                        ambiguous='N')
                consensus = kralign.consensus(aln, threshold=0.4, unknown='N')

                cons_id = '$'.join(cons_ids)
                old_record_id_split = consensus_list[0].id.split('|')
                sequence_record_id = (
                    'consensus$' + cons_id + '|' +
                    old_record_id_split[1] + '|' +
                    old_record_id_split[2] + '|' +
                    old_record_id_split[3] + '|' +
                    old_record_id_split[4])

                sequence_record = SeqRecord.SeqRecord(
                    seq=consensus, id=sequence_record_id, name='',
                    description='')

                cons_originals = list()
                for cid in cons_ids:
                    r = all_records[cid]
                    cons_originals.append(r)

                krbioio.write_sequence_file(cons_originals, taxa_dir + old_record_id_split[3] + '_originals.fasta', 'fasta')

                # Check if the names are not the same
                aln_name_set = set()
                for i, a in enumerate(aln):
                    if a.id in aln_name_set:
                        a.id = str(i) + '_' + a.id
                    aln_name_set.add(a.id)
                # if len(aln_name_set) == 1:
                #     krbioio.write_sequence_file(sequence_record, taxa_dir + old_record_id_split[3] + '.fasta', 'fasta')
                # else:
                AlignIO.write(aln, taxa_dir + old_record_id_split[3] + '.phy', "phylip-relaxed")
                results.append(sequence_record)
            else:
                krbioio.write_sequence_file(consensus_list, taxa_dir + consensus_list[0].id.split('|')[3] + '.fasta', 'fasta')
                results.append(consensus_list[0])

        for result in results:
            desc = result.id.split('|')

            acc_in_result = desc[0].split('$')
            if acc_in_result[0] == 'consensus':
                acc_in_result.remove('consensus')
            acc_in_result = ' '.join(acc_in_result)
            log_handle.write(desc[-2] + '\t' +
                             acc_in_result + '\n')

            desc = desc[-2]
            result.id = desc
            result.name = ''
            result.description = ''
        krbioio.write_sequence_file(results, output_file, 'fasta')

        print('\n\tAccepted sequences from', len(results), 'taxa.')

        log_handle.close()

    # END LOOP ALL NAME1 ------------------------------------------------------

    shutil.rmtree(temp_dir)

    krcl.show_cursor()

    if review_needed:
        print('\tWARNING! Please review accessions in review directory.')
        if not proceed_without_review:
            exit(0)


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
    number_of_gaps_between_loci
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
        matrix_output_file = output_dir + ps + 'presence' + '.csv'
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
