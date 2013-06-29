#!/usr/bin/env python

from __future__ import print_function
#from __future__ import unicode_literals


# Hacks!

def sol_species_from_voucher(voucher):

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
    # print(voucher)
    voucher = re.findall('LA\d+', voucher)
    if len(voucher) == 0:
        return(None)
    voucher = voucher[0]
    v_split_2 = re.findall('(\d+|[a-zA-Z]+)', voucher)
    if len(v_split_2) > 1:
        voucher = v_split_2[0] + "%04d" % (int(v_split_2[1]),)
    # print(voucher)
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

    # print(old, '->', new)
    # print(o_parsed, '->', n_parsed)

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
    print('\tPreparing output directory "', output_dir, '"', sep='')

    krio.prepare_directory(output_dir)

    for query_dict in queries:

        name1 = query_dict['name1']
        name2 = query_dict['name2']
        #locus = query_dict['locus']
        #minlen = query_dict['minlen']
        #feature_type = query_dict['ncbi_feature_type']
        #qualifier_label = query_dict['ncbi_qualifier_label']
        db = query_dict['database']
        query = query_dict['query']

        '''
        File name is genrated based on the search query (periods will be
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
                     #locus + file_name_sep +
                     #minlen + file_name_sep +
                     #feature_type + file_name_sep +
                     #qualifier_label + file_name_sep +
                     #db +
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
                   ### KR ###
                   # keeplist_records_file,
                   cutlist_taxonomy_file
                   ### KR ###
                   # keeplist_taxonomy_file
                   ):

    import os
    import krio
    import krbioio
    import krcl
    import krseq
    import krncbi

    ps = os.path.sep

    print('\nFiltering records.')

    print('\tPreparing output directory "', output_dir, '"', sep='')
    krio.prepare_directory(output_dir)

    cutlist_records = None
    ### KR ###
    # keeplist_records = None
    cutlist_taxonomy = None
    ### KR ###
    # keeplist_taxonomy = None

    if cutlist_records_file:

        cutlist_records = krio.read_table_file(
            path=cutlist_records_file,
            has_headers=False,
            headers=None,
            delimiter=',',
            quotechar='"',
            rettype='set')

    ### KR ###
    # if keeplist_records_file:

    #     keeplist_records = krio.read_table_file(
    #         path=keeplist_records_file,
    #         has_headers=False,
    #         headers=None,
    #         delimiter=',',
    #         quotechar='"',
    #         rettype='set')

    if cutlist_taxonomy_file:

        cutlist_taxonomy = krio.read_table_file(
            path=cutlist_taxonomy_file,
            has_headers=False,
            headers=None,
            delimiter=',',
            quotechar='"',
            rettype='set')

    ### KR ###
    # if keeplist_taxonomy_file:

    #     keeplist_taxonomy = krio.read_table_file(
    #         path=keeplist_taxonomy_file,
    #         has_headers=False,
    #         headers=None,
    #         delimiter=',',
    #         quotechar='"',
    #         rettype='set')

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

        print('\n\t\tFiltering records for ids.')
        ### KR ###
        # if keeplist_records or cutlist_records:
        if cutlist_records:
            records_count = len(records)
            for i, record in enumerate(records):
                krcl.print_progress(i + 1, records_count, 50, '\t')

                ### KR ###
                # if keeplist_records:
                #     if ((record.id in keeplist_records) or
                #         (krseq.get_annotation(record, 'gi') in
                #             keeplist_records)):
                #         records_filtered_id.append(record)
                if cutlist_records:
                    if ((record.id not in cutlist_records) and
                        (krseq.get_annotation(record, 'gi') not in
                            cutlist_records)):
                        records_filtered_id.append(record)
        else:
            records_filtered_id = records
        print('\t\tAccepted', len(records_filtered_id), 'records.')

        print('\n\t\tFiltering records for taxonomy.')
        ### KR ###
        # if keeplist_taxonomy or cutlist_taxonomy:
        if cutlist_taxonomy:
            records_count = len(records_filtered_id)
            for i, record in enumerate(records_filtered_id):
                krcl.print_progress(i + 1, records_count, 50, '\t')

                ### KR ###
                # if keeplist_taxonomy:
                #     if ((krncbi.get_ncbi_tax_id(record) in
                #         keeplist_taxonomy) or
                #         krseq.get_annotation(record, 'organism') in
                #             keeplist_taxonomy):
                #         records_filtered_tax.append(record)
                if cutlist_taxonomy:
                    if ((krncbi.get_ncbi_tax_id(record) not in
                        cutlist_taxonomy) and
                        krseq.get_annotation(record, 'organism') not in
                            cutlist_taxonomy):
                        records_filtered_tax.append(record)
        else:
            records_filtered_tax = records_filtered_id
        print('\t\tAccepted', len(records_filtered_tax), 'records.')

        # Write results
        ### KR ###
        # if keeplist_taxonomy or cutlist_taxonomy:
        if cutlist_taxonomy:
            krbioio.write_sequence_file(records_filtered_tax, output_file,
                                        'gb')
        else:
            krbioio.write_sequence_file(records_filtered_id, output_file, 'gb')

        krcl.show_cursor()


def extract_loci(search_results_dir, output_dir, queries,

                 ### KR ###
                 # sequence_samples,

                 ncbi_names_table,
                 ### KR ###
                 # min_similarity,
                 temp_dir, file_name_sep,
                 synonymy_table=None, auth_file=None, hacks=None,
                 hacks_data_location=None):

    '''
    Extract relevant loci from the search results, do some filtering by length
    and similarity.
    '''

    print('\nExtracting relevant loci.')

    import os
    from Bio import SeqRecord
    import krio
    import krbioio
    import krseq
    import krncbi
    import krcl
    import krbionames
    ### KR ###
    # import krusearch

    ps = os.path.sep

    # HACKS ###################################################################

    hack_sol_genus = None
    hack_sol_species = None

    found_previously = dict()
    do_not_repeat = list()

    if hacks and 'solanum' in hacks:
        hack_sol_species_set = krio.read_table_file(
            # path='..' + ps + 'data' + ps + 'sol_sp_with_vouchers',
            path=hacks_data_location,
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

    print('\tPreparing output directory "', output_dir, '"', sep='')
    krio.prepare_directory(output_dir)
    print('\tPreparing temporary directory "', temp_dir, '"', sep='')
    krio.prepare_directory(temp_dir)

    file_list = krio.parse_directory(search_results_dir, file_name_sep)

    # Iterate over search results
    for f in file_list:

        if not f['ext'].startswith('gb'):
            continue

        # Stores all loci that have passed quality control
        loci = list()
        # Stores all loci that have failed quality control
        #loci_excluded = list()

        # Read search results
        records = krbioio.read_sequence_file(f['path'], 'genbank')

        # Get information from the search result file name
        file_name = f['name']
        name1 = f['split'][0]
        name2 = f['split'][1]
        #locus = f['split'][2]
        #minlen = int(f['split'][3])
        #feature_type = f['split'][4]
        #qualifier_label = f['split'][5]

        locus = ''
        minlen = 0
        feature_type = ''
        qualifier_label = ''
        match_stringency = ''

        for query_dict in queries:
            if name1 == query_dict['name1'] and name2 == query_dict['name2']:
                #name1 = query_dict['name1']
                #name2 = query_dict['name2']
                force_rev_comp = query_dict['force_rev_comp']
                locus = query_dict['locus']
                minlen = int(query_dict['minlen'])
                feature_type = query_dict['ncbi_feature_type']
                qualifier_label = query_dict['ncbi_qualifier_label']
                match_stringency = query_dict['match_stringency']
                break
                #db = query_dict['database']
                #query = query_dict['query']

        log_file = output_dir + ps + file_name + '.log'
        log_handle = open(log_file, 'w')

        output_file = output_dir + ps + file_name + '.fasta'
        #output_file_excluded = (output_dir + ps + file_name + file_name_sep +
        #    'excluded.fasta')

        records_count = len(records)

        print('\n\tProcessing: ', name1, ' ', name2, ' / ', locus, sep='')

        krcl.hide_cursor()

        # Locus may have different names
        locus = locus.split(',')

        for i, record in enumerate(records):
            # print(record.id)
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
            # ToDo: this should never occur, and occured only once
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
            if hacks and 'solanum' in hacks:

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
                        # print(voucher)
                        if voucher not in do_not_repeat:
                            if voucher in found_previously.keys():
                                hack_species_found = found_previously[voucher]
                            else:
                                hack_species_found = sol_species_from_voucher(
                                    voucher)
                            if hack_species_found:
                                # print('FOUND')
                                found_previously[voucher] = hack_species_found
                                acc_name = hack_species_found[1]
                                acc_name['id'] = (
                                    'tgrc-' + hack_species_found[2])
                                acc_name['status'] = 'TGRC'
                                acc_name_flat = (
                                    krbionames.flatten_organism_name(
                                        acc_name, '_'))
                            else:
                                do_not_repeat.append(voucher)
            ###################################################################

            loose_mode = False

            if not hack_species_found and synonymy_table and auth_file:
                #print('--- --- --- --- --- ---')
                #print('O', organism)

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
                        if acc_name and acc_name['status'].lower() == 'acc':
                            loose_mode = True
                            found_match = True
                            break

                # If no, then look for the next best thing "prov"
                if not found_match:
                    for oa in organism_authority:
                        acc_name = krbionames.accepted_name(
                            name=oa,
                            synonymy_table=synonymy_table,
                            auth_file=auth_file,
                            allow_loose_matching=True)
                        if acc_name and acc_name['status'].lower() == 'prov':
                            loose_mode = True
                            found_match = True
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
                            acc_name['status'].lower() == 'as' or
                            acc_name['status'].lower() == 'nn' or
                            acc_name['status'].lower() == 'unc' or
                                acc_name['status'].lower() == 'unr')):
                            loose_mode = True
                            found_match = True
                            break

                if not found_match:
                    log_handle.write(
                        name1 + '_' + name2 + '\t' +
                        record.id + '\t' +
                        organism.replace('_', ' ') + '\t' +
                        tax_id + '\t'
                        'No taxonomic match.\n')
                    #print('___!!!___!!!_______ NO MATCH _______!!!___!!!___')
                    continue

                acc_name_flat = krbionames.flatten_organism_name(acc_name, '_')
                #print('A', acc_name_flat)
                if loose_mode:
                    #print('___!!!___!!!_______ LOOSE _______!!!___!!!___')
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

            ### KR ###
            # # We will try to cluster this sequence with a sample at the
            # # relatively low similarity treshold, to weed out sequences
            # # that have nothing to do with what we are looking for.
            # to_cluster = [sequence_record, sequence_samples[name1]]
            # cluster_dict = krusearch.cluster_records(
            #     to_cluster,
            #     min_similarity, temp_dir)

            if acc_name['status'] == '':
                log_handle.write(
                    name1 + '_' + name2 + '\t' +
                    record.id + '\t' +
                    organism.replace('_', ' ') + '\t' +
                    tax_id + '\t' +
                    'No taxonomic match.\n')
                #loci_excluded.append(sequence_record)
            elif len(seq) <= minlen:
                log_handle.write(
                    name1 + '_' + name2 + '\t' +
                    record.id + '\t' +
                    organism.replace('_', ' ') + '\t' +
                    tax_id + '\t' +
                    'Sequence is too short.\n')
                #loci_excluded.append(sequence_record)

            ### KR ###
            # # If the sequences are similar enough, there will be only one
            # # cluster
            # elif len(cluster_dict.keys()) != 1:
            #     log_handle.write(
            #         name1 + '_' + name2 + '\t' +
            #         record.id + '\t' +
            #         organism.replace('_', ' ') + '\t' +
            #         tax_id + '\t' +
            #         'Sequence is too dissimilar from a sample sequence.\n')
            #     #loci_excluded.append(sequence_record)

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
                    ###
                    # print('TGRC', record.id, tax_id,
                    #       organism.replace('_', ' '), '->',
                    #       acc_name_flat.replace('_', ' '))
                    ###
                    log_handle.write(
                        name1 + '_' + name2 + '\t' +
                        record.id + '\t' +
                        organism.replace('_', ' ') + '\t' +
                        tax_id + '\t' +
                        'Match TGRC' + '\t' +
                        acc_name_flat.replace('_', ' ') + '\n')
                else:
                    ###
                    # print(record.id, acc_name_flat)
                    ###
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
        #krbioio.write_sequence_file(loci_excluded, output_file_excluded,
        #    'fasta')

        print('\n\tAccepted', len(loci), 'sequences.')
        #print('\n\tRejected', len(loci_excluded), 'sequences.')

        krcl.show_cursor()

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

    os.removedirs(temp_dir)


def one_locus_per_organism(
    extracted_results_dir,
    output_dir,
    queries,
    min_similarity,
    cutlist_records_file,
    keeplist_records_file,
    temp_dir,
    file_name_sep,
    aln_program,
    threads
):

    '''
    Produce single sequence per locus, per organism.
    '''

    print('\nProduce single sequence per locus, per organism.')

    import os
    from Bio.Align import AlignInfo
    from Bio import SeqRecord
    import krio
    import krbioio
    import krusearch
    import krcl
    import kralign

    ps = os.path.sep

    print('\tPreparing output directory "', output_dir, '"', sep='')
    krio.prepare_directory(output_dir)

    print('\tPreparing temporary directory "', temp_dir, '"', sep='')
    krio.prepare_directory(temp_dir)

    name_dict = dict()

    file_list = krio.parse_directory(extracted_results_dir, file_name_sep)

    ### KR ###
    cutlist_records = None
    if cutlist_records_file:

        cutlist_records = krio.read_table_file(
            path=cutlist_records_file,
            has_headers=False,
            headers=None,
            delimiter=',',
            quotechar='"',
            rettype='set')

    keeplist_records = None
    if keeplist_records_file:

        keeplist_records = krio.read_table_file(
            path=keeplist_records_file,
            has_headers=False,
            headers=None,
            delimiter=',',
            quotechar='"',
            rettype='set')
    ### KR ###

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
    #       name_dict: dictionary with key name1, each value is a list of
    #       dictionaries with keys:
    #           name2:      specific name (string)
    #           records:    list of BioSeqRecord objects (list)
    #           lrp:        locus relative position (int or 'X')
    #
    for f in file_list:

        if not f['ext'].startswith('fasta'):
            continue

        # Not sure if this will be used
        #if f['split'][-1].startswith('excluded'):
        #    continue

        name1 = f['split'][0]
        name2 = f['split'][1]
        locus_relative_position = 0

        for query_dict in queries:
            if name1 == query_dict['name1'] and name2 == query_dict['name2']:
                locus_relative_position = query_dict['locus_relative_position']
                break

        if not name1 in name_dict:
            name_dict[name1] = list()

        records = krbioio.read_sequence_file(f['path'], 'fasta')
        name_dict[name1].append({
            'name2': name2,
            'records': records,
            'lrp': locus_relative_position})

    # Because for each locus, there could be several records from the same
    # species, we need to pick the best representative sequence. Only one
    # per locus, per organism. Each name1 will have one and only one file
    # associated with it.

    # LOOP ALL NAME1 ----------------------------------------------------------
    for name1 in name_dict.keys():
        print('\n\tProcessing', name1)

        results = list()

        log_file = output_dir + ps + name1 + '.log'
        log_handle = open(log_file, 'w')
        output_file = output_dir + ps + name1 + '.fasta'

        krcl.hide_cursor()

        # Each row will contain a dictionary with a list of results and a
        # relative locus position:
        name1_results = list()

        # All name2 associated with name1
        list_of_name2 = name_dict[name1]

        # LOOP ALL NAME2 ------------------------------------------------------
        for name2_dict in list_of_name2:

            name2 = name2_dict['name2']
            records = name2_dict['records']
            lrp = name2_dict['lrp']

            if lrp.lower() != 'x':
                lrp = int(lrp)

            ### KR ###
            # Cluster top X percent of records within a locus within NAME2
            # Get a representative consensus sequence
            records.sort(key=lambda x: len(x), reverse=True)
            records_count = len(records)
            fraction_to_cluster = 0.01
            number_to_cluster = int(records_count * fraction_to_cluster)
            number_to_cluster = max(10, number_to_cluster)
            records_to_cluster = list()
            for i, record in enumerate(records):
                if i == number_to_cluster:
                    break
                records_to_cluster.append(record)
            cluster_dict = krusearch.cluster_records(
                records=records_to_cluster,
                identity_threshold=0.75,  ### KR ###
                temp_dir=temp_dir,
                sorted_input=False,
                algorithm='smallmem',  # fast smallmem
                strand='both',  # plus both
                threads=1,
                quiet=True,
                program='usearch6.1.544_i86linux32',
                heuristics=False,
                query_coverage=0.1,
                target_coverage=0.1,
                sizein=False,
                sizeout=False,
                usersort=False
            )
            clusters = list()
            for k in cluster_dict.keys():
                clusters.append(cluster_dict[k])
            clusters.sort(key=lambda x: len(x), reverse=True)
            largest_cluster = clusters[0]
            largest_cluster_names = list()
            for n in largest_cluster:
                largest_cluster_names.append(n[1])
            # print(largest_cluster_names)
            records_to_align = list()
            for record in records:
                # print(record.id)
                if record.id in largest_cluster_names:
                    records_to_align.append(record)
            ### KR ###
            aln = kralign.align(records_to_align, 'mafft')
            # from Bio import AlignIO
            # AlignIO.write(aln, '/home/karolis/Dropbox/Code/test/sol-out-NEW/temp/xxx.fasta', "fasta")
            summary_aln = AlignInfo.SummaryInfo(aln)
            consensus = summary_aln.dumb_consensus(
                threshold=0.001, ambiguous='N')
            rep_cons_record = SeqRecord.SeqRecord(
                seq=consensus, id='representative', name='',
                description='')
            # print(consensus)

            print('\n\t\tProcessing', name1, name2)

            # Keys will be taxid and value will be the list of all sequences
            # for name2 with that taxid
            ### ### ### ### ### ### ### ###
            ### CANNOT USE TAXID BECAUSE OF SOLANUM HACK !!!
            ### ### ### ### ### ### ### ###
            tax_records_dict_name2 = dict()

            for record in records:
                # print(record.id.split('|')[0], len(record))
                # taxid is acc_name_flat
                taxid = record.description.split('|')[3]
                if not taxid in tax_records_dict_name2:
                    tax_records_dict_name2[taxid] = list()
                tax_records_dict_name2[taxid].append(record)

            # Will hold one sequence record per taxid, keys are taxid
            all_taxid_for_name2_dict = dict()

            records_count = len(tax_records_dict_name2.keys())

            for i, taxid in enumerate(tax_records_dict_name2.keys()):

                krcl.print_progress(i + 1, records_count, 50, '\t\t')

                ### KR ### ???
                ### Make sure to remove id field from the log file or put
                ### multiple ids as each species can have multiple ids now

                tax_records_name2_prefiltered = tax_records_dict_name2[taxid]
                tax_records_name2 = list()
                ### KR ###
                for tr in tax_records_name2_prefiltered:
                    if tr.id.split('|')[0] not in cutlist_records:
                        tax_records_name2.append(tr)
                ### KR ###

                if len(tax_records_name2) > 1:
                    # ...we cluster these sequences and hope for only one
                    # cluster

                    tax_records_name2_to_cluster = list()
                    ### KR ###
                    for tr in tax_records_name2:
                        if tr.id.split('|')[0] not in keeplist_records:
                            tax_records_name2_to_cluster.append(tr)
                    ### KR ###
                    tax_records_name2_to_cluster.append(rep_cons_record)

                    cluster_dict = None
                    if len(tax_records_name2_to_cluster) > 1:
                        cluster_dict = krusearch.cluster_records(
                            records=tax_records_name2_to_cluster,
                            # identity_threshold=min_similarity,
                            identity_threshold=0.75,
                            temp_dir=temp_dir,
                            sorted_input=False,
                            algorithm='smallmem',  # fast smallmem
                            strand='both',  # plus both
                            threads=1,
                            quiet=True,
                            program='usearch6.1.544_i86linux32',
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
                        aln = kralign.align(tax_records_name2, aln_program)
                        summary_aln = AlignInfo.SummaryInfo(aln)

                        ### KR ### Consensus!
                        consensus = summary_aln.dumb_consensus(
                            threshold=0.001, ambiguous='N')

                        # We want the consensus record description to reflect
                        # all the sequence ids that were combined to form it
                        cons_ids = list()
                        for n in tax_records_name2:
                            cons_ids.append(n.id.split('|')[0])
                        new_cons_id = '_'.join(cons_ids)

                        old_record_id_split = (
                            tax_records_name2[0].id.split('|'))

                        sequence_record_id = (
                            'consensus_' + new_cons_id + '|' +
                            old_record_id_split[1] + '|' +
                            old_record_id_split[2] + '|' +
                            old_record_id_split[3] + '|' +
                            old_record_id_split[4])

                        # Produce a biopython sequence record object
                        sequence_record = SeqRecord.SeqRecord(
                            seq=consensus, id=sequence_record_id, name='',
                            description='')

                        all_taxid_for_name2_dict[taxid] = sequence_record
                    else:
                        ### KR ###
                        # What to do if there is more than one cluster?

                        review_dir = output_dir + '-review' + ps + name1 + ps
                        krio.prepare_directory(review_dir)

                        review_name = tax_records_name2_to_cluster[0].id.split('|')[0]

                        review_ids_handle = open(review_dir+review_name+'.tsv', 'w')

                        for r in tax_records_name2_to_cluster:
                            label = r.id.split('|')[0]
                            if label is not 'representative':
                                review_ids_handle.write(label+'\n')

                        review_ids_handle.close()

                        from Bio import AlignIO
                        aln = kralign.align(tax_records_name2_to_cluster, 'mafft')
                        AlignIO.write(aln, review_dir + review_name + '.phy', "phylip-relaxed")

                        ### KR ###
                        # clusters = list()
                        # for k in cluster_dict.keys():
                        #     cluster = cluster_dict[k]
                        #     count = len(cluster)
                        #     clusters.append((count, cluster))

                        # clusters = cluster_dict.values()
                        # clusters.sort(key=lambda x: len(x), reverse=True)

                        # # # if len(clusters) > 1:
                        # for x, c in enumerate(clusters):
                        #     print(x, c[0][1])
                        #     for s in c:
                        #         print(s[1])

                        # print('----------------------------------------------')

                        # if len(clusters) > 1:
                        #     for c in clusters:
                        #         print(c)

                        #     print('==============================================')
                        ### KR ###

                else:
                    all_taxid_for_name2_dict[taxid] = tax_records_name2[0]

            name1_results.append({'results': all_taxid_for_name2_dict,
                                 'lrp': lrp})

            print('\n')

        # END LOOP ALL NAME2 --------------------------------------------------

        # Keys will be taxid, values will be a list of dictionaries with keys:
        #   record
        #   lrp
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

            #print(i+1)

            # -----------------------------------------------------------------
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
                        desc_0_split = desc_split[0].split('_')
                        if desc_0_split[0] == 'consensus':
                            desc_0_split.remove('consensus')
                            for acc in desc_0_split:
                                cons_ids_1.append(acc)
                        else:
                            cons_ids_1.append(desc_split[0])

                        # for new record
                        cons_ids_2 = list()
                        desc_split = record.id.split('|')
                        desc_0_split = desc_split[0].split('_')
                        if desc_0_split[0] == 'consensus':
                            desc_0_split.remove('consensus')
                            for acc in desc_0_split:
                                cons_ids_2.append(acc)
                        else:
                            cons_ids_2.append(desc_split[0])

                        cons_ids = list(set(cons_ids_1 + cons_ids_2))
                        consensus_list[0] = consensus_list[0] + record

                        cons_id = '_'.join(cons_ids)
                        old_record_id_split = record.id.split('|')
                        sequence_record_id = (
                            'consensus_' + cons_id + '|' +
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
            # -----------------------------------------------------------------

            #    print(record.id, lrp, len(r['record']))
            #print('--- --- --- --- --- --- --- --- --- --- --- ---')
            #for c in consensus_list:
            #    print(c.id, len(c))
            #print('--- --- --- --- --- --- --- --- --- --- --- ---')

            cons_ids = list()
            for c in consensus_list:
                desc_split = c.id.split('|')
                desc_0_split = desc_split[0].split('_')
                if desc_0_split[0] == 'consensus':
                    desc_0_split.remove('consensus')
                    for acc in desc_0_split:
                        cons_ids.append(acc)
                else:
                    cons_ids.append(desc_split[0])
            cons_ids = list(set(cons_ids))

            if len(consensus_list) > 1:
                aln = kralign.align(consensus_list, aln_program)
                summary_aln = AlignInfo.SummaryInfo(aln)
                consensus = summary_aln.dumb_consensus(threshold=0.001,
                                                       ambiguous='N')

                cons_id = '_'.join(cons_ids)
                old_record_id_split = consensus_list[0].id.split('|')
                sequence_record_id = (
                    'consensus_' + cons_id + '|' +
                    old_record_id_split[1] + '|' +
                    old_record_id_split[2] + '|' +
                    old_record_id_split[3] + '|' +
                    old_record_id_split[4])

                sequence_record = SeqRecord.SeqRecord(
                    seq=consensus, id=sequence_record_id, name='',
                    description='')

                results.append(sequence_record)
            #    print(sequence_record.id, len(sequence_record))
            #    print('=== === === === === === === === === === === ===')

            else:
                results.append(consensus_list[0])
            #    print(consensus_list[0].id, len(consensus_list[0]))
            #    print('=== === === === === === === === === === === ===')

        krcl.show_cursor()

        for result in results:
            desc = result.id.split('|')

            acc_in_result = desc[0].split('_')
            if acc_in_result[0] == 'consensus':
                acc_in_result.remove('consensus')
            acc_in_result = ' '.join(acc_in_result)
            log_handle.write(desc[-2] + '\t' +
                             # desc[-1] + '\t' +
                             acc_in_result + '\n')

            desc = desc[-2]  # + '|' + desc[-1]
            result.id = desc
            result.name = ''
            result.description = ''
        krbioio.write_sequence_file(results, output_file, 'fasta')

        print('\n\tAccepted', len(results), 'sequences.')

        log_handle.close()

    # END LOOP ALL NAME1 ------------------------------------------------------

    os.removedirs(temp_dir)


def align_loci(processed_results_dir, output_dir, program, options,
               spacing, temp_dir, order):

    '''
    Align individual loci, then concatenate alignments.

    spacing - the number of gaps to insert between individual alignments.
    '''

    print('\nAlign individual loci, then concatenate alignments.')

    import os
    import krio
    import krbioio
    import kralign
    import copy

    ps = os.path.sep

    #alignments = []

    print('\tPreparing output directory "', output_dir, '"', sep='')
    krio.prepare_directory(output_dir)
    print('\tPreparing temporary directory "', temp_dir, '"', sep='')
    krio.prepare_directory(temp_dir)

    file_list = krio.parse_directory(processed_results_dir, ' ')

    order_list = [x.strip() for x in order.split(',')]
    alignments = [x.strip() for x in order.split(',')]

    for f in file_list:

        if not f['ext'].startswith('fasta'):
            continue

        file_name = f['name']
        output_file = output_dir + ps + file_name + '.fasta'
        records = krbioio.read_sequence_file(f['path'], 'fasta')

        # Align each locus individually first.
        print('\n\tAligning', file_name)
        aln = kralign.align(records, program, options=options)
        if aln:
            krbioio.write_alignment_file(aln, output_file, 'fasta')
            #alignments.append(aln)
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
        concatenated = kralign.concatenate(raw_alignments, int(spacing))
        concatenated_output_file = output_dir + ps + 'concatenated' + '.fasta'
        krbioio.write_alignment_file(concatenated, concatenated_output_file,
                                     'fasta')

    os.removedirs(temp_dir)


# End pipeline functions ------------------------------------------------------


# if __name__ == '__main__':

#     import sys
#     import argparse
#     import os
#     import krio
#     import krbioio

#     parser = argparse.ArgumentParser()
#     parser.add_argument('-t', '--test', action='store_true', help='Run tests.')
#     parser.add_argument(
#         '-c', '--command', type=unicode,
#         choices=['search_and_download', 'filter_records', 'extract_loci',
#         'one_locus_per_organism', 'align_loci'], help='Run a command.')
#     parser.add_argument('-q', '--query', type=unicode, help='Query file path.')
#     parser.add_argument(
#         '-o', '--output', type=unicode,
#         help='Output directory path.')
#     parser.add_argument(
#         '-s', '--sep', type=unicode,
#         help='Output file name separator.')
#     parser.add_argument('-e', '--email', type=unicode, help='Email.')
#     parser.add_argument(
#         '-i', '--input', type=unicode,
#         help='Input directory path.')
#     parser.add_argument(
#         '--ncbinames', type=unicode,
#         help='NCBI organism names file.')
#     parser.add_argument('--synonymy', type=unicode, help='Synonymy file.')
#     parser.add_argument(
#         '--authority', type=unicode,
#         help='Authority alternates file.')
#     parser.add_argument(
#         '--samples', type=unicode,
#         help='Sequence samples file in FASTA format.')
#     parser.add_argument(
#         '--similarity', type=float,
#         help='Minimum sequence similarity.')
#     parser.add_argument(
#         '--tempdir', type=unicode,
#         help='Temporary directory path.')
#     parser.add_argument(
#         '--alignprog', type=unicode,
#         choices=['mafft', 'einsi', 'linsi', 'muscle'],
#         help='Alignment program to be used.')
#     parser.add_argument(
#         '--threads', type=int,
#         help='Number of CPU cores to use.')
#     parser.add_argument(
#         '--alignspacing', type=int,
#         help='Number of gaps to add between alignments \
#         in concatenated alignment.')
#     parser.add_argument('--clrec', type=unicode, help='Cutlist records file.')
#     parser.add_argument('--klrec', type=unicode, help='Keeplist records file.')
#     parser.add_argument('--cltax', type=unicode, help='Cutlist taxonomy file.')
#     parser.add_argument('--kltax', type=unicode,
#                         help='Keeplist taxonomy file.')
#     parser.add_argument('--alignorder', type=unicode,
#                         help='Order of loci in the concatenated alignment.')

#     args = parser.parse_args()

#     PS = os.path.sep

#     if args.test:

#         print('Running tests.')

#         # Tests

#         # parse_directory
#         # t_parse_directory = parse_directory('testdata' + PS +
#         #                                     'parse_directory', '$')

#         # t_parse_directory = parse_directory(
#         #     '/data/gbs-new/02-demultiplexed-fastq', ' ')
#         # for d in t_parse_directory:
#         #     print(d)

#         v1 = 'LA1964'
#         v2 = 'LA2744'

#         s1 = sol_species_from_voucher(v1)
#         print(s1)

#         s2 = sol_species_from_voucher(v2)
#         print(s2)

#         sys.exit(0)

#     else:

#         queries = krio.read_table_file(args.query, has_headers=True,
#                                        headers=None, delimiter='\t')

#         if args.command == 'search_and_download':

#             search_and_download(queries, args.output + PS, args.sep,
#                                 args.email)

#         if args.command == 'filter_records':
#             filter_records(
#                 search_results_dir=args.input,
#                 output_dir=args.output,
#                 cutlist_records_file=args.clrec,
#                 keeplist_records_file=args.klrec,
#                 cutlist_taxonomy_file=args.cltax,
#                 keeplist_taxonomy_file=args.kltax
#             )

#         if args.command == 'extract_loci':

#             ncbi_names = None
#             synonymy_table = None
#             sequence_samples = None

#             if args.ncbinames:

#                 ncbi_names = krio.read_table_file(
#                     path=args.ncbinames,
#                     has_headers=False,
#                     headers=('tax_id', 'name_txt', 'unique_name',
#                              'name_class'),
#                     delimiter='\t|',
#                     quotechar=None,
#                     stripchar='"',
#                     rettype='dict')

#             if args.authority and args.synonymy:

#                 synonymy_table = krio.read_table_file(
#                     args.synonymy, has_headers=True, headers=None,
#                     delimiter=',')

#             if args.samples:

#                 sequence_samples = krbioio.read_sequence_file(
#                     args.samples, 'fasta', ret_type='dict')

#             extract_loci(
#                 search_results_dir=args.input,
#                 output_dir=args.output,
#                 queries=queries,
#                 sequence_samples=sequence_samples,
#                 ncbi_names_table=ncbi_names,
#                 synonymy_table=synonymy_table,
#                 auth_file=args.authority,
#                 min_similarity=args.similarity,
#                 temp_dir=args.tempdir,
#                 file_name_sep=args.sep
#             )

#         if args.command == 'one_locus_per_organism':

#             one_locus_per_organism(
#                 extracted_results_dir=args.input,
#                 output_dir=args.output,
#                 queries=queries,
#                 min_similarity=args.similarity,
#                 temp_dir=args.tempdir,
#                 file_name_sep=args.sep,
#                 aln_program=args.alignprog,
#                 threads=args.threads
#             )

#         if args.command == 'align_loci':

#             align_loci(
#                 processed_results_dir=args.input,
#                 output_dir=args.output,
#                 program=args.alignprog,
#                 threads=args.threads,
#                 spacing=args.alignspacing,
#                 temp_dir=args.tempdir,
#                 order=args.alignorder
#             )
