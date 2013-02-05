#!/usr/bin/env python

from __future__ import print_function
#from __future__ import unicode_literals

# Utility functions used by the pipeline functions.

def parse_directory(path, file_name_sep):

    '''
    Will parse a directory at a given path and return a list of dictionary
    objects with keys:
        name: name of a file without extension
        ext: file extension
        full: file name with extension
        path: full file path, relative to the input path
        split: file name split using file_name_sep input variable
    '''
    
    import os
    
    ps = os.path.sep
    
    file_list = os.listdir(path)
    
    return_list = list()
    
    for f in file_list:
    
        file_name = os.path.splitext(f)[0]
        file_ext = os.path.splitext(f)[1].split('.')[1]
        file_name_split = file_name.split(file_name_sep)
    
        file_dict = dict()
    
        file_dict['name'] = file_name
        file_dict['ext'] = file_ext
        file_dict['full'] = f
        file_dict['path'] = path + ps + f
        file_dict['split'] = file_name_split
    
        return_list.append(file_dict)
    
    return return_list

# Pipeline functions ----------------------------------------------------------

def search_and_download(queries, output_dir, file_name_sep, email):

    '''
    This will search NCBI and download sequences.
    '''

    import krio
    import krncbi

    print('\nSearching NCBI.')
    print('\tPreparing output directory "', output_dir, '"', sep='')

    krio.prepare_directory(output_dir)

    for query_dict in queries:

        name1 = query_dict['name1']
        name2 = query_dict['name2']
        locus = query_dict['locus']
        minlen = query_dict['minlen']
        feature_type = query_dict['ncbi_feature_type']
        qualifier_label = query_dict['ncbi_qualifier_label']
        db = query_dict['database']
        query = query_dict['query']

        print('\n\tSearching for:', name1, name2, 'Locus:', locus, 'in', db,
            'database.\n', sep=' ')

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

        file_name = (output_dir +
                     name1 + file_name_sep +
                     name2 + file_name_sep +
                     locus + file_name_sep +
                     minlen + file_name_sep +
                     feature_type + file_name_sep +
                     qualifier_label + file_name_sep +
                     db + '.gb')
        # Search NCBI.
        result_uids = krncbi.esearch(query, db, email)
        # Download records.
        krncbi.download_sequence_records(file_name, result_uids, db, email)

def extract_loci(search_results_dir, output_dir, sequence_samples,
    ncbi_names_table, min_similarity, temp_dir, file_name_sep,
    synonymy_table=None, auth_file=None):

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
    import krusearch

    ps = os.path.sep

    print('\tPreparing output directory "', output_dir, '"', sep='')
    krio.prepare_directory(output_dir)
    print('\tPreparing temporary directory "', temp_dir, '"', sep='')
    krio.prepare_directory(temp_dir)

    file_list = parse_directory(search_results_dir, file_name_sep)

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
        locus = f['split'][2]
        minlen = int(f['split'][3])
        feature_type = f['split'][4]
        qualifier_label = f['split'][5]
        
        log_file = output_dir + ps + file_name + '.log'
        log_handle = open(log_file, 'w')

        output_file = output_dir + ps + file_name + '.fasta'
        output_file_excluded = (output_dir + ps + file_name + file_name_sep +
            'excluded.fasta')
        
        records_count = len(records)

        print('\n\tProcessing: ', name1, ' ', name2, ' / ', locus, sep='')
        
        krcl.hide_cursor()

        # Locus may have different names
        locus = locus.split(',')

        for i, record in enumerate(records):
            krcl.print_progress(i+1, records_count, 50, '\t')
            
            # genbank records contain "features" which conatin annotation
            #   information. Here we look for the feature that contains our
            #   target locus and note its index

            feature_indexes = list()

            for l in locus:
                fi = krseq.get_features_with_qualifier(record=record,
                    qualifier_label=qualifier_label, qualifier=l.strip(),
                    feature_type=feature_type)
                feature_indexes = feature_indexes + fi
            # ToDo: this should never occur, and occured only once
            # Same gene annotated more than once
            if len(feature_indexes) > 1:
                log_handle.write(record.id + '\tMore than one locus annotation.\n')
                continue
            if len(feature_indexes) == 0:
                log_handle.write(record.id + '\tNo locus annotation.\n')
                continue
            # There should be only one matching index.
            feature = record.features[feature_indexes[0]]
            start = int(feature.location.start)
            end = int(feature.location.end)
            strand = int(feature.location.strand)
            # Extract relevant region
            seq = record.seq[start:end]
            # If the feature is in reverse orientation, reverse-complement.
            if strand == -1:
                seq = seq.reverse_complement()
            
            # Deal with the organism name
            tax_id = krncbi.get_ncbi_tax_id(record)

            ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
            
            organism = krseq.get_annotation(record, 'organism')
            organism = organism.replace(' ', '_')

            acc_name = None
            acc_name_flat = None
            organism_authority = None

            if synonymy_table and auth_file:

                # A list of organism names based on NCBI taxid. This is a sorted
                #   list with the most complete names at lower indexes.
                organism_authority = krncbi.names_for_ncbi_taxid(tax_id,
                    ncbi_names_table, sorting='authority')

                # Iterate over the list of NCBI names and try to resolve accepted
                #   name.

                # First look if there is a best possible match "acc"

                found_match = False

                for oa in organism_authority:
                    acc_name = krbionames.accepted_name(oa, synonymy_table,
                        auth_file)
                    if acc_name['status'].lower() == 'acc':
                        found_match = True
                        break

                # If no, then look for the next best thing "prov"

                if not found_match:
                    for oa in organism_authority:
                        acc_name = krbionames.accepted_name(oa, synonymy_table,
                            auth_file)
                        if acc_name['status'].lower() == 'prov':
                            found_match = True
                            break

                # Otherwise, let's find something that isn't blank

                if not found_match:
                    for oa in organism_authority:
                        acc_name = krbionames.accepted_name(oa, synonymy_table,
                            auth_file)
                        if (acc_name['status'].lower() == 'as' or
                            acc_name['status'].lower() == 'nn' or
                            acc_name['status'].lower() == 'unc' or
                            acc_name['status'].lower() == 'unr'):
                            found_match = True
                            break

                if not found_match:
                    log_handle.write(record.id + '\t' + organism.replace('_', ' ') + '\t' + tax_id + '\tNo taxonomic match.\n')
                    continue

                acc_name_flat = krbionames.flatten_organism_name(acc_name, '_')

            else:

                acc_name = dict()
                acc_name['id'] = tax_id
                acc_name['status'] = 'NA'
                organism_authority = krncbi.names_for_ncbi_taxid(tax_id,
                    ncbi_names_table, sorting='class')
                acc_name_flat = krbionames.flatten_organism_name(organism_authority[0], '_')

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
            # We will try to cluster this sequence with a sample at the
            # relatively low similarity treshold, to weed out sequences
            # that have nothing to do with what we are looking for.
            to_cluster = [sequence_record, sequence_samples[name1]]
            cluster_dict = krusearch.cluster_records(to_cluster,
                min_similarity, temp_dir)

            if acc_name['status'] == '':
                log_handle.write(record.id + '\t' + organism.replace('_', ' ') + '\t' + tax_id + '\tNo taxonomic match.\n')
                #loci_excluded.append(sequence_record)
            elif len(seq) <= minlen:
                log_handle.write(record.id + '\t' + organism.replace('_', ' ') + '\t' + tax_id + '\tSequence is too short.\n')
                #loci_excluded.append(sequence_record)
            # If the sequences are similar enough, there will be only one
            # cluster
            elif len(cluster_dict.keys()) != 1:
                log_handle.write(record.id + '\t' + organism.replace('_', ' ') + '\t' + tax_id + '\tSequence is too dissimilar from a sample sequence.\n')
                #loci_excluded.append(sequence_record)
            else:
                loci.append(sequence_record)

        log_handle.close()

        # Write results
        krbioio.write_sequence_file(loci, output_file, 'fasta')
        #krbioio.write_sequence_file(loci_excluded, output_file_excluded,
        #    'fasta')

        print('\n\tAccepted', len(loci), 'sequences.')
        #print('\n\tRejected', len(loci_excluded), 'sequences.')

        krcl.show_cursor()

    os.removedirs(temp_dir)

def one_locus_per_organism(extracted_results_dir, output_dir, min_similarity,
    temp_dir, file_name_sep):

    '''
    Produce single sequence per locus, per organism.
    '''

    print('\nProduce single sequence per locus, per organism.')

    import os
    import krio
    import krbioio
    import krusearch
    import krcl

    ps = os.path.sep

    print('\tPreparing output directory "', output_dir, '"', sep='')
    krio.prepare_directory(output_dir)
    
    print('\tPreparing temporary directory "', temp_dir, '"', sep='')
    krio.prepare_directory(temp_dir)

    locus_dict = dict()

    file_list = parse_directory(extracted_results_dir, file_name_sep)
    
    for f in file_list:

        if not f['ext'].startswith('fasta'):
            continue

        if f['split'][-1].startswith('excluded'):
            continue

        file_name = f['name']
        name = f['split'][0]

        if not locus_dict.has_key(name):
            locus_dict[name] = list()

        records = krbioio.read_sequence_file(f['path'], 'fasta')
        locus_dict[name] = locus_dict[name] + records

    # Because for each locus, there could be several records from the same
    # species, we need to pick the best representative sequence. Only one
    # per locus, per organism.
    for locus_name in locus_dict.keys():
        print('\n\tProcessing', locus_name)

        output_file = output_dir + ps + locus_name + '.fasta'
        results = list()
        records = locus_dict[locus_name]
        records_dict = dict()
        
        log_file = output_dir + ps + locus_name + '.log'
        log_handle = open(log_file, 'w')

        for record in records:
            taxid = record.description.split('|')[4]
            if not records_dict.has_key(taxid):
                records_dict[taxid] = list()
            records_dict[taxid].append(record)
        records_count = len(records_dict.keys())
        
        krcl.hide_cursor()

        for i, taxid in enumerate(records_dict.keys()):
            tax_records = records_dict[taxid]
            krcl.print_progress(i+1, records_count, 50, '\t')
            # If there is more than one sequence for particular locus and
            # particular organism.
            if len(tax_records) > 1:
                # ...we cluster these sequences and hope for only one cluster
                cluster_dict = krusearch.cluster_records(tax_records, min_similarity,
                                               temp_dir)
                # ...if there is only one cluster
                if len(cluster_dict.keys()) == 1:
                    # ...we pick the longest available sequence
                    tax_records.sort(key=lambda x:len(x), reverse=True)
                    results.append(tax_records[0])
                else:
                    # ToDo: What happens if the sequences from the same
                    # organism and the same locus are not similar enough?
                    # Note: This happened only twice in my expirience so far,
                    # so I don't think this is top priority.
                    log_handle.write(tax_records[0].description.split('|')[-2] + '\t' + taxid + '\tSequences are too dissimilar.\n')
            else:
                results.append(tax_records[0])

        krcl.show_cursor()

        for result in results:
            desc = result.description.split('|')
            desc = desc[-2] + '|' + desc[-1]
            result.id = desc
            result.name = ''
            result.description = ''
        krbioio.write_sequence_file(results, output_file, 'fasta')

        log_handle.close()

        print('\n\tAccepted', len(results), 'sequences.')

    os.removedirs(temp_dir)

def align_loci(processed_results_dir, output_dir, program, threads, spacing,
    temp_dir):

    '''
    Align individual loci, then concatenate alignments.

    spacing - the number of gaps to insert between individual alignments.
    '''

    print('\nAlign individual loci, then concatenate alignments.')

    import os
    import krio
    import krbioio
    import kralign

    ps = os.path.sep

    alignments = []
    
    print('\tPreparing output directory "', output_dir, '"', sep='')
    krio.prepare_directory(output_dir)
    print('\tPreparing temporary directory "', temp_dir, '"', sep='')
    krio.prepare_directory(temp_dir)

    gene_dict = dict()

    file_list = parse_directory(processed_results_dir, ' ')
    
    for f in file_list:

        if not f['ext'].startswith('fasta'):
            continue
        
        file_name = f['name']
        output_file = output_dir + ps + file_name + '.fasta'
        records = krbioio.read_sequence_file(f['path'], 'fasta')

        # Align each locus individually first.
        print('\n\tAligning', file_name)
        aln = kralign.align(records, program, threads, temp_dir)
        krbioio.write_alignment_file(aln, output_file, 'fasta')
        alignments.append(aln)

    print('\n\tProducing concatenated alignment.')
    concatenated = kralign.concatenate(alignments, spacing)
    concatenated_output_file = output_dir + ps + 'concatenated' + '.fasta'
    krbioio.write_alignment_file(concatenated, concatenated_output_file, 'fasta')

    os.removedirs(temp_dir)

# End pipeline functions ------------------------------------------------------

if __name__ == '__main__':
    
    import sys
    import argparse
    import os
    import krio
    import krbioio
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--test', action='store_true', help='Run tests.')
    parser.add_argument('-c', '--command', type=unicode, choices=['search_and_download', 'extract_loci', 'one_locus_per_organism', 'align_loci'], help='Run a command.')
    parser.add_argument('-q', '--query', type=unicode, help='Query file path.')
    parser.add_argument('-o', '--output', type=unicode, help='Output directory path.')
    parser.add_argument('-s', '--sep', type=unicode, help='Output file name separator.')
    parser.add_argument('-e', '--email', type=unicode, help='Email.')
    parser.add_argument('-i', '--input', type=unicode, help='Input directory path.')
    parser.add_argument('--ncbinames', type=unicode, help='NCBI organism names file.')
    parser.add_argument('--synonymy', type=unicode, help='Synonymy file.')
    parser.add_argument('--authority', type=unicode, help='Authority alternates file.')
    parser.add_argument('--samples', type=unicode, help='Sequence samples file in FASTA format.')
    parser.add_argument('--similarity', type=float, help='Minimum sequence similarity.')
    parser.add_argument('--tempdir', type=unicode, help='Temporary directory path.')
    parser.add_argument('--alignprog', type=unicode, choices=['mafft', 'einsi', 'linsi', 'muscle'], help='Alignment program to be used.')
    parser.add_argument('--threads', type=int, help='Number of CPU cores to use.')
    parser.add_argument('--alignspacing', type=int, help='Number of gaps to add between alignments in concatenated alignment.')

    args = parser.parse_args()

    PS = os.path.sep

    if args.test:

        print('Running tests.')

        # Tests

        # parse_directory
        t_parse_directory = parse_directory('testdata'+PS+'parsedirectory',
            '$')
        for d in t_parse_directory:
            print(d)

        sys.exit(0)

    else:

        if args.command == 'search_and_download':

            queries = krio.read_table_file(args.query, has_headers=True,
                          headers=None, delimiter='\t')

            search_and_download(queries, args.output+PS, args.sep, args.email)

        if args.command == 'extract_loci':

            ncbi_names = None
            synonymy_table = None
            sequence_samples = None

            if args.ncbinames:

                ncbi_names = krio.read_table_file(args.ncbinames,
                    has_headers=False,
                    headers=('tax_id', 'name_txt', 'unique_name', 'name_class'),
                    delimiter='\t|')

            if args.authority and args.synonymy:

                synonymy_table = krio.read_table_file(args.synonymy,
                    has_headers=True, headers=None, delimiter=',')

            if args.samples:

                sequence_samples = krbioio.read_sequence_file(args.samples, 'fasta',
                    ret_type='dict')
            
            extract_loci(
                search_results_dir = args.input,
                output_dir = args.output,
                sequence_samples = sequence_samples,
                ncbi_names_table = ncbi_names,
                synonymy_table = synonymy_table,
                auth_file = args.authority,
                min_similarity = args.similarity,
                temp_dir = args.tempdir,
                file_name_sep = args.sep
                )

        if args.command == 'one_locus_per_organism':

            one_locus_per_organism(
                extracted_results_dir = args.input,
                output_dir = args.output,
                min_similarity = args.similarity,
                temp_dir = args.tempdir,
                file_name_sep = args.sep
                )

        if args.command == 'align_loci':

            align_loci(
                processed_results_dir = args.input,
                output_dir = args.output,
                program = args.alignprog,
                threads = args.threads,
                spacing = args.alignspacing,
                temp_dir = args.tempdir
                )
