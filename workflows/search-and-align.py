#!/usr/bin/env python

from __future__ import print_function

'''
This pipeline expects these files to be located in the input directory:

    config
'''

if __name__ == '__main__':

    import os
    import sys

    import ConfigParser

    import argparse

    import krseqsearch
    import krio
    import krbioio

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_dir', type=unicode,
                        help='Input directory path.')
    parser.add_argument('-o', '--output_dir', type=unicode,
                        help='Output directory path.')
    parser.add_argument('-c', '--commands', type=unicode,
                        help='Commands to run.')

    args = parser.parse_args()

    input_dir = None
    output_dir = None

    ps = os.path.sep

    if args.input_dir:
        input_dir = args.input_dir.rstrip(ps) + ps
    if args.output_dir:
        output_dir = args.output_dir.rstrip(ps) + ps

    if (not input_dir) or (not output_dir):
        print('Input and output directories are required.')
        sys.exit(1)
    else:

        file_name_sep = '$'

        query_file = input_dir + 'search_queries'
        synonymy_file = input_dir + 'synonymy.csv'
        ncbi_names_file = input_dir + 'ncbi_tax_names'
        authority_file = input_dir + 'authority_alternates.dat'
        cutlist_records_file = input_dir + 'cutlist_records'
        cutlist_records_auto_file = input_dir + 'cutlist_records_auto'
        # keeplist_records_file = input_dir + 'keeplist_records'
        cutlist_taxonomy_file = input_dir + 'cutlist_taxonomy'
        unresolvable_taxonomy_file = input_dir + 'unresolvable.csv'
        keeplist_taxonomy_file = input_dir + 'keeplist_taxonomy'
        force_reverse_complement_file = input_dir + 'force_reverse_complement_records'
        good_sequences_dir = input_dir + 'good_sequences'
        additional_sequences_file = input_dir + 'additional_sequences.fasta'
        taxa_mappings_file = input_dir + 'taxa_mappings.tsv'

        temp_dir = output_dir + 'temp'
        log_dir = output_dir + '99-logs'

        search_results_dir = output_dir + '01-search-results'
        filtered_results_dir = output_dir + '02-filtered-results'
        extract_loci_dir = output_dir + '03-extracted-loci'
        one_locus_per_organism_dir = output_dir + '04-one-locus-per-organism'
        aligned_dir = output_dir + '05-aligned'
        concatenated_dir = output_dir + '06-concatenated'

        # Create output directory
        krio.prepare_directory(output_dir)
        # Create log directory
        krio.prepare_directory(log_dir)

        # Get configuration information
        config_file_path = input_dir + 'config'
        config = ConfigParser.SafeConfigParser()
        config.read(config_file_path)

        hacks_items = config.items('Hacks')
        hacks = list()
        for hack in hacks_items:
            h = hack[0].split('_')
            if h[-1] != 'data':
                if config.getboolean('Hacks', '_'.join(h)):
                    hacks.append('_'.join(h))

        hacks_data_location = dict()

        for hack in hacks:
            hacks_data_location[hack] = config.get('Hacks', hack + '_data')

        # Parse queries file
        queries = krio.read_table_file(path=query_file, has_headers=True,
                                       headers=None, delimiter='\t',
                                       quotechar="'")

        commands = None
        if args.commands:
            commands = set([x.strip() for x in args.commands.split(',')])

        if not commands or 'search' in commands:

            # Search and download
            krseqsearch.search_and_download(
                queries=queries,
                output_dir=search_results_dir,
                file_name_sep=file_name_sep,
                email=config.get('General', 'email'),
                input_dir=input_dir,
                log_dir=log_dir)

        if not commands or 'filter' in commands:

            # Filter records
            krseqsearch.filter_results(
                search_results_dir=search_results_dir,
                output_dir=filtered_results_dir,
                cutlist_records_file=cutlist_records_file,
                cutlist_taxonomy_file=cutlist_taxonomy_file
            )

        if not commands or 'extract' in commands:

            ncbi_names = krio.read_table_file(
                path=ncbi_names_file,
                has_headers=False,
                headers=('tax_id', 'name_txt', 'unique_name', 'name_class'),
                delimiter='\t|',
                quotechar=None,
                stripchar='"',
                rettype='dict')

            synonymy_table = krio.read_table_file(
                path=synonymy_file,
                has_headers=True, headers=None, delimiter=',')

            unresolvable_taxonomy_list = krio.read_table_file(
                path=unresolvable_taxonomy_file,
                has_headers=True,
                headers=None,
                delimiter=',',
                quotechar=None,
                stripchar='"',
                rettype='dict')

            keeplist_taxonomy_list = krio.read_table_file(
                path=keeplist_taxonomy_file,
                has_headers=False,
                headers=None,
                delimiter=',',
                quotechar=None,
                stripchar='"',
                rettype='set')

            force_reverse_complement_list = krio.read_table_file(
                path=force_reverse_complement_file,
                has_headers=False,
                headers=None,
                delimiter=',',
                quotechar=None,
                stripchar='"',
                rettype='set')

            taxa_mappings_list = krio.read_table_file(
                path=taxa_mappings_file,
                has_headers=False,
                headers=('accession', 'taxon'),
                delimiter='\t',
                quotechar=None,
                stripchar='"',
                rettype='dict')

            # Extract loci
            krseqsearch.extract_loci(
                search_results_dir=filtered_results_dir,
                output_dir=extract_loci_dir,
                queries=queries,
                ncbi_names_table=ncbi_names,
                temp_dir=temp_dir,
                file_name_sep=file_name_sep,
                synonymy_table=synonymy_table,
                auth_file=authority_file,
                hacks=hacks,
                hacks_data_location=hacks_data_location,
                unresolvable_taxonomy_list=unresolvable_taxonomy_list,
                keeplist_taxonomy_list=keeplist_taxonomy_list,
                force_reverse_complement_list=force_reverse_complement_list,
                taxa_mappings_list=taxa_mappings_list,
                log_dir=log_dir)

        if not commands or 'flatten' in commands:

            ref_aln_program = config.get('Flatten', 'ref_align_program')
            ref_aln_program_exe = config.get('General', ref_aln_program + '_executable')

            locus_aln_program = config.get('Flatten', 'locus_align_program')
            locus_aln_program_exe = config.get('General', locus_aln_program + '_executable')

            cutlist_records = krio.read_table_file(
                path=cutlist_records_file,
                has_headers=False,
                headers=None,
                delimiter=',',
                quotechar='"',
                rettype='set')

            cutlist_records_auto = krio.read_table_file(
                path=cutlist_records_auto_file,
                has_headers=False,
                headers=None,
                delimiter=',',
                quotechar='"',
                rettype='set')

            # keeplist_records = krio.read_table_file(
            #     path=keeplist_records_file,
            #     has_headers=False,
            #     headers=None,
            #     delimiter=',',
            #     quotechar='"',
            #     rettype='set')

            additional_sequences = krbioio.read_sequence_file(additional_sequences_file, 'fasta')

            #print(additional_sequences)

            #if len(additional_sequences) == 0:
            #    additional_sequences = list()

            # One locus per organism
            new_cutlist_records_auto = krseqsearch.one_locus_per_organism(
                extracted_results_dir=extract_loci_dir,
                output_dir=one_locus_per_organism_dir,
                queries=queries,
                cutlist_records=cutlist_records,
                cutlist_records_auto=cutlist_records_auto,
                # keeplist_records=keeplist_records,
                temp_dir=temp_dir,
                file_name_sep=file_name_sep,
                usearch_exe=config.get('General', 'usearch_executable'),
                ref_aln_program_exe=ref_aln_program_exe,
                ref_aln_program=ref_aln_program,
                ref_aln_program_options=config.get('Flatten', 'ref_align_program_options'),
                locus_aln_program_exe=locus_aln_program_exe,
                locus_aln_program=locus_aln_program,
                locus_aln_program_options=config.get('Flatten', 'locus_align_program_options'),
                good_sequences_dir=good_sequences_dir,
                log_dir=log_dir,
                additional_sequences=additional_sequences)

            f = open(cutlist_records_auto_file, 'wb')
            for clr in new_cutlist_records_auto:
                f.write(clr + '\n')
            f.close()

        if not commands or 'align' in commands:

            aln_program = config.get('Align', 'align_program')
            aln_program_exe = config.get('General', aln_program + '_executable')
            aln_options = config.get('Align', 'align_program_options')

            # Align loci
            krseqsearch.align_loci(
                processed_results_dir=one_locus_per_organism_dir,
                output_dir=aligned_dir,
                aln_program_exe=aln_program_exe,
                aln_program=aln_program,
                aln_options=aln_options,
                temp_dir=temp_dir)

        if not commands or 'concatenate' in commands:

            order_of_loci = config.get('Concatenate', 'order_of_loci').split(',')
            number_of_gaps_between_loci = config.getint('Concatenate', 'number_of_gaps_between_loci')

            # Produce a concatenated alignment
            krseqsearch.concatenate(
                aligned_loci_dir=aligned_dir,
                output_dir=concatenated_dir,
                order_of_loci=order_of_loci,
                number_of_gaps_between_loci=number_of_gaps_between_loci,
                log_dir=log_dir)
