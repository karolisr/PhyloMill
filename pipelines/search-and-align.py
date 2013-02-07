#!/usr/bin/env python

from __future__ import print_function

'''
This pipeline expects these files to be located in the input directory:

    config
'''

if __name__ == '__main__':

    import os
    import sys
    import argparse
    import krpipe
    import krio
    import krbioio

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_dir', type=unicode,
        help='Input directory path.')
    parser.add_argument('-o', '--output_dir', type=unicode,
        help='Output directory path.')
    parser.add_argument('--steps', type=unicode,
        help='Which steps to perform 1, 2, 3, 4')

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

        query_file = input_dir + 'search-queries'
        sequence_samples_file = input_dir + 'sequence-samples'
        synonymy_file = input_dir + 'synonymy'
        ncbi_names_file = input_dir + 'ncbi-names'
        authority_file = input_dir + 'authority'

        temp_dir = output_dir + 'temp'

        search_results_dir = output_dir + '01-search-results'
        extract_loci_dir = output_dir + '02-extracted-loci'
        one_locus_per_organism_dir = output_dir + '03-one-locus-per-organism'
        aligned_dir = output_dir + '04-aligned'

        # Create output directory
        krio.prepare_directory(output_dir)

        # Get configuration information
        config = None
        config_dict = dict()
        config_file_path = input_dir + 'config'
        if os.path.exists(config_file_path):
            config = krio.read_table_file(config_file_path, has_headers=False,
                headers=['name', 'value'], delimiter='\t', quotechar='"')
            for l in config:
                config_dict[l['name']] = l['value']

        # Parse queries file
        queries = krio.read_table_file(query_file, has_headers=True,
            headers=None, delimiter='\t')

        steps = None

        if args.steps:
            steps = set([int(x.strip()) for x in args.steps.split(',')])

        if not steps or 1 in steps:

            # Search and download
            krpipe.search_and_download(
                queries=queries,
                output_dir=search_results_dir,
                file_name_sep=file_name_sep,
                email=config_dict['email'])

        if not steps or 2 in steps:

            ncbi_names = krio.read_table_file(ncbi_names_file,
                has_headers=False,
                headers=('tax_id', 'name_txt', 'unique_name', 'name_class'),
                delimiter='\t|')
            synonymy_table = krio.read_table_file(synonymy_file,
                has_headers=True, headers=None, delimiter=',')
            sequence_samples = krbioio.read_sequence_file(
                sequence_samples_file, 'fasta', ret_type='dict')

            # Extract loci
            krpipe.extract_loci(
                search_results_dir=search_results_dir,
                output_dir=extract_loci_dir,
                sequence_samples=sequence_samples,
                ncbi_names_table=ncbi_names,
                min_similarity=config_dict['asim'],
                temp_dir=temp_dir,
                file_name_sep=file_name_sep,
                synonymy_table=synonymy_table,
                auth_file=authority_file)

        if not steps or 3 in steps:

            # One locus per organism
            krpipe.one_locus_per_organism(
                extracted_results_dir=extract_loci_dir,
                output_dir=one_locus_per_organism_dir,
                min_similarity=config_dict['wsim'],
                temp_dir=temp_dir,
                file_name_sep=file_name_sep)

        if not steps or 4 in steps:

            # Align loci and produce a concatenated alignment
            krpipe.align_loci(
                processed_results_dir=one_locus_per_organism_dir,
                output_dir=aligned_dir,
                program=config_dict['alnprg'],
                threads=config_dict['threads'],
                spacing=config_dict['alngaps'],
                temp_dir=temp_dir)
