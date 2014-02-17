#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
# from __future__ import unicode_literals

if __name__ == '__main__':

    import os
    import argparse

    from krpy import krseqsearch
    from krpy import krio
    from krpy import krbioio
    from krpy import krseq
    from krpy import krbionames
    from krpy import krcl

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_file', type=unicode,
                        help='')
    parser.add_argument('-o', '--output_file', type=unicode,
                        help='')
    parser.add_argument('-l', '--log_dir', type=unicode,
                        help='')

    parser.add_argument('-n', '--ncbi_names_file', type=unicode,
                        help='')
    parser.add_argument('-s', '--synonymy_file', type=unicode,
                        help='')
    parser.add_argument('-u', '--unresolvable_taxonomy_file', type=unicode,
                        help='')
    parser.add_argument('-k', '--keeplist_taxonomy_file', type=unicode,
                        help='')
    parser.add_argument('-t', '--taxa_mappings_file', type=unicode,
                        help='')
    parser.add_argument('-a', '--authority_file', type=unicode,
                        help='')
    parser.add_argument('-c', '--hacks', type=unicode,
                        help='')
    parser.add_argument('-d', '--hacks_data_location', type=unicode,
                        help='')

# record, ncbi_names_table, synonymy_table, auth_file,
# hacks, hacks_data_location, unresolvable_taxonomy_list,
# keeplist_taxonomy_list, taxa_mappings_list, log_dir

    args = parser.parse_args()

    hacks = None
    hacks_data_location = None

    if args.hacks:
        hacks = args.hacks.split(',')
        if args.hacks_data_location:
            hacks_data_location = dict()
            for i, hack in enumerate(hacks):
                hacks_data_location[hack] = args.hacks_data_location.split(',')[i]

    ncbi_names_table = None
    if args.ncbi_names_file:
        ncbi_names_table = krio.read_table_file(
            path=args.ncbi_names_file,
            has_headers=False,
            headers=('tax_id', 'name_txt', 'unique_name', 'name_class'),
            delimiter='\t|',
            quotechar=None,
            stripchar='"',
            rettype='dict')

    synonymy_table = None
    if args.synonymy_file:
        synonymy_table = krio.read_table_file(
            path=args.synonymy_file,
            has_headers=True, headers=None, delimiter=',')

    unresolvable_taxonomy_list = None
    if args.unresolvable_taxonomy_file:
        unresolvable_taxonomy_list = krio.read_table_file(
            path=args.unresolvable_taxonomy_file,
            has_headers=True,
            headers=None,
            delimiter=',',
            quotechar=None,
            stripchar='"',
            rettype='dict')

    keeplist_taxonomy_list = None
    if args.keeplist_taxonomy_file:
        keeplist_taxonomy_list = krio.read_table_file(
            path=args.keeplist_taxonomy_file,
            has_headers=False,
            headers=None,
            delimiter=',',
            quotechar=None,
            stripchar='"',
            rettype='set')

    taxa_mappings_list = None
    if args.taxa_mappings_file:
        taxa_mappings_list = krio.read_table_file(
            path=args.taxa_mappings_file,
            has_headers=False,
            headers=('accession', 'taxon'),
            delimiter='\t',
            quotechar=None,
            stripchar='"',
            rettype='dict')

    input_file = None
    output_file = None
    authority_file = None
    log_dir = None

    if args.input_file:
        input_file = args.input_file
    if args.output_file:
        output_file = args.output_file
    if args.authority_file:
        authority_file = args.authority_file
    if args.log_dir:
        log_dir = args.log_dir

    records = krbioio.read_sequence_file(input_file, 'gb', ret_type='list')

    ps = os.path.sep
    tax_log_handle = krseqsearch.__tax_log_open(log_dir, ps)
    tax_log_html_handle = krseqsearch.__tax_log_html_open(log_dir, ps)

    #########
    krcl.hide_cursor()

    for i, record in enumerate(records):

        krcl.print_progress(i, len(records), 50, '')

        name = krseqsearch.check_organism_name(
            record,
            ncbi_names_table,
            synonymy_table,
            authority_file,
            hacks,
            hacks_data_location,
            unresolvable_taxonomy_list,
            keeplist_taxonomy_list,
            taxa_mappings_list,
            tax_log_handle,
            tax_log_html_handle)

        # tn = name[0]
        an = name[1]

        an_flat = krbionames.flatten_organism_name(an, ' ')

        record.annotations['organism_old'] = record.annotations['organism']
        record.annotations['organism'] = an_flat
        record.annotations['source'] = an_flat

        record.description = record.description.replace(record.annotations['organism'], '')
        record.description = record.description.strip()

    krcl.show_cursor()
    #########

    krseqsearch.__tax_log_close(tax_log_handle)
    krseqsearch.__tax_log_html_close(tax_log_html_handle)

    krbioio.write_sequence_file(records, output_file, 'gb')
