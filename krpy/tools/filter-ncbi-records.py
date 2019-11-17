#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import unicode_literals

if __name__ == '__main__':

    import argparse

    from krpy import krseqsearch
    from krpy import krio
    from krpy import krbioio

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_file', type=unicode,
                        help='Records to download, one per line.')
    parser.add_argument('-o', '--output_file', type=unicode,
                        help='Results output file.')
    parser.add_argument('-r', '--cutlist_records_file', type=unicode,
                        help='')
    parser.add_argument('-t', '--cutlist_taxonomy_file', type=unicode,
                        help='')

    input_file = None
    output_file = None
    cutlist_records_file = None
    cutlist_taxonomy_file = None

    args = parser.parse_args()

    if args.input_file:
        input_file = args.input_file
    if args.output_file:
        output_file = args.output_file
    if args.cutlist_records_file:
        cutlist_records_file = args.cutlist_records_file
    if args.cutlist_taxonomy_file:
        cutlist_taxonomy_file = args.cutlist_taxonomy_file

    if input_file and output_file and (cutlist_records_file or cutlist_taxonomy_file):
        records = krbioio.read_sequence_file(input_file, 'gb', ret_type='list')
        filtered = krseqsearch.filter_records(records, cutlist_records_file, cutlist_taxonomy_file)
        krbioio.write_sequence_file(filtered, output_file, 'gb')
