#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import unicode_literals

if __name__ == '__main__':

    import argparse

    from krpy import krncbi
    from krpy import krio

    parser = argparse.ArgumentParser()

    parser.add_argument('-d', '--db', type=unicode,
                        help='Database.')
    parser.add_argument('-o', '--output_file', type=unicode,
                        help='Results output file.')
    parser.add_argument('-e', '--email', type=unicode,
                        help='Email.')
    parser.add_argument('-i', '--input_file', type=unicode,
                        help='Records to download, one per line.')
    parser.add_argument('-r', '--records', type=unicode,
                        help='Records to download, separated by commas.')

    input_file = None
    records = None
    output_file = None

    args = parser.parse_args()

    if args.input_file:
        input_file = args.input_file
    if args.records:
        records = args.records
    if args.output_file:
        output_file = args.output_file
    if args.email:
        email = args.email
    if args.db:
        db = args.db

    if (input_file or records) and output_file and email and db:
        if records:
            records = records.split(',')
        elif input_file:
            records = krio.read_table_file(
                handle=None,
                path=input_file,
                has_headers=False,
                headers=None,
                delimiter=',',
                quotechar='"',
                stripchar='',
                commentchar="#",
                rettype='set'  # dict, list, set
            )

        records = [int(x) for x in records]
        records = sorted(list(records))
        records = [str(x) for x in records]

        krncbi.download_sequence_records(output_file, records, db, email)
