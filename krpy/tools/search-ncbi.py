#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import unicode_literals

if __name__ == '__main__':

    import argparse

    from krpy import krncbi

    parser = argparse.ArgumentParser()

    parser.add_argument('-q', '--query', type=unicode,
                        help='Search query.')
    parser.add_argument('-d', '--db', type=unicode,
                        help='Database.')
    parser.add_argument('-o', '--output_file', type=unicode,
                        help='Results output file.')
    parser.add_argument('-e', '--email', type=unicode,
                        help='Email.')

    input_file = None
    output_file = None

    args = parser.parse_args()

    if args.query:
        query = args.query
    if args.output_file:
        output_file = args.output_file
    if args.email:
        email = args.email
    if args.db:
        db = args.db

    if query and output_file and email and db:
        results = krncbi.esearch(query, db, email)
        results = [int(x) for x in results]
        results = sorted(list(results))
        handle = open(output_file, 'w')
        for uid in results:
            handle.write(str(uid)+'\n')
        handle.close()
