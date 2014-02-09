#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import unicode_literals

if __name__ == '__main__':

    import argparse
    from Bio import AlignIO
    from Bio.Align import MultipleSeqAlignment
    from Bio.SeqRecord import SeqRecord

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_file', type=unicode,
                        help='Input alignment file path.')
    parser.add_argument('-o', '--output_file', type=unicode,
                        help='Output alignment file path.')
    parser.add_argument('-f', '--format', type=unicode,
                        help='Alignment file format.')

    # fasta, phylip-relaxed
    # http://biopython.org/wiki/AlignIO

    input_file = None
    output_file = None
    format = None

    args = parser.parse_args()

    if args.input_file:
        input_file = args.input_file
    if args.output_file:
        output_file = args.output_file
    if args.format:
        format = args.format

    if input_file and output_file and format:

        alignment = AlignIO.read(input_file, format)
        alignment.sort()
        AlignIO.write(alignment, output_file, format)
