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
                        help='Newick output file path.')
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

    genus_delimiter = '_'

    if input_file and output_file and format:

        alignment = AlignIO.read(input_file, format)

        genus_dict = dict()

        for a in alignment:
            full_name = a.id
            genus = full_name.split(genus_delimiter)[0]
            if genus not in genus_dict.keys():
                genus_dict[genus] = list()
            genus_dict[genus].append(full_name)

        handle = open(output_file, 'w')
        handle.write('(')
        first_genus = True
        for genus in genus_dict.keys():
            if first_genus:
                handle.write('(')
                first_genus = False
            else:
                handle.write(',')
            first_species = True
            number_of_tips = len(genus_dict[genus])
            handle.write('('*(min(1, number_of_tips-1)))
            for i,full_name in enumerate(genus_dict[genus]):
                # if i != 0 and i%2 == 0:
                #     handle.write(')')
                if not first_species:
                    handle.write(',')


                # else:
                handle.write(''+full_name+'')
                if not first_species or number_of_tips == 1:
                    handle.write(')')
                first_species = False
            # if i%2 == 0:
            #     handle.write(')')
        handle.write('));')
        handle.close()

