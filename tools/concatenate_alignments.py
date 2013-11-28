#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import unicode_literals

def concatenate(
    aligned_loci_dir,
    output_dir,
    order_of_loci,
    number_of_gaps_between_loci,
    log_dir
):

    '''
    Concatenate the alignments of individual loci.

    number_of_gaps_between_loci - the number of gaps to insert between individual alignments.
    '''

    print('\nConcatenate the alignments of individual loci.')

    import os

    from Bio import AlignIO

    import krio
    import krbioio
    import kralign
    import copy

    ps = os.path.sep

    print('\n\tPreparing output directory "', output_dir, '"', sep='')
    krio.prepare_directory(output_dir)

    file_list = krio.parse_directory(aligned_loci_dir, ' ')

    order_list = [x.strip() for x in order_of_loci]
    alignments = [x.strip() for x in order_of_loci]

    for f in file_list:

        if not f['ext'].startswith('phy'):
            continue

        file_name = f['name']

        aln = AlignIO.read(f['path'], "phylip-relaxed")
        if aln:
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
        matrix_output_file = log_dir + ps + '06-locus-presence' + '.csv'
        f = open(matrix_output_file, 'wb')
        f.write('taxon' + ',' + 'count' + ',' + ','.join(order_list) + '\n')
        f.write('' + ',' + '' + ',' + ','.join(length_list) + '\n')
        for key in matrix.keys():
            f.write(key + ',' + str(matrix[key].count('1')) + ',' +
                    ','.join(matrix[key]) + '\n')
        f.close()

        # Concatenate
        partitions_output_file = log_dir + ps + '06-locus-partitions' + '.csv'
        raxml_partitions_output_file = log_dir + ps + '06-locus-partitions-raxml'
        f_part = open(partitions_output_file, 'wb')
        f_part_raxml = open(raxml_partitions_output_file, 'wb')
        raw_alignments = list()
        for a in alignments:
            raw_alignments.append(a[0])
        concatenated = kralign.concatenate(raw_alignments, int(number_of_gaps_between_loci))
        cat_aln = concatenated[0]
        cat_partitions = concatenated[1]
        f_part.write('locus,start,end\n')
        for i, part in enumerate(cat_partitions):
            raxml_part_line = 'DNA, ' + order_list[i] + ' = ' + str(part[0]) + '-' + str(part[1]) + '\n'
            f_part_raxml.write(raxml_part_line)
            part_line = order_list[i] + ',' + str(part[0]) + ',' + str(part[1]) + '\n'
            f_part.write(part_line)
        concatenated_output_file = output_dir + ps + 'concatenated' + '.phy'
        krbioio.write_alignment_file(cat_aln, concatenated_output_file,
                                     'phylip-relaxed')
        f_part.close()
        f_part_raxml.close()

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', type=unicode,
                        help='')
    parser.add_argument('-o', '--output_dir', type=unicode,
                        help='')
    parser.add_argument('--order', type=unicode,
                        help='')
    parser.add_argument('--gaps', type=int,
                        help='')

    args = parser.parse_args()

    order = args.order.split(',')

    concatenate(args.input,args.output_dir,order,args.gaps,args.output_dir)
