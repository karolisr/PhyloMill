#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
# from __future__ import unicode_literals

if __name__ == '__main__':

    import argparse

    from krpy import krbioio
    from krpy import krncbi

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_file', type=unicode,
                        help='Input alignment file path.')
    parser.add_argument('-o', '--output_file', type=unicode,
                        help='Output alignment file path.')
    parser.add_argument('-m', '--min_taxa', type=int,
                        help='Minimum number of taxa represented.')

    # fasta, phylip-relaxed
    # http://biopython.org/wiki/AlignIO

    input_file = None
    output_file = None
    min_taxa = 0

    args = parser.parse_args()

    if args.input_file:
        input_file = args.input_file
    if args.output_file:
        output_file = args.output_file
    if args.min_taxa:
        min_taxa = args.min_taxa

    if input_file and output_file and min_taxa:
        records = krbioio.read_sequence_file(input_file, 'gb', ret_type='list')

        # print('Found', len(records), 'in', input_file)

        excluded_qualifiers = ['translation', 'db_xref', 'exception',
        'rpt_unit_seq', 'gene_synonym', 'rpt_type', 'satellite', 'transl_table',
        'replace', 'rpt_unit_range', 'protein_id', 'codon_recognized',
        'EC_number', 'function', 'estimated_length', 'mobile_element_type',
        'codon_start', 'transl_except', 'number', 'standard_name', 'allele',
        'inference']

        feature_dict = dict()
        taxa_dict = dict()

        for record in records:
            txid = krncbi.get_ncbi_tax_id(record)
            for feature in record.features:
                if feature.type != 'source':
                    for qualifier in feature.qualifiers:
                        if qualifier not in excluded_qualifiers:
                            key = feature.type + '.' + qualifier
                            qualifier_label = feature.qualifiers[qualifier][0]
                            qualifier_label_key = key + '.' + qualifier_label
                            if key not in feature_dict.keys():
                                feature_dict[key] = list()
                            if qualifier_label_key not in taxa_dict.keys():
                                taxa_dict[qualifier_label_key] = list()
                            feature_dict[key].append(qualifier_label)
                            taxa_dict[qualifier_label_key].append(txid)
                            # print(key)
                            # print(feature_dict[key])
                            # print(qualifier_label_key)
                            # print(taxa_dict[qualifier_label_key])
                            # raw_input("Press Enter to continue...")

        with open(output_file, 'w') as output:
            for key in feature_dict.keys():

                print(key)
                output.write(key + '\n')

                for q_label in set(feature_dict[key]):
                    qualifier_label_key = key + '.' + q_label
                    taxa = taxa_dict[qualifier_label_key]
                    taxa_count = len(set(taxa))
                    count = feature_dict[key].count(q_label)
                    if taxa_count >= min_taxa:
                        print('\t' + q_label + ' : ' + str(count) + ' : ' + str(taxa_count))
                        output.write('\t' + q_label + ' : ' + str(count) + ' : ' + str(taxa_count) + '\n')

                # value_dict = feature_dict[key]
                # value = set(value_dict['q'])
                # taxa_count = len(value_dict['t'])
                # count = value_dict['q'].count(value)
                # if count >= min_occur:
                #     print('\t' + value + ' : ' + str(count) + ' : ' + str(taxa_count))
                #     output.write('\t' + value + ' : ' + str(count) + ' : ' + str(taxa_count))

