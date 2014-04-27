#!/usr/bin/env python

from __future__ import print_function
#from __future__ import unicode_literals

def search_genbank(ncbi_db, query_terms, ncbi_tax_ids, max_seq_length, email):

    from krpy import krncbi

    tax_ncbi_query_strings = list()
    for t in ncbi_tax_ids:
        tnqs = 'txid' + str(t) + '[Organism]'
        tax_ncbi_query_strings.append(tnqs)

    taxa_query_str = ' OR '.join(tax_ncbi_query_strings)
    taxa_query_str = '(' + taxa_query_str + ')'

    # query_term_str = '[Gene Name] OR '.join(query_terms)
    # query_term_str = query_term_str + '[Gene Name]'

    query_term_str = '[Title] OR '.join(query_terms)
    query_term_str = query_term_str + '[Title]'

    # query_term_str = ' OR '.join(query_terms)

    query_term_str = '(' + query_term_str + ')'
    seq_length_str = '0:' + str(max_seq_length) + '[Sequence Length]'

    query_str = taxa_query_str + ' AND ' + query_term_str + ' AND ' + \
                seq_length_str

    # print(query_str)

    result_uids = krncbi.esearch(esearch_terms=query_str, db=ncbi_db,
        email=krncbi)

    uid_list = list()

    for uid in result_uids:
        gi = int(uid)
        # print(gi, type(gi))
        uid_list.append(gi)

    return uid_list


def accept_records(records, temp_dir, min_clust_size=10):

    from krpy import krusearch

    clusters = krusearch.cluster_records(
        records=records,
        identity_threshold=0.55,
        temp_dir=temp_dir,
        sorted_input=False,
        algorithm='smallmem',  # fast smallmem
        strand='both',  # plus both aa
        threads=1,
        quiet=True,
        program='usearch',
        heuristics=True,
        query_coverage=0.2,
        target_coverage=0.2,
        sizein=False,
        sizeout=False,
        usersort=False,
        seq_id='gi',
        cluster_key='centroid'  # clust_number centroid
        )

    good = list()
    bad = list()

    for key in clusters.keys():
        clust_size = len(clusters[key])
        if clust_size >= min_clust_size:
            good = good + clusters[key]
        else:
            bad = bad + clusters[key]

        # print(key, clust_size)

    return (good, bad)


# def summarize_features(records):

#     from krpy import krncbi

#     min_taxa = 10

#     excluded_qualifiers = ['translation', 'db_xref', 'exception',
#     'rpt_unit_seq', 'gene_synonym', 'rpt_type', 'satellite', 'transl_table',
#     'replace', 'rpt_unit_range', 'protein_id', 'codon_recognized',
#     'EC_number', 'function', 'estimated_length', 'mobile_element_type',
#     'codon_start', 'transl_except', 'number', 'standard_name', 'allele',
#     'inference']

#     feature_dict = dict()
#     taxa_dict = dict()

#     for record in records:
#         txid = krncbi.get_ncbi_tax_id(record)
#         for feature in record.features:
#             if feature.type != 'source':
#                 for qualifier in feature.qualifiers:
#                     if qualifier not in excluded_qualifiers:
#                         key = feature.type + '.' + qualifier
#                         qualifier_label = feature.qualifiers[qualifier][0]
#                         qualifier_label_key = key + '.' + qualifier_label
#                         if key not in feature_dict.keys():
#                             feature_dict[key] = list()
#                         if qualifier_label_key not in taxa_dict.keys():
#                             taxa_dict[qualifier_label_key] = list()
#                         feature_dict[key].append(qualifier_label)
#                         taxa_dict[qualifier_label_key].append(txid)
#                         # print(key)
#                         # print(feature_dict[key])
#                         # print(qualifier_label_key)
#                         # print(taxa_dict[qualifier_label_key])
#                         # raw_input("Press Enter to continue...")

#     for key in feature_dict.keys():

#         print(key)

#         for q_label in set(feature_dict[key]):
#             qualifier_label_key = key + '.' + q_label
#             taxa = taxa_dict[qualifier_label_key]
#             taxa_count = len(set(taxa))
#             count = feature_dict[key].count(q_label)
#             if taxa_count >= min_taxa:
#                 print('\t' + q_label + ' : ' + str(count) + ' : ' + str(taxa_count))


# def generate_search_terms(ncbi_db, query_terms, ncbi_tax_ids, max_seq_length, email, temp_dir):

#     import os

#     from krpy import krncbi
#     from krpy import krbioio
#     from krpy import krusearch

#     ps = os.path.sep

#     temp_dir = temp_dir.rstrip(ps) + ps

#     gb_file_path = temp_dir + 'gst.gb'

#     gis = search_genbank(ncbi_db, query_terms, ncbi_tax_ids, max_seq_length, email)

#     krncbi.download_sequence_records(
#         file_path=gb_file_path,
#         uids=gis,
#         db=ncbi_db,
#         entrez_email=email)

#     records = krbioio.read_sequence_file(
#         file_path=gb_file_path,
#         file_format='gb',
#         ret_type='list',
#         key='gi')

#     clusters = krusearch.cluster_records(
#         records=records,
#         identity_threshold=0.55,
#         temp_dir=temp_dir,
#         sorted_input=False,
#         algorithm='smallmem',  # fast smallmem
#         strand='both',  # plus both aa
#         threads=1,
#         quiet=True,
#         program='usearch',
#         heuristics=True,
#         query_coverage=0.2,
#         target_coverage=0.2,
#         sizein=False,
#         sizeout=False,
#         usersort=False,
#         seq_id='gi',
#         cluster_key='centroid'  # clust_number centroid
#         )

#     for key in clusters.keys():
#         print(key, len(clusters[key]))
