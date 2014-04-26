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

    query_term_str = ' OR '.join(query_terms)
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
