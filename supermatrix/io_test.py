'''
    Created on Dec 27, 2013
    @author: Karolis Ramanauskas
    @copyright: 2013 Karolis Ramanauskas. All rights reserved.
    @license: GPLv3
'''

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

__updated__ = '2013-12-27'


def test_parse_search_queries_file():
    '''
    test_parse_search_queries_file
    '''
    from krpy import supermatrix
    with open('test_data/search_queries.tsv', 'r') as search_queries_handle:
        search_queries = supermatrix.io.parse_search_queries_file(
                                                        search_queries_handle)
    assert 'GBSSI' in search_queries.keys()
    assert 'matK' in search_queries.keys()

    for key in search_queries.keys():
        if key == 'GBSSI':
            assert search_queries['GBSSI'][0]['query'] == \
            '(waxy[Gene Name] OR GBSSI[Gene Name]) '\
            'AND (txid197382[Organism] OR txid197388[Organism] OR '\
            'txid4121[Organism] OR txid4120[Organism] OR txid4123[Organism] '\
            'OR txid4070[Organism])'
        if key == 'matK':
            for name_2_line in search_queries['matK']:
                if name_2_line['name2'] == '2':
                    assert name_2_line['query'] == 'DONOTDOWNLOAD'


def test_parse_input_file():
    '''
    test_parse_input_file
    '''
    from krpy import supermatrix
    out_dict = supermatrix.io.parse_input_file('test_data/sm_input.krsm')
    assert out_dict['search_queries'] is not None
