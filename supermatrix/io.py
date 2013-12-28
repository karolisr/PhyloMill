'''
Created on Dec 26, 2013
@author: Karolis Ramanauskas
@copyright: 2013 Karolis Ramanauskas. All rights reserved.
@license: GPLv3
'''

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

__all__ = []
__version__ = 0.1
__updated__ = '2013-12-27'

import zipfile
import krpy


def parse_search_queries_file(search_queries_handle):
    '''
    Parse tab-separated search queries file.
    TODO: describe search queries file structure.
    '''
    search_queries_raw = krpy.io.read_table_file(handle=search_queries_handle,
        has_headers=True, headers=None, delimiter='\t', quotechar="'")

    # Produce a dictionary with name1 as keys and a list of dictionaries as
    # value. Each item in a list will have a unique name2.
    search_queries = dict()
    for query_line in search_queries_raw:
        if query_line['name1'] not in search_queries:
            search_queries[query_line['name1']] = list()
        search_queries[query_line['name1']].append(query_line)

    if krpy.debug.RUN_DEBUG_CODE:
        for key in search_queries.keys():
            for key_2_entry in search_queries[key]:
                krpy.debug.message(
                    'name1: ' + key + ', name2: ' + key_2_entry['name2'],
                    parse_search_queries_file)

    return search_queries


def parse_input_file(file_path):
    '''
    Read an input (zip) file and return correct representations of all input
    files.
    '''
    search_queries_file_name = 'search_queries.tsv'

    search_queries_handle = None

    with zipfile.ZipFile(file_path, 'r') as zip_file:
        krpy.debug.message('Reading file: ' + file_path + '.',
                           parse_input_file)
        for file_name in zip_file.namelist():
            krpy.debug.message('Reading file: ' + file_name + ' in ' +
                               file_path + '.', parse_input_file)
            if file_name == search_queries_file_name:
                search_queries_handle = zip_file.open(file_name, 'rU')

            # TODO: check that all files exist in the input file.

    search_queries = parse_search_queries_file(search_queries_handle)

    search_queries_handle.close()

    return {'search_queries': search_queries}

if __name__ == '__main__':
    if krpy.debug.RUN_DEBUG_CODE:
        parse_input_file('test_data/sm_input.krsm')
