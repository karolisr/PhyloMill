# -*- coding: utf-8 -*-

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
__updated__ = '2013-12-28'

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
    search_queries_dict = dict()
    for query_line in search_queries_raw:
        if query_line['name1'] not in search_queries_dict:
            search_queries_dict[query_line['name1']] = list()
        search_queries_dict[query_line['name1']].append(query_line)

    # Construct hierarchical structure.
    search_queries_tree = krpy.Node('loci')
    for key in search_queries_dict.keys():
        name_1_root = krpy.Node(name=key, data=None,
                                parent=search_queries_tree)
        for key_2_entry in search_queries_dict[key]:
            krpy.Node(
                name=key_2_entry['name2'],
                data=key_2_entry,
                parent=name_1_root)

    if krpy.debug.RUN_DEBUG_CODE:
        krpy.debug.message('Search queries:', parse_search_queries_file)
        print(search_queries_tree)

    return search_queries_tree


def write_search_queries_file(search_queries_tree, file_handle):
    '''
    Write tab-separated search queries file.
    '''
    headers = None
    for name1_node in search_queries_tree.children():
        for name2_node in name1_node.children():
            if not headers:
                headers = sorted(name2_node.data().keys())
                file_handle.write('\t'.join(headers) + '\n')
            fields = list()
            for header in headers:
                fields.append(name2_node.data()[header])
            file_handle.write('\t'.join(fields) + '\n')


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
            krpy.debug.message('Listing file: ' + file_name + ' in ' +
                               file_path + '.', parse_input_file)
            if file_name == search_queries_file_name:
                search_queries_handle = zip_file.open(file_name, 'rU')

            # TODO: check that all files exist in the input file.

    search_queries = parse_search_queries_file(search_queries_handle)

    search_queries_handle.close()

    return {'search_queries': search_queries}

if __name__ == '__main__':
    if krpy.debug.RUN_DEBUG_CODE:
        INPUT_FILE_DICT = parse_input_file('test_data/sm_input.krsm')
        with open('test_data/search_queries_write_test.tsv', 'w') as \
        FILE_HANDLE:
            write_search_queries_file(INPUT_FILE_DICT['search_queries'],
                                      FILE_HANDLE)
