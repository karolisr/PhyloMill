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
__updated__ = '2013-12-29'

import os
import zipfile
import tempfile
import shutil
import krpy

SEARCH_QUERIES_FILE_NAME = 'search_queries.tsv'
SEARCH_QUERIES_KEY = 'search_queries'


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


def write_search_queries_file(search_queries_tree, file_path):
    '''
    Write tab-separated search queries file.
    '''
    with open(file_path, 'w') as file_handle:
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

    search_queries_handle = None

    with zipfile.ZipFile(file_path, 'r') as zip_file:
        krpy.debug.message('Reading file: ' + file_path + '.',
                           parse_input_file)
        for file_name in zip_file.namelist():
            krpy.debug.message('Listing file: ' + file_name + ' in ' +
                               file_path + '.', parse_input_file)
            if file_name == SEARCH_QUERIES_FILE_NAME:
                search_queries_handle = zip_file.open(file_name, 'rU')

            # TODO: check that all files exist in the input file.

    search_queries_tree = parse_search_queries_file(search_queries_handle)

    search_queries_handle.close()

    return {SEARCH_QUERIES_KEY: search_queries_tree}


def write_input_file(data, file_path):
    '''
    Will write an input file to disk.
    data is a dictionary where each item is a separate file that needs to be
    written. Keys are:

        search_queries

    '''
    for key in data.keys():
        if key == SEARCH_QUERIES_KEY:
            krpy.debug.message('Found: \'' + SEARCH_QUERIES_KEY + \
                               '\' in data.', write_input_file)
            search_queries_tree = data[key]
            if os.path.exists(file_path):
                krpy.debug.message('Removing: ' + SEARCH_QUERIES_FILE_NAME + \
                               ' from ' + file_path + '.', write_input_file)
                krpy.io.remove_file_from_zip(file_path,
                                             SEARCH_QUERIES_FILE_NAME)
            temp_dir = tempfile.mkdtemp()
            temp_file_path = os.path.join(temp_dir, 'temp_' + \
                                          SEARCH_QUERIES_FILE_NAME)
            write_search_queries_file(search_queries_tree, temp_file_path)
            zip_file_mode = 'a'
            if not os.path.exists(file_path):
                zip_file_mode = 'w'
            with zipfile.ZipFile(file_path, zip_file_mode) as zip_file:
                krpy.debug.message('Writing: ' + temp_file_path + ' to ' + \
                                   file_path + '.', write_input_file)
                zip_file.write(temp_file_path, SEARCH_QUERIES_FILE_NAME)
            shutil.rmtree(temp_dir)

if __name__ == '__main__':
    if krpy.debug.RUN_DEBUG_CODE:

        INPUT_FILE_DICT = parse_input_file('test_data/sm_input.krsm')

        write_search_queries_file(INPUT_FILE_DICT['search_queries'],
                                  'test_data/search_queries_write_test.tsv')

        write_input_file(INPUT_FILE_DICT, 'test_data/sm_input_test.zip')
