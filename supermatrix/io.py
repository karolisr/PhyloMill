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
__date__ = '2013-12-26'
__updated__ = '2013-12-26'


def parse_input_file(file_path):

    '''
        Read an input (zip) file and return correct representations of all
        input files.
    '''

    import zipfile
    import krpy

    search_queries_file_name = 'search_queries.tsv'

    with zipfile.ZipFile(file_path, 'r') as zip_file:
        krpy.debug.message('Reading file: ' + file_path + '.',
                           parse_input_file)
        for file_name in zip_file.namelist():
            krpy.debug.message('Reading file: ' + file_name + ' in ' +
                               file_path + '.', parse_input_file)
            if file_name == search_queries_file_name:
                pass

    return()

if __name__ == '__main__':

    parse_input_file('test_data/sm_input.krsm')
