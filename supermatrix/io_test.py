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


def test_parse_input_file():
    from krpy import supermatrix
    out_dict = supermatrix.io.parse_input_file('test_data/sm_input.krsm')
    assert out_dict['search_queries'] is not None
