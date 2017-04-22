#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import with_statement

import unittest

# import datetime

from krpy import PYTHON_VERSION_STRING
py_ver_msg = '\nPython version: {pv}\n'.format(pv=PYTHON_VERSION_STRING)
print(py_ver_msg)

from krtests import *
# from krtests import krpySeqTests, krpyIUPACTests

def main():

    unittest.main()

    # start_time = datetime.datetime.now()
    # end_time = datetime.datetime.now()
    # time_taken = end_time - start_time
    # print('Time taken by tests:', time_taken)

if __name__ == '__main__':
    main()
