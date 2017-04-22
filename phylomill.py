#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import with_statement

# import sys
# import os

from krpy import PYTHON_VERSION_STRING

__version__ = '0.0.0'


def main():

    py_ver_msg = '\nPython version: {pv}\n'.format(pv=PYTHON_VERSION_STRING)
    print(py_ver_msg)

if __name__ == '__main__':
    main()
