# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import with_statement

import os
import sys

PS = os.path.sep

PYTHON_VERSION_MAJOR = sys.version_info[0]
PYTHON_VERSION_MINOR = sys.version_info[1]
PYTHON_VERSION_MICRO = sys.version_info[2]

PYTHON_VERSION_STRING = '{v0}.{v1}.{v2}'.format(
    v0=PYTHON_VERSION_MAJOR,
    v1=PYTHON_VERSION_MINOR,
    v2=PYTHON_VERSION_MICRO)

# Python 2 and 3 differences
STRING_TYPE = str
FILE_OPEN_MODE_READ = 'r'
if sys.hexversion < 0x03000000:
    STRING_TYPE = basestring
    FILE_OPEN_MODE_READ = 'rU'
