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

from krpy import io
from krpy import PS

user_home_dir = os.path.expanduser('~')
user_home_dir = io.clean_path(path=user_home_dir)
DATA_DIR_PATH = user_home_dir + 'phylomill_test_data' + PS

TEST_DATA_DIR_PATH = os.path.dirname(
    os.path.abspath(sys.argv[0])) + PS + 'krtests' + PS + 'data' + PS
