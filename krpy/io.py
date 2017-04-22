# -*- coding: utf-8 -*-

"""

This module fascilitates basic input/output operations.

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import with_statement

import os
import hashlib

from krpy import PS
from krpy import Error


def clean_path(path):
    return path.rstrip(PS) + PS


def prepare_directory(path):

    if not os.path.exists(path):
        os.makedirs(path)
    elif os.path.isfile(path):
        message = ('Provided path is a file: {f}')
        message = message.format(f=path)
        raise Error(message)

    path = clean_path(path)

    return path


def generate_md5_hash_for_file(file_path):
    return_value = None
    with open(file_path, 'rb') as f:
        return_value = hashlib.md5(f.read()).hexdigest()
    return return_value
