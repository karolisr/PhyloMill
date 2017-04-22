# -*- coding: utf-8 -*-

"""

This module defines exception classes for krpy package.

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import with_statement

from krpy import Root


class Error(Exception, Root):

    """

    This is an abstract exception class.

    """

    def __init__(self, message):
        self._message = message

    def __str__(self):
        return self._message
