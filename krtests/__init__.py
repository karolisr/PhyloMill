# -*- coding: utf-8 -*-

"""

This package contains unit tests for packages krpy and krpm.

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import with_statement

from .config import DATA_DIR_PATH
from .config import TEST_DATA_DIR_PATH

from .krpy_entrez import krpyEntrezTests
from .krpy_seq import krpySeqTests
from .krpy_bioio import krpyBioioTests
from .krpy_iupac import krpyIUPACTests
from .krpy_ncbi import krpyNCBITests
