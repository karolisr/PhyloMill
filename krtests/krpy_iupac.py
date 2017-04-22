# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import with_statement

import unittest

from krpy import iupac


def setUpModule():
    print('\nsetUpModule krpyIUPACTests')


def tearDownModule():
    print('\n\ntearDownModule krpyIUPACTests')


class krpyIUPACTests(unittest.TestCase):

    # print('krpyIUPACTests')

    def test_iupac(self):

        print('\ntest_iupac')

        self.assertEqual(
            iupac.NT_SHARED_CHARS,
            set(['A', 'C', 'G']))
        self.assertEqual(
            iupac.DNA_ONLY_CHARS,
            set(['T']))
        self.assertEqual(
            iupac.RNA_ONLY_CHARS,
            set(['U']))
        self.assertEqual(
            iupac.NT_AMBIGUOUS_CHARS,
            set(['B', 'D', 'H', 'K', 'M', 'N', 'S', 'R', 'W', 'V',
                 'Y']))
        self.assertEqual(
            iupac.NT_GAPS_CHARS,
            set(['-']))
        self.assertEqual(
            iupac.DNA_UNAMBIGUOUS,
            set(['A', 'C', 'T', 'G']))
        self.assertEqual(
            iupac.DNA_UNAMBIGUOUS_GAPS,
            set(['A', 'C', '-', 'T', 'G']))
        self.assertEqual(
            iupac.DNA_AMBIGUOUS,
            set(['A', 'C', 'B', 'D', 'G', 'H', 'K', 'M', 'N', 'S', 'R',
                 'T', 'W', 'V', 'Y']))
        self.assertEqual(
            iupac.DNA_AMBIGUOUS_GAPS,
            set(['A', 'C', 'B', 'D', 'G', 'H', 'K', '-', 'M', 'N', 'S',
                 'R', 'T', 'W', 'V', 'Y']))
        self.assertEqual(
            iupac.RNA_UNAMBIGUOUS,
            set(['A', 'C', 'U', 'G']))
        self.assertEqual(
            iupac.RNA_UNAMBIGUOUS_GAPS,
            set(['A', '-', 'C', 'U', 'G']))
        self.assertEqual(
            iupac.RNA_AMBIGUOUS,
            set(['A', 'C', 'B', 'D', 'G', 'H', 'K', 'M', 'N', 'S', 'R',
                 'U', 'W', 'V', 'Y']))
        self.assertEqual(
            iupac.RNA_AMBIGUOUS_GAPS,
            set(['A', 'C', 'B', 'D', 'G', 'H', 'K', '-', 'M', 'N', 'S',
                 'R', 'U', 'W', 'V', 'Y']))
        self.assertEqual(
            iupac.NT_UNAMBIGUOUS,
            set(['A', 'C', 'G', 'U', 'T']))
        self.assertEqual(
            iupac.NT_UNAMBIGUOUS_GAPS,
            set(['A', 'C', 'G', '-', 'U', 'T']))
        self.assertEqual(
            iupac.NT_AMBIGUOUS,
            set(['A', 'C', 'B', 'D', 'G', 'H', 'K', 'M', 'N', 'S', 'R',
                 'U', 'T', 'W', 'V', 'Y']))
        self.assertEqual(
            iupac.NT_AMBIGUOUS_GAPS,
            set(['A', 'C', 'B', 'D', 'G', 'H', 'K', 'M', 'N', 'S', 'R',
                 'U', 'T', 'W', 'V', 'Y', '-']))

        self.assertEqual(
            iupac.AA_CHARS,
            set(['G', 'A', 'V', 'L', 'I', 'P', 'F', 'Y', 'C', 'M', 'H',
                 'K', 'R', 'W', 'S', 'T', 'D', 'E', 'N', 'Q', '.']))
        self.assertEqual(
            iupac.AA_AMBIGUOUS_CHARS,
            set(['B', 'Z', 'X']))
        self.assertEqual(
            iupac.AA_GAPS_CHARS,
            set(['-']))
        self.assertEqual(
            iupac.AA_UNAMBIGUOUS,
            set(['G', 'A', 'V', 'L', 'I', 'P', 'F', 'Y', 'C', 'M', 'H',
                 'K', 'R', 'W', 'S', 'T', 'D', 'E', 'N', 'Q', '.']))
        self.assertEqual(
            iupac.AA_UNAMBIGUOUS_GAPS,
            set(['G', 'A', 'V', 'L', 'I', 'P', 'F', 'Y', 'C', 'M', 'H',
                 'K', 'R', 'W', 'S', 'T', 'D', 'E', 'N', 'Q', '.', '-']))
        self.assertEqual(
            iupac.AA_AMBIGUOUS,
            set(['G', 'A', 'V', 'L', 'I', 'P', 'F', 'Y', 'C', 'M', 'H',
                 'K', 'R', 'W', 'S', 'T', 'D', 'E', 'N', 'Q', '.', 'B',
                 'Z', 'X']))
        self.assertEqual(
            iupac.AA_AMBIGUOUS_GAPS,
            set(['G', 'A', 'V', 'L', 'I', 'P', 'F', 'Y', 'C', 'M', 'H',
                 'K', 'R', 'W', 'S', 'T', 'D', 'E', 'N', 'Q', '.', 'B',
                 'Z', 'X', '-']))

        # print('NT_SHARED_CHARS', iupac.NT_SHARED_CHARS)
        # print('DNA_ONLY_CHARS', iupac.DNA_ONLY_CHARS)
        # print('RNA_ONLY_CHARS', iupac.RNA_ONLY_CHARS)
        # print('NT_AMBIGUOUS_CHARS', iupac.NT_AMBIGUOUS_CHARS)
        # print('NT_GAPS_CHARS', iupac.NT_GAPS_CHARS)

        # print('DNA_UNAMBIGUOUS', iupac.DNA_UNAMBIGUOUS)
        # print('DNA_UNAMBIGUOUS_GAPS', iupac.DNA_UNAMBIGUOUS_GAPS)
        # print('DNA_AMBIGUOUS', iupac.DNA_AMBIGUOUS)
        # print('DNA_AMBIGUOUS_GAPS', iupac.DNA_AMBIGUOUS_GAPS)

        # print('RNA_UNAMBIGUOUS', iupac.RNA_UNAMBIGUOUS)
        # print('RNA_UNAMBIGUOUS_GAPS', iupac.RNA_UNAMBIGUOUS_GAPS)
        # print('RNA_AMBIGUOUS', iupac.RNA_AMBIGUOUS)
        # print('RNA_AMBIGUOUS_GAPS', iupac.RNA_AMBIGUOUS_GAPS)

        # print('NT_UNAMBIGUOUS', iupac.NT_UNAMBIGUOUS)
        # print('NT_UNAMBIGUOUS_GAPS', iupac.NT_UNAMBIGUOUS_GAPS)
        # print('NT_AMBIGUOUS', iupac.NT_AMBIGUOUS)
        # print('NT_AMBIGUOUS_GAPS', iupac.NT_AMBIGUOUS_GAPS)

        # print('AA_CHARS', iupac.AA_CHARS)
        # print('AA_AMBIGUOUS_CHARS', iupac.AA_AMBIGUOUS_CHARS)
        # print('AA_GAPS_CHARS', iupac.AA_GAPS_CHARS)

        # print('AA_UNAMBIGUOUS', iupac.DNA_UNAMBIGUOUS)
        # print('AA_UNAMBIGUOUS_GAPS', iupac.DNA_UNAMBIGUOUS_GAPS)
        # print('AA_AMBIGUOUS', iupac.DNA_AMBIGUOUS)
        # print('AA_AMBIGUOUS_GAPS', iupac.DNA_AMBIGUOUS_GAPS)
