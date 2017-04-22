# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import with_statement

import unittest

from krpy import PS
from krpy import FILE_OPEN_MODE_READ
from krpy.bioio import _parse_gbseq_xml_handle

from krtests import TEST_DATA_DIR_PATH


def setUpModule():
    print('\nsetUpModule krpyBioioTests')


def tearDownModule():
    print('\n\ntearDownModule krpyBioioTests')


class krpyBioioTests(unittest.TestCase):

    # print('krpyBioioTests')

    def test_parse_gbseq_xml_handle(self):

        print('\ntest_parse_gbseq_xml_handle')

        with open(TEST_DATA_DIR_PATH + 'gbseq_xml_sample.xml',
                  FILE_OPEN_MODE_READ) as f:
            results = _parse_gbseq_xml_handle(f)

        self.assertEqual(len(results), 6)

        for r in results:

            self.assertEqual(len(r.keys()), 15)

            self.assertTrue('accession' in r)
            self.assertTrue('date_create' in r)
            self.assertTrue('date_update' in r)
            self.assertTrue('definition' in r)
            self.assertTrue('division' in r)
            self.assertTrue('features' in r)
            # self.assertTrue('gi' in r)
            self.assertTrue('length' in r)
            self.assertTrue('mol_type' in r)
            self.assertTrue('organism' in r)
            self.assertTrue('seq' in r)
            self.assertTrue('strandedness' in r)
            self.assertTrue('taxid' in r)
            self.assertTrue('taxonomy' in r)
            self.assertTrue('topology' in r)
            self.assertTrue('version' in r)

            # # Print raw data
            # for k in r.keys():
            #     if k == 'features':
            #         pass
            #         features = r[k]
            #         print('\nfeatures:')
            #         for fk in features.keys():
            #             print(fk, '::', features[fk])
            #         print()
            #     elif k == 'seq':
            #         pass
            #     elif k == 'taxonomy':
            #         pass
            #     else:
            #         print(k, '::', r[k])
            # print('\n')
