# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import with_statement

import unittest
from datetime import datetime

from krpy.seq import Seq
from krpy.seq import SeqRecord

from krpy.seq import NTSeq
from krpy.seq import DNASeq
from krpy.seq import RNASeq
from krpy.seq import AASeq

from krpy.seq import SEQ_TYPE_NT
from krpy.seq import SEQ_TYPE_DNA
from krpy.seq import SEQ_TYPE_RNA
from krpy.seq import SEQ_TYPE_AA

from krpy.seq import SEQ_TYPES
from krpy.seq import MOL_TO_SEQ_TYPE_MAP

from krpy import PS
from krpy import FILE_OPEN_MODE_READ
from krpy.bioio import _parse_gbseq_xml_handle

from krtests import TEST_DATA_DIR_PATH


def setUpModule():
    print('\nsetUpModule krpySeqTests')


def tearDownModule():
    print('\n\ntearDownModule krpySeqTests')


class krpySeqTests(unittest.TestCase):

    # print('krpySeqTests')

    def test_Seq(self):

        print('\ntest_Seq')

        seq = Seq(seq='ACGTUBDHKMNRSVWY', seq_type=SEQ_TYPE_NT)
        self.assertTrue(isinstance(seq, NTSeq))
        self.assertEqual(seq.length, 16)

        seq = Seq(seq='ACGTBDHKMNRSVWY', seq_type=SEQ_TYPE_DNA)
        self.assertTrue(isinstance(seq, DNASeq))
        self.assertEqual(seq.length, 15)

        seq = Seq(seq='ACGUBDHKMNRSVWY', seq_type=SEQ_TYPE_RNA)
        self.assertTrue(isinstance(seq, RNASeq))
        self.assertEqual(seq.length, 15)

        seq = Seq(seq='GAVLIPFYCMHKRWSTDENQ.BZX', seq_type=SEQ_TYPE_AA)
        self.assertTrue(isinstance(seq, AASeq))
        self.assertEqual(seq.length, 24)

    def test_SeqRecord(self):

        print('\ntest_SeqRecord')

        seq_record = SeqRecord(seq='ACGTUBDHKMNRSVWY', mol_type=SEQ_TYPE_NT)
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, NTSeq))
        self.assertEqual(seq_record.seq.length, 16)

        seq_record = SeqRecord(seq='ACGTBDHKMNRSVWY', mol_type=SEQ_TYPE_DNA)
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, DNASeq))
        self.assertEqual(seq_record.seq.length, 15)

        seq_record = SeqRecord(seq='ACGUBDHKMNRSVWY', mol_type=SEQ_TYPE_RNA)
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, RNASeq))
        self.assertEqual(seq_record.seq.length, 15)

        seq_record = SeqRecord(seq='ACGTBDHKMNRSVWYU', mol_type=SEQ_TYPE_DNA)
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, DNASeq))
        self.assertEqual(seq_record.seq.length, 16)

        seq_record = SeqRecord(seq='ACGUBDHKMNRSVWYT', mol_type=SEQ_TYPE_RNA)
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, RNASeq))
        self.assertEqual(seq_record.seq.length, 16)

        seq_record = SeqRecord(seq='GAVLIPFYCMHKRWSTDENQ.BZX', mol_type=SEQ_TYPE_AA)
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, AASeq))
        self.assertEqual(seq_record.seq.length, 24)

        seq_record = SeqRecord(seq='ACGTBDHKMNRSVWY', mol_type='DNA')
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, DNASeq))
        self.assertEqual(seq_record.seq.length, 15)

        seq_record = SeqRecord(seq='ACGUBDHKMNRSVWY', mol_type='mRNA')
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, RNASeq))
        self.assertEqual(seq_record.seq.length, 15)

        seq_record = SeqRecord(seq='ACGUBDHKMNRSVWY', mol_type='rRNA')
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, RNASeq))
        self.assertEqual(seq_record.seq.length, 15)

        seq_record = SeqRecord(seq='ACGUBDHKMNRSVWY', mol_type='tRNA')
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, RNASeq))
        self.assertEqual(seq_record.seq.length, 15)

        seq_record = SeqRecord(seq='GAVLIPFYCMHKRWSTDENQ.BZX', mol_type='AA')
        self.assertTrue(isinstance(seq_record, SeqRecord))
        self.assertTrue(isinstance(seq_record.seq, AASeq))
        self.assertEqual(seq_record.seq.length, 24)

        with open(TEST_DATA_DIR_PATH + 'gbseq_xml_sample.xml',
                  FILE_OPEN_MODE_READ) as f:
            results = _parse_gbseq_xml_handle(f)

        self.assertEqual(len(results), 6)

        for r in results:

            self.assertEqual(len(r.keys()), 15)

            # # Print raw data
            # for k in r.keys():
            #     if k == 'features':
            #         pass
            #         features = r[k]
            #         print('\nfeatures:')
            #         for fk in features.keys():
            #             print(fk, '::', features[fk])
            #         print()
            #     # elif k == 'seq':
            #     #     pass
            #     # elif k == 'taxonomy':
            #     #     pass
            #     else:
            #         print(k, '::', r[k])
            # print('\n')

            seq_record = SeqRecord(
                seq=r['seq'],
                mol_type=r['mol_type'],
                accession=r['accession'],
                version=r['version'],
                description=r['definition'],
                strandedness=r['strandedness'],
                topology=r['topology'],
                division=r['division'],
                date_create=(datetime.strptime(r['date_create'], '%d-%b-%Y')).strftime('%Y-%m-%d'),
                date_update=(datetime.strptime(r['date_update'], '%d-%b-%Y')).strftime('%Y-%m-%d'),
                taxid=r['taxid'],
                organism=r['organism'],
                taxonomy=r['taxonomy'],
                features=r['features']
                )

            # self.assertTrue('accession' in r)
            # self.assertTrue('date_create' in r)
            # self.assertTrue('date_update' in r)
            # self.assertTrue('definition' in r)
            # self.assertTrue('division' in r)
            # self.assertTrue('features' in r)
            # self.assertTrue('length' in r)
            # self.assertTrue('mol_type' in r)
            # self.assertTrue('organism' in r)
            # self.assertTrue('seq' in r)
            # self.assertTrue('strandedness' in r)
            # self.assertTrue('taxid' in r)
            # self.assertTrue('taxonomy' in r)
            # self.assertTrue('topology' in r)
            # self.assertTrue('version' in r)


            self.assertTrue(isinstance(seq_record, SeqRecord))

            seq_type = None

            if r['mol_type'] in SEQ_TYPES:
                seq_type = r['mol_type']
            elif r['mol_type'] in MOL_TO_SEQ_TYPE_MAP:
                seq_type = MOL_TO_SEQ_TYPE_MAP[r['mol_type']]

            if seq_type is SEQ_TYPE_NT:
                self.assertTrue(isinstance(seq_record.seq, NTSeq))
            elif seq_type is SEQ_TYPE_DNA:
                self.assertTrue(isinstance(seq_record.seq, DNASeq))
            elif seq_type is SEQ_TYPE_RNA:
                self.assertTrue(isinstance(seq_record.seq, RNASeq))
            elif seq_type is SEQ_TYPE_AA:
                self.assertTrue(isinstance(seq_record.seq, AASeq))

            # for attr in dir(seq_record):
            #     if not attr.startswith('_'):
            #         print(attr, getattr(seq_record, attr))

            # print('\n----------------------------------------------------------------------\n')