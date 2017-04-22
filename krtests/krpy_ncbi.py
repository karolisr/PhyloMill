# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import with_statement

import unittest

from krpy import Error
from krpy import Taxonomy

from krtests import DATA_DIR_PATH


def setUpModule():
    print('\nsetUpModule krpyNCBITests')


def tearDownModule():
    print('\n\ntearDownModule krpyNCBITests')


class krpyNCBITests(unittest.TestCase):

    # print('krpyNCBITests')

    taxonomy = Taxonomy(
        data_dir_path=DATA_DIR_PATH + 'taxonomy',
        check_for_updates=False)

    def test_codons(self):

        run = True
        if run:

            print('\ntest_codons\n')

            return_value = self.taxonomy.codons

            self.assertEqual(len(return_value), 64)
            self.assertEqual(''.join(return_value[0]), 'TTT')
            self.assertEqual(''.join(return_value[63]), 'GGG')

    def test_name_classes(self):

        run = True
        if run:

            print('\ntest_name_classes\n')

            return_value = self.taxonomy.name_classes
            self.assertEqual(
                return_value,
                ['acronym', 'anamorph', 'authority', 'blast name',
                 'common name', 'equivalent name', 'genbank acronym',
                 'genbank anamorph', 'genbank common name',
                 'genbank synonym', 'in-part', 'includes', 'misnomer',
                 'misspelling', 'scientific name', 'synonym',
                 'teleomorph', 'type material'])

    def test_taxid_deleted(self):

        run = True
        if run:

            print('\ntest_taxid_deleted\n')

            # deleted
            return_value = self.taxonomy.taxid_deleted('736973')
            self.assertTrue(return_value)

            # merged
            return_value = self.taxonomy.taxid_deleted('360')
            self.assertFalse(return_value)

            # not deleted or merged
            return_value = self.taxonomy.taxid_deleted('707246')
            self.assertFalse(return_value)

    def test_taxid_merged(self):

        run = True
        if run:

            print('\ntest_taxid_merged\n')

            # deleted
            return_value = self.taxonomy.taxid_merged('736973')
            self.assertFalse(return_value)

            # merged
            return_value = self.taxonomy.taxid_merged('360')
            self.assertTrue(return_value)

            # not deleted or merged
            return_value = self.taxonomy.taxid_merged('707246')
            self.assertFalse(return_value)

    def test_updated_taxid(self):

        run = True
        if run:

            print('\ntest_updated_taxid\n')

            # deleted
            return_value = self.taxonomy.updated_taxid('736973')
            self.assertEqual(return_value, None)

            # merged
            return_value = self.taxonomy.updated_taxid('360')
            self.assertEqual(return_value, '359')

            # not deleted or merged
            return_value = self.taxonomy.updated_taxid('707246')
            self.assertEqual(return_value, '707246')

    def test_names_for_taxid(self):

        run = True
        if run:

            print('\ntest_names_for_taxid\n')

            # deleted
            return_dict = self.taxonomy.names_for_taxid('736973')
            self.assertEqual(return_dict['old_taxid'], '736973')
            self.assertEqual(return_dict['new_taxid'], None)
            self.assertEqual(return_dict['names'], None)

            # merged
            return_dict = self.taxonomy.names_for_taxid('360')
            self.assertEqual(return_dict['old_taxid'], '360')
            self.assertEqual(return_dict['new_taxid'], '359')
            expected_name_list = [

                {'name_class': 'type material',
                 'name': 'ATCC 11325',
                 'unique_name': ''},

                {'name_class': 'synonym',
                 'name': 'Agrobacterium biovar 2',
                 'unique_name': ''},

                {'name_class': 'synonym',
                 'name': 'Agrobacterium genomic group 10',
                 'unique_name': ''},

                {'name_class': 'synonym',
                 'name': 'Agrobacterium genomic species 10',
                 'unique_name': ''},

                {'name_class': 'synonym',
                 'name': 'Agrobacterium genomosp. 10',
                 'unique_name': ''},

                {'name_class': 'scientific name',
                 'name': 'Agrobacterium rhizogenes',
                 'unique_name': ''},

                {'name_class': 'includes',
                 'name': 'Agrobacterium rhizogenes (RI plasmid PRI1724)',
                 'unique_name': ''},

                {'name_class': 'includes',
                 'name': 'Agrobacterium rhizogenes (RI plasmid PRI8196)',
                 'unique_name': ''},

                {'name_class': 'includes',
                 'name': 'Agrobacterium rhizogenes (RI plasmid PRIA4B)',
                 'unique_name': ''},

                {'name_class': 'authority',
                 'name': ('Agrobacterium rhizogenes (Riker et al. 1930) Conn '
                          '1942 (Approved Lists 1980) emend. '
                          'Sawada et al. 1993'),
                 'unique_name': ''},

                {'name_class': 'type material',
                 'name': 'CFBP 5520',
                 'unique_name': ''},

                {'name_class': 'type material',
                 'name': 'CIP 104328',
                 'unique_name': ''},

                {'name_class': 'type material',
                 'name': 'DSM 30148',
                 'unique_name': ''},

                {'name_class': 'type material',
                 'name': 'HAMBI 1816',
                 'unique_name': ''},

                {'name_class': 'type material',
                 'name': 'ICMP 5794',
                 'unique_name': ''},

                {'name_class': 'type material',
                 'name': 'IFO 13257',
                 'unique_name': ''},

                {'name_class': 'type material',
                 'name': 'JCM 20919',
                 'unique_name': ''},

                {'name_class': 'type material',
                 'name': 'LMG 150',
                 'unique_name': ''},

                {'name_class': 'type material',
                 'name': 'NBRC 13257',
                 'unique_name': ''},

                {'name_class': 'type material',
                 'name': 'NCPPB 2991',
                 'unique_name': ''},

                {'name_class': 'genbank synonym',
                 'name': 'Rhizobium rhizogenes',
                 'unique_name': ''},

                {'name_class': 'authority',
                 'name': ('Rhizobium rhizogenes (Riker et al. 1930) '
                          'Young et al. 2001'),
                 'unique_name': ''},

                {'name_class': 'includes',
                 'name': 'Rhizobium sp. LMG 9509',
                 'unique_name': ''}]

            self.assertEqual(return_dict['names'], expected_name_list)

            # not deleted or merged
            return_dict = self.taxonomy.names_for_taxid('707246')
            self.assertEqual(return_dict['old_taxid'], '707246')
            self.assertEqual(return_dict['new_taxid'], '707246')

            expected_name_list = [

                {'name_class': 'scientific name',
                 'name': 'Inocybe rhodiola',
                 'unique_name': ''},

                {'name_class': 'authority',
                 'name': 'Inocybe rhodiola Bres. 1887',
                 'unique_name': ''}]

            self.assertEqual(return_dict['names'], expected_name_list)

    def test_name_class_for_taxid(self):

        run = True
        if run:

            print('\ntest_name_class_for_taxid\n')

            self.assertRaises(
                Error,
                self.taxonomy.name_class_for_taxid,
                taxid='1',
                name_class='strange')

            return_dict = self.taxonomy.name_class_for_taxid(
                taxid='736973', name_class='scientific name')
            self.assertEqual(
                return_dict,
                {'name': None, 'old_taxid': '736973', 'new_taxid': None})

            return_dict = self.taxonomy.name_class_for_taxid(
                taxid='707246', name_class='scientific name')
            self.assertEqual(
                return_dict,
                {'name': ['Inocybe rhodiola'], 'old_taxid': '707246',
                 'new_taxid': '707246'})

            return_dict = self.taxonomy.name_class_for_taxid(
                taxid='4070', name_class='scientific name')
            self.assertEqual(
                return_dict,
                {'name': ['Solanaceae'], 'old_taxid': '4070',
                 'new_taxid': '4070'})

            return_dict = self.taxonomy.name_class_for_taxid(
                taxid='4070', name_class='common name')
            self.assertEqual(
                return_dict,
                {'name': ['nightshade family'], 'old_taxid': '4070',
                 'new_taxid': '4070'})

            return_dict = self.taxonomy.name_class_for_taxid(
                taxid='4081', name_class='scientific name')
            self.assertEqual(
                return_dict,
                {'name': ['Solanum lycopersicum'], 'old_taxid': '4081',
                 'new_taxid': '4081'})

            return_dict = self.taxonomy.name_class_for_taxid(
                taxid='4081', name_class='common name')
            self.assertEqual(
                return_dict,
                {'name': None, 'old_taxid': '4081', 'new_taxid': '4081'})

            return_dict = self.taxonomy.name_class_for_taxid(
                taxid='1', name_class='scientific name')
            self.assertEqual(
                return_dict,
                {'name': ['root'], 'old_taxid': '1', 'new_taxid': '1'})

            return_dict = self.taxonomy.name_class_for_taxid(
                taxid='1', name_class='common name')
            self.assertEqual(
                return_dict,
                {'name': None, 'old_taxid': '1', 'new_taxid': '1'})

    def test_parent_taxid(self):

        run = True
        if run:

            print('\ntest_parent_taxid\n')

            return_value = self.taxonomy.parent_taxid('1')
            self.assertEqual(return_value, '1')

            return_value = self.taxonomy.parent_taxid('2')
            self.assertEqual(return_value, '131567')

            return_value = self.taxonomy.parent_taxid('1343429')
            self.assertEqual(return_value, '212542')

    def test_lineage_for_taxid(self):

        run = True
        if run:

            print('\ntest_lineage_for_taxid\n')

            return_value = self.taxonomy.lineage_for_taxid('1')
            self.assertEqual(
                return_value,
                {'new_taxid': '1', 'lineage': ['1'], 'old_taxid': '1'})

            return_value = self.taxonomy.lineage_for_taxid('2')
            self.assertEqual(
                return_value,
                {'new_taxid': '2',
                 'old_taxid': '2',
                 'lineage': ['1', '131567', '2']})

            return_value = self.taxonomy.lineage_for_taxid('1343429')
            self.assertEqual(
                return_value,
                {'new_taxid': '1343429',
                 'old_taxid': '1343429',
                 'lineage': ['1', '131567', '2759', '33154', '33208',
                             '6072', '33213', '33317', '1206794',
                             '88770', '6656', '197563', '197562',
                             '6960', '50557', '85512', '7496', '33340',
                             '33392', '7041', '41071', '535382',
                             '41073', '27450', '1324880', '27452',
                             '212542', '1343429']})

            return_value = self.taxonomy.lineage_for_taxid('4070')
            self.assertEqual(
                return_value,
                {'new_taxid': '4070',
                 'old_taxid': '4070',
                 'lineage': ['1', '131567', '2759', '33090', '35493',
                             '131221', '3193', '58023', '78536',
                             '58024', '3398', '1437183', '71240',
                             '91827', '1437201', '71274', '91888',
                             '4069', '4070']})

    def test_examine_ncbi_taxonomy_names_dump(self):

        run = True
        if run:

            print('\ntest_examine_ncbi_taxonomy_names_dump\n')

            return_dict = self.taxonomy.examine_names_dump()

            column_count_set_expected = set([4])
            column_count_set_returned = return_dict['column_count_set']
            self.assertEqual(column_count_set_returned,
                             column_count_set_expected)

            name_class_set_expected = set(
                ['acronym', 'teleomorph', 'scientific name', 'synonym',
                 'blast name', 'misspelling', 'in-part', 'includes',
                 'genbank acronym', 'misnomer', 'common name',
                 'equivalent name', 'type material', 'authority',
                 'genbank common name', 'genbank synonym', 'anamorph',
                 'genbank anamorph'])
            name_class_set_returned = return_dict['name_class_set']
            self.assertEqual(name_class_set_returned, name_class_set_expected)

            tiwmsns = return_dict['tax_ids_with_multiple_scientific_names_set']
            self.assertEqual(tiwmsns, set())

            print()

    def test_examine_ncbi_taxonomy_nodes_dump(self):

        run = True
        if run:

            print('\ntest_examine_ncbi_taxonomy_nodes_dump\n')

            return_dict = self.taxonomy.examine_nodes_dump()

            column_count_set_expected = set([12, 13])
            column_count_set_returned = return_dict['column_count_set']
            self.assertEqual(column_count_set_returned,
                             column_count_set_expected)

            rank_set_expected = set(['subfamily', 'species group',
                                     'class', 'subgenus', 'subtribe',
                                     'phylum', 'forma', 'superkingdom',
                                     'superorder', 'kingdom', 'superphylum',
                                     'subclass', 'subkingdom', 'no rank',
                                     'superclass', 'subspecies',
                                     'superfamily', 'parvorder', 'infraorder',
                                     'order', 'family', 'subphylum',
                                     'infraclass', 'species subgroup',
                                     'genus', 'species', 'suborder', 'tribe',
                                     'varietas'])
            rank_set_returned = return_dict['rank_set']
            self.assertEqual(rank_set_returned, rank_set_expected)

            genetic_code_id_set_expected = set(['3', '11', '12', '10', '1',
                                                '4', '25', '6', '5', '26'])
            genetic_code_id_set_returned = return_dict['genetic_code_id_set']
            self.assertEqual(genetic_code_id_set_returned,
                             genetic_code_id_set_expected)

            mitochondrial_genetic_code_id_set_expected = set(
                ['3', '5', '21', '11', '16', '0', '13', '14', '1', '9',
                 '23', '4', '22', '2', '24'])
            mitochondrial_genetic_code_id_set_returned = \
                return_dict['mitochondrial_genetic_code_id_set']
            self.assertEqual(mitochondrial_genetic_code_id_set_returned,
                             mitochondrial_genetic_code_id_set_expected)

            print()

    def test_examine_ncbi_taxonomy_delnodes_dump(self):

        run = True
        if run:

            print('\ntest_examine_ncbi_taxonomy_delnodes_dump\n')

            return_dict = self.taxonomy.examine_delnodes_dump()

            column_count_set_expected = set([1])
            column_count_set_returned = return_dict['column_count_set']
            self.assertEqual(column_count_set_returned,
                             column_count_set_expected)

            tax_id_set_returned = return_dict['tax_id_set']
            row_count_returned = return_dict['row_count']
            self.assertEqual(len(tax_id_set_returned),
                             row_count_returned)

    def test_examine_ncbi_taxonomy_merged_dump(self):

        run = True
        if run:

            print('\ntest_examine_ncbi_taxonomy_merged_dump\n')

            return_dict = self.taxonomy.examine_merged_dump()

            column_count_set_expected = set([2])
            column_count_set_returned = return_dict['column_count_set']
            self.assertEqual(column_count_set_returned,
                             column_count_set_expected)

            old_tax_id_set_returned = return_dict['old_tax_id_set']
            # new_tax_id_set_returned = return_dict['new_tax_id_set']
            row_count_returned = return_dict['row_count']
            self.assertEqual(len(old_tax_id_set_returned),
                             row_count_returned)

    def test_examine_ncbi_taxonomy_gencode_dump(self):

        run = True
        if run:

            print('\ntest_examine_ncbi_taxonomy_gencode_dump\n')

            return_dict = self.taxonomy.examine_gencode_dump()

            column_count_set_expected = set([5])
            column_count_set_returned = return_dict['column_count_set']
            self.assertEqual(column_count_set_returned,
                             column_count_set_expected)

            genetic_code_id_set_returned = return_dict['genetic_code_id_set']
            name_set_returned = return_dict['name_set']
            row_count_returned = return_dict['row_count']
            self.assertEqual(len(genetic_code_id_set_returned),
                             row_count_returned)
            self.assertEqual(len(name_set_returned),
                             row_count_returned)

            codon_count_set_expected = set([0, 64])
            codon_count_set_returned = return_dict['codon_count_set']
            self.assertEqual(codon_count_set_returned,
                             codon_count_set_expected)
