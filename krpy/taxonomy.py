# -*- coding: utf-8 -*-

"""

This module implements Taxonomy class which uses NCBI databases to
provide various taxonomi name resolution services (TNRS).

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import with_statement

# import os.path
# from os import remove
# import hashlib
# import zipfile

from krpy import PS
# from krpy import FILE_OPEN_MODE_READ

from krpy import Root
from krpy import Error
# from krpy import internet
from krpy import io

# from krpy.ncbi import TAX_BASE_URL
# from krpy.ncbi import TAXDMP_FILES
# from krpy.ncbi import TAXDMP_ARCHIVE
# from krpy.ncbi import TAXCAT_FILES
# from krpy.ncbi import TAXCAT_ARCHIVE
from krpy.ncbi import NAME_CLASS_SET

from krpy.ncbi import _update_ncbi_taxonomy_data
from krpy.ncbi import _parse_codons

from krpy.ncbi import _parse_names_dump
from krpy.ncbi import _parse_delnodes_dump
from krpy.ncbi import _parse_merged_dump
from krpy.ncbi import _parse_nodes_dump
from krpy.ncbi import _parse_gencode_dump

from krpy.ncbi_examine import _examine_ncbi_taxonomy_names_dump
from krpy.ncbi_examine import _examine_ncbi_taxonomy_nodes_dump
from krpy.ncbi_examine import _examine_ncbi_taxonomy_delnodes_dump
from krpy.ncbi_examine import _examine_ncbi_taxonomy_merged_dump
from krpy.ncbi_examine import _examine_ncbi_taxonomy_gencode_dump


class Taxonomy(Root):

    def __init__(self, data_dir_path, check_for_updates=True):

        self._data_dir_path = io.prepare_directory(path=data_dir_path)

        self._tax_dmp_path = self._data_dir_path + 'taxdmp' + PS
        io.prepare_directory(path=self._tax_dmp_path)

        self._tax_cat_path = self._data_dir_path + 'taxcat' + PS
        io.prepare_directory(path=self._tax_cat_path)

        self._tax_names_dmp_path = self._tax_dmp_path + 'names.dmp'
        self._tax_nodes_dmp_path = self._tax_dmp_path + 'nodes.dmp'
        self._tax_delnodes_dmp_path = self._tax_dmp_path + 'delnodes.dmp'
        self._tax_merged_dmp_path = self._tax_dmp_path + 'merged.dmp'
        self._tax_gencode_dmp_path = self._tax_dmp_path + 'gencode.dmp'

        self._tax_gencode_prt_path = self._tax_dmp_path + 'gc.prt'

        self.update(check_for_updates=check_for_updates)

    def update(self, check_for_updates=True):

        print('Updating NCBI taxonomy data if necessary or requested.')

        _update_ncbi_taxonomy_data(
            taxdmp_path=self._tax_dmp_path,
            taxcat_path=self._tax_cat_path,
            force_redownload=False,
            check_for_updates=check_for_updates)

        print('Loading NCBI taxonomy data.')

        self._codons = _parse_codons(
            tax_gencode_prt_path=self._tax_gencode_prt_path)

        self._taxids_names_dict = _parse_names_dump(
            file_path=self._tax_names_dmp_path)

        self._taxids_deleted_set = _parse_delnodes_dump(
            file_path=self._tax_delnodes_dmp_path)

        self._taxids_merged_dict = _parse_merged_dump(
            file_path=self._tax_merged_dmp_path)

        parse_nodes_dump_result = _parse_nodes_dump(
            file_path=self._tax_nodes_dmp_path)

        self._taxids_child_parent_dict = parse_nodes_dump_result[0]
        self._taxids_rank_dict = parse_nodes_dump_result[1]
        self._taxids_genetic_code_id_dict = parse_nodes_dump_result[2]
        self._taxids_mito_genetic_code_id_dict = parse_nodes_dump_result[3]

        parse_gencode_dump_result = _parse_gencode_dump(
            file_path=self._tax_gencode_dmp_path)

        self._gen_code_id_name_dict = parse_gencode_dump_result[0]
        self._gen_code_id_translation_table_dict = parse_gencode_dump_result[1]
        self._gen_code_id_start_codons_dict = parse_gencode_dump_result[2]

        self._name_classes = list(NAME_CLASS_SET)
        self._name_classes.sort()

    # class properties =================================================
    @property
    def codons(self):
        return self._codons

    @property
    def name_classes(self):
        return self._name_classes

    # class methods ====================================================
    def taxid_deleted(self, taxid):
        taxid_deleted = False
        if taxid in self._taxids_deleted_set:
            taxid_deleted = True
        return taxid_deleted

    def taxid_merged(self, taxid):
        taxid_merged = False
        if taxid in self._taxids_merged_dict:
            taxid_merged = self._taxids_merged_dict[taxid]
        return taxid_merged

    def updated_taxid(self, taxid):
        return_taxid = taxid
        taxid_deleted = self.taxid_deleted(taxid=taxid)
        taxid_merged = self.taxid_merged(taxid=taxid)

        if taxid_deleted is False:
            if taxid_merged is not False:
                return_taxid = taxid_merged
        else:
            return_taxid = None

        return return_taxid

    def names_for_taxid(self, taxid):
        return_dict = dict()
        return_dict['old_taxid'] = taxid
        new_taxid = self.updated_taxid(taxid=taxid)
        return_dict['new_taxid'] = new_taxid
        return_dict['names'] = None

        if new_taxid is not None:
            return_dict['names'] = self._taxids_names_dict[new_taxid]

        return return_dict

    def name_class_for_taxid(self, taxid, name_class='scientific name'):

        if name_class not in NAME_CLASS_SET:
            message = 'Name class \'{n}\' is not valid.'
            message = message.format(n=name_class)
            raise Error(message)

        all_names_dict = self.names_for_taxid(taxid=taxid)

        return_dict = dict()
        return_dict['old_taxid'] = all_names_dict['old_taxid']
        return_dict['new_taxid'] = all_names_dict['new_taxid']
        return_dict['name'] = None

        names = all_names_dict['names']

        if names is not None:
            for n in names:
                if n['name_class'] == name_class:
                    if return_dict['name'] is None:
                        return_dict['name'] = list()
                    return_dict['name'].append(n['name'])

        return return_dict

    def parent_taxid(self, taxid):
        return_taxid = None
        taxid = self.updated_taxid(taxid=taxid)
        if taxid in self._taxids_child_parent_dict:
            return_taxid = self._taxids_child_parent_dict[taxid]

        return return_taxid

    def lineage_for_taxid(self, taxid):
        return_dict = dict()
        return_dict['old_taxid'] = taxid
        new_taxid = self.updated_taxid(taxid=taxid)
        return_dict['new_taxid'] = new_taxid
        return_dict['lineage'] = None

        def recurse_lineage(taxid, lineage):
            lineage.append(taxid)
            if taxid != '1':
                taxid = self.parent_taxid(taxid=taxid)
                return recurse_lineage(taxid, lineage)
            else:
                return lineage

        lineage = list()
        if new_taxid is not None:
            lineage = recurse_lineage(
                taxid=new_taxid, lineage=lineage)

        lineage.reverse()
        return_dict['lineage'] = lineage

        return return_dict

    # examine methods ==================================================
    def examine_names_dump(self):
        return _examine_ncbi_taxonomy_names_dump(
            file_path=self._tax_names_dmp_path)

    def examine_nodes_dump(self):
        return _examine_ncbi_taxonomy_nodes_dump(
            file_path=self._tax_nodes_dmp_path)

    def examine_delnodes_dump(self):
        return _examine_ncbi_taxonomy_delnodes_dump(
            file_path=self._tax_delnodes_dmp_path)

    def examine_merged_dump(self):
        return _examine_ncbi_taxonomy_merged_dump(
            file_path=self._tax_merged_dmp_path)

    def examine_gencode_dump(self):
        return _examine_ncbi_taxonomy_gencode_dump(
            file_path=self._tax_gencode_dmp_path)
