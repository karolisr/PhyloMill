# -*- coding: utf-8 -*-

"""

This module is for internal use only. It fascilitates NCBI database
consistency testing and examination of table dump files. This code could
reside in ncbi module itself, but is moved here instea to make ncbi
module cleaner.

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

# from krpy import PS
# from krpy import FILE_OPEN_MODE_READ

# from krpy import Root
# from krpy import Error
# from krpy import internet
# from krpy import io

# from krpy.ncbi import TAX_BASE_URL
# from krpy.ncbi import TAXDMP_FILES
# from krpy.ncbi import TAXDMP_ARCHIVE
# from krpy.ncbi import TAXCAT_FILES
# from krpy.ncbi import TAXCAT_ARCHIVE

from krpy.ncbi import _parse_ncbi_taxonomy_dump_file


def _examine_ncbi_taxonomy_names_dump(file_path):
    rows = _parse_ncbi_taxonomy_dump_file(file_path=file_path)

    row_count = len(rows)

    column_count_set = set()

    tax_id_set = set()
    name_set = set()
    unique_name_set = set()
    name_class_set = set()

    tax_id_keyed_dict = dict()

    for r in rows:
        column_count = len(r)
        column_count_set.add(column_count)

        tax_id = r[0]
        name = r[1]
        unique_name = r[2]
        name_class = r[3]

        # if unique_name:
        #     print(name, '::', unique_name, '::', name_class)

        tax_id_set.add(tax_id)
        name_set.add(name)
        unique_name_set.add(unique_name)
        name_class_set.add(name_class)

        if tax_id not in tax_id_keyed_dict:
            tax_id_keyed_dict[tax_id] = list()

        dict_entry = {
            'name': name,
            'unique_name': unique_name,
            'name_class': name_class}

        tax_id_keyed_dict[tax_id].append(dict_entry)

    tax_ids_with_multiple_scientific_names_set = set()
    entries_per_tax_id_set = set()
    for k in tax_id_keyed_dict.keys():
        entries = tax_id_keyed_dict[k]
        entry_count = len(entries)
        entries_per_tax_id_set.add(entry_count)

        found_scientific_name = False
        for entry in entries:
            name = entry['name']
            unique_name = entry['unique_name']
            name_class = entry['name_class']

            name_class_is_scientific_name = (name_class == 'scientific name')

            if found_scientific_name and name_class_is_scientific_name:
                tax_ids_with_multiple_scientific_names_set.add(k)
                break

            if name_class_is_scientific_name:
                found_scientific_name = True

    print('row_count', row_count)
    print('column_count_set', column_count_set)

    print('len(tax_id_set)', len(tax_id_set))
    print('len(name_set)', len(name_set))
    print('len(unique_name_set)', len(unique_name_set))
    print('name_class_set', name_class_set)

    print('len(tax_id_keyed_dict)', len(tax_id_keyed_dict))
    print('entries_per_tax_id_set', entries_per_tax_id_set)
    print('tax_ids_with_multiple_scientific_names_set',
          tax_ids_with_multiple_scientific_names_set)

    return_dict = dict()

    return_dict['column_count_set'] = column_count_set
    return_dict['name_class_set'] = name_class_set
    return_dict['tax_ids_with_multiple_scientific_names_set'] = \
        tax_ids_with_multiple_scientific_names_set

    return return_dict


def _examine_ncbi_taxonomy_nodes_dump(file_path):
    rows = _parse_ncbi_taxonomy_dump_file(file_path=file_path)

    row_count = len(rows)

    column_count_set = set()

    tax_id_set = set()
    parent_tax_id_set = set()
    rank_set = set()
    embl_code_set = set()
    division_id_set = set()
    genetic_code_id_set = set()
    mitochondrial_genetic_code_id_set = set()

    genetic_code_id_list = list()
    mitochondrial_genetic_code_id_list = list()

    for r in rows:
        column_count = len(r)
        column_count_set.add(column_count)

        tax_id = r[0]
        parent_tax_id = r[1]
        rank = r[2]
        embl_code = r[3]
        division_id = r[4]
        # inherited_div_flag = r[5]
        genetic_code_id = r[6]
        # inherited_GC_flag = r[7]
        mitochondrial_genetic_code_id = r[8]
        # inherited_MGC_flag = r[9]
        # GenBank_hidden_flag = r[10]
        # hidden_subtree_root_flag = r[11]

        # not every row contains this column, and it is not needed
        # comments = r[12]

        tax_id_set.add(tax_id)
        parent_tax_id_set.add(parent_tax_id)
        rank_set.add(rank)
        embl_code_set.add(embl_code)
        division_id_set.add(division_id)
        genetic_code_id_set.add(genetic_code_id)
        mitochondrial_genetic_code_id_set.add(mitochondrial_genetic_code_id)

        genetic_code_id_list.append(genetic_code_id)
        mitochondrial_genetic_code_id_list.append(
            mitochondrial_genetic_code_id)

    print('row_count', row_count)
    print('column_count_set', column_count_set)

    print('len(tax_id_set)', len(tax_id_set))
    print('len(parent_tax_id_set)', len(parent_tax_id_set))

    print('len(rank_set)', len(rank_set))
    print('rank_set', rank_set)

    print('len(embl_code_set)', len(embl_code_set))

    print('len(division_id_set)', len(division_id_set))
    print('division_id_set', division_id_set)

    print('len(genetic_code_id_list)', len(genetic_code_id_list))
    print('len(genetic_code_id_set)', len(genetic_code_id_set))
    print('genetic_code_id_set', genetic_code_id_set)

    print('len(mitochondrial_genetic_code_id_list)',
          len(mitochondrial_genetic_code_id_list))
    print('len(mitochondrial_genetic_code_id_set)',
          len(mitochondrial_genetic_code_id_set))
    print('mitochondrial_genetic_code_id_set',
          mitochondrial_genetic_code_id_set)

    return_dict = dict()

    return_dict['column_count_set'] = column_count_set
    return_dict['rank_set'] = rank_set
    return_dict['genetic_code_id_set'] = genetic_code_id_set
    return_dict['mitochondrial_genetic_code_id_set'] = \
        mitochondrial_genetic_code_id_set

    return return_dict


def _examine_ncbi_taxonomy_delnodes_dump(file_path):
    rows = _parse_ncbi_taxonomy_dump_file(file_path=file_path)

    row_count = len(rows)

    column_count_set = set()

    tax_id_set = set()

    for r in rows:
        column_count = len(r)
        column_count_set.add(column_count)

        tax_id = r[0]
        tax_id_set.add(tax_id)

    print('row_count', row_count)
    print('column_count_set', column_count_set)

    print('len(tax_id_set)', len(tax_id_set))

    return_dict = dict()

    return_dict['column_count_set'] = column_count_set
    return_dict['tax_id_set'] = tax_id_set
    return_dict['row_count'] = row_count

    return return_dict


def _examine_ncbi_taxonomy_merged_dump(file_path):
    rows = _parse_ncbi_taxonomy_dump_file(file_path=file_path)

    row_count = len(rows)

    column_count_set = set()

    old_tax_id_set = set()
    new_tax_id_set = set()

    for r in rows:
        column_count = len(r)
        column_count_set.add(column_count)

        old_tax_id = r[0]
        new_tax_id = r[1]

        old_tax_id_set.add(old_tax_id)
        new_tax_id_set.add(new_tax_id)

    print('row_count', row_count)
    print('column_count_set', column_count_set)

    print('len(old_tax_id_set)', len(old_tax_id_set))
    print('len(new_tax_id_set)', len(new_tax_id_set))

    return_dict = dict()

    return_dict['column_count_set'] = column_count_set
    return_dict['old_tax_id_set'] = old_tax_id_set
    return_dict['new_tax_id_set'] = new_tax_id_set
    return_dict['row_count'] = row_count

    return return_dict


def _examine_ncbi_taxonomy_gencode_dump(file_path):
    rows = _parse_ncbi_taxonomy_dump_file(file_path=file_path)

    row_count = len(rows)

    column_count_set = set()

    genetic_code_id_set = set()
    abbreviation_set = set()
    name_set = set()
    translation_table_set = set()
    start_codons_set = set()

    codon_count_set = set()

    for r in rows:
        column_count = len(r)
        column_count_set.add(column_count)

        genetic_code_id = r[0]
        abbreviation = r[1]
        name = r[2]
        translation_table = r[3].strip()
        start_codons = r[4].strip()

        # print([genetic_code_id, abbreviation, name, translation_table,
        #        start_codons])

        genetic_code_id_set.add(genetic_code_id)
        abbreviation_set.add(abbreviation)
        name_set.add(name)
        translation_table_set.add(translation_table)
        start_codons_set.add(start_codons)

        codon_count_set.add(len(translation_table))
        codon_count_set.add(len(start_codons))

    print('row_count', row_count)
    print('column_count_set', column_count_set)

    print('len(genetic_code_id_set)', len(genetic_code_id_set))
    print('len(abbreviation_set)', len(abbreviation_set))
    print('len(name_set)', len(name_set))
    print('len(translation_table_set)', len(translation_table_set))
    print('len(start_codons_set)', len(start_codons_set))

    print('len(codon_count_set)', len(codon_count_set))
    print('codon_count_set', codon_count_set)

    return_dict = dict()

    return_dict['column_count_set'] = column_count_set
    return_dict['genetic_code_id_set'] = genetic_code_id_set
    return_dict['name_set'] = name_set
    return_dict['codon_count_set'] = codon_count_set
    return_dict['row_count'] = row_count

    return return_dict
