# -*- coding: utf-8 -*-

"""

This module fascilitates interaction with NCBI databases outside of
ENTREZ system.

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import with_statement

import os.path
from os import remove
import hashlib
import zipfile

from krpy import PS
from krpy import FILE_OPEN_MODE_READ

from krpy import Root
from krpy import Error
from krpy import internet
from krpy import io

TAX_BASE_URL = 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/'

# Files expected to be in taxdmp.zip MD5 file should be first
TAXDMP_FILES = ['taxdmp.zip.md5', 'readme.txt', 'nodes.dmp',
                'names.dmp', 'merged.dmp', 'gencode.dmp', 'gc.prt',
                'division.dmp', 'delnodes.dmp', 'citations.dmp']
TAXDMP_ARCHIVE = 'taxdmp.zip'

# Files expected to be in taxcat.zip MD5 file should be first
TAXCAT_FILES = ['taxcat.zip.md5', 'categories.dmp']
TAXCAT_ARCHIVE = 'taxcat.zip'


NAME_CLASS_SET = set([
    'acronym', 'teleomorph', 'scientific name', 'synonym', 'blast name',
    'misspelling', 'in-part', 'includes', 'genbank acronym', 'misnomer',
    'common name', 'equivalent name', 'type material', 'authority',
    'genbank common name', 'genbank synonym', 'anamorph', 'genbank anamorph'])


def _extract_md5_hash(file_path):
    md5_reported = None
    with open(file_path, FILE_OPEN_MODE_READ) as f_md5:
        line = f_md5.readline()
        md5_reported = line.split(' ')[0]
    return md5_reported


def _download_ncbi_taxonomy_data(directory_path,
                                 archive_url,
                                 md5_url,
                                 archive_path,
                                 md5_path):

    internet.urlretrieve(url=archive_url, filename=archive_path)
    internet.urlretrieve(url=md5_url, filename=md5_path)

    md5_reported = _extract_md5_hash(file_path=md5_path)
    print('md5_reported:', md5_reported)
    md5_actual = io.generate_md5_hash_for_file(file_path=archive_path)
    print('  md5_actual:', md5_actual)

    if md5_reported != md5_actual:
        message = (
            'MD5 hash of the file {f} does not match the reported '
            'hash in file {r}')
        message = message.format(f=archive_path, r=md5_path)
        raise Error(message)
    else:
        z = zipfile.ZipFile(file=archive_path, mode='r')
        z.extractall(path=directory_path)
        return True


def _update_ncbi_taxonomy_data(taxdmp_path, taxcat_path,
                               force_redownload=False,
                               check_for_updates=True):

    download_taxdmp = False
    download_taxcat = False

    taxdmp_archive_path = taxdmp_path + TAXDMP_ARCHIVE
    taxcat_archive_path = taxcat_path + TAXCAT_ARCHIVE

    taxdmp_md5_url = TAX_BASE_URL + TAXDMP_FILES[0]
    taxcat_md5_url = TAX_BASE_URL + TAXCAT_FILES[0]

    taxdmp_md5_path = taxdmp_path + TAXDMP_FILES[0]
    taxcat_md5_path = taxcat_path + TAXCAT_FILES[0]

    taxdmp_md5_path_new = taxdmp_md5_path + '_new'
    taxcat_md5_path_new = taxcat_md5_path + '_new'

    if not force_redownload:

        for file_name in TAXDMP_FILES:
            file_path = taxdmp_path + file_name
            if not os.path.exists(file_path):
                download_taxdmp = True
                break

        for file_name in TAXCAT_FILES:
            file_path = taxcat_path + file_name
            if not os.path.exists(file_path):
                download_taxcat = True
                break

        if check_for_updates:
            if not download_taxdmp:
                internet.urlretrieve(
                    url=taxdmp_md5_url,
                    filename=taxdmp_md5_path_new)

                old_md5 = _extract_md5_hash(file_path=taxdmp_md5_path)
                new_md5 = _extract_md5_hash(file_path=taxdmp_md5_path_new)

                remove(taxdmp_md5_path_new)

                print('old_md5:', old_md5)
                print('new_md5:', new_md5)

                if old_md5 != new_md5:
                    download_taxdmp = True

            if not download_taxcat:
                internet.urlretrieve(
                    url=taxcat_md5_url,
                    filename=taxcat_md5_path_new)

                old_md5 = _extract_md5_hash(file_path=taxcat_md5_path)
                new_md5 = _extract_md5_hash(file_path=taxcat_md5_path_new)

                remove(taxcat_md5_path_new)

                print('old_md5:', old_md5)
                print('new_md5:', new_md5)

                if old_md5 != new_md5:
                    download_taxcat = True

    else:
        download_taxdmp = True
        download_taxcat = True

    # print('download_taxdmp', download_taxdmp)
    # print('download_taxcat', download_taxcat)

    if download_taxdmp:
        _download_ncbi_taxonomy_data(
            directory_path=taxdmp_path,
            archive_url=TAX_BASE_URL+TAXDMP_ARCHIVE,
            md5_url=taxdmp_md5_url,
            archive_path=taxdmp_archive_path,
            md5_path=taxdmp_md5_path)

        remove(taxdmp_archive_path)

    if download_taxcat:
        _download_ncbi_taxonomy_data(
            directory_path=taxcat_path,
            archive_url=TAX_BASE_URL+TAXCAT_ARCHIVE,
            md5_url=taxcat_md5_url,
            archive_path=taxcat_archive_path,
            md5_path=taxcat_md5_path)

        remove(taxcat_archive_path)

    return [download_taxdmp, download_taxcat]


def _parse_ncbi_taxonomy_dump_file(file_path):

    field_terminator = '\t|\t'
    row_terminator = '\t|\n'

    f = open(file_path, FILE_OPEN_MODE_READ)

    lines = list()

    for l in f:
        l = l.strip(row_terminator)
        l = l.split(field_terminator)
        lines.append(l)

    f.close()

    return lines


def _parse_codons(tax_gencode_prt_path):

    base1 = ''
    base2 = ''
    base3 = ''

    base1_start = '  -- Base1  '
    base2_start = '  -- Base2  '
    base3_start = '  -- Base3  '

    with open(tax_gencode_prt_path, FILE_OPEN_MODE_READ) as f:
        for l in f:
            l = l.strip('\n')

            if l.startswith(base1_start):
                base1 = l.strip(base1_start)
            elif l.startswith(base2_start):
                base2 = l.strip(base2_start)
            elif l.startswith(base3_start):
                base3 = l.strip(base3_start)
                break

    return_value = list(zip(base1, base2, base3))

    return return_value


def _parse_names_dump(file_path):

    rows = _parse_ncbi_taxonomy_dump_file(file_path=file_path)
    tax_id_keyed_dict = dict()

    for r in rows:

        tax_id = r[0]
        name = r[1]
        unique_name = r[2]
        name_class = r[3]

        if tax_id not in tax_id_keyed_dict:
            tax_id_keyed_dict[tax_id] = list()

        dict_entry = {
            'name': name,
            'unique_name': unique_name,
            'name_class': name_class}

        tax_id_keyed_dict[tax_id].append(dict_entry)

    return tax_id_keyed_dict


def _parse_delnodes_dump(file_path):

    rows = _parse_ncbi_taxonomy_dump_file(file_path=file_path)
    tax_id_set = set()

    for r in rows:

        tax_id = r[0]
        tax_id_set.add(tax_id)

    return tax_id_set


def _parse_merged_dump(file_path):

    rows = _parse_ncbi_taxonomy_dump_file(file_path=file_path)
    new_to_old_tax_id_mapping_dict = dict()

    for r in rows:

        old_tax_id = r[0]
        new_tax_id = r[1]

        new_to_old_tax_id_mapping_dict[old_tax_id] = new_tax_id

    return new_to_old_tax_id_mapping_dict


def _parse_nodes_dump(file_path):

    rows = _parse_ncbi_taxonomy_dump_file(file_path=file_path)

    child_to_parent_tax_id_mapping_dict = dict()
    taxid_rank_dict = dict()
    taxid_genetic_code_id_dict = dict()
    taxid_mitochondrial_genetic_code_id_dict = dict()

    for r in rows:

        tax_id = r[0]
        parent_tax_id = r[1]
        rank = r[2]

        # embl_code = r[3]
        # division_id = r[4]
        # inherited_div_flag = r[5]
        genetic_code_id = r[6]
        # inherited_GC_flag = r[7]
        mitochondrial_genetic_code_id = r[8]
        # inherited_MGC_flag = r[9]
        # GenBank_hidden_flag = r[10]
        # hidden_subtree_root_flag = r[11]

        # not every row contains this column, and it is not needed
        # comments = r[12]

        child_to_parent_tax_id_mapping_dict[tax_id] = parent_tax_id
        taxid_rank_dict[tax_id] = rank
        taxid_genetic_code_id_dict[tax_id] = genetic_code_id
        taxid_mitochondrial_genetic_code_id_dict[tax_id] = \
            mitochondrial_genetic_code_id

    return_value = [
        child_to_parent_tax_id_mapping_dict,
        taxid_rank_dict,
        taxid_genetic_code_id_dict,
        taxid_mitochondrial_genetic_code_id_dict]

    return return_value


def _parse_gencode_dump(file_path):

    rows = _parse_ncbi_taxonomy_dump_file(file_path=file_path)

    genetic_code_id_to_name_dict = dict()
    genetic_code_id_to_translation_table_dict = dict()
    genetic_code_id_to_start_codons_dict = dict()

    for r in rows:

        genetic_code_id = r[0]
        # abbreviation = r[1]
        name = r[2]
        translation_table = r[3].strip()
        start_codons = r[4].strip()

        genetic_code_id_to_name_dict[genetic_code_id] = name
        genetic_code_id_to_translation_table_dict[genetic_code_id] = \
            translation_table
        genetic_code_id_to_start_codons_dict[genetic_code_id] = start_codons

    return_value = [
        genetic_code_id_to_name_dict,
        genetic_code_id_to_translation_table_dict,
        genetic_code_id_to_start_codons_dict]

    return return_value
