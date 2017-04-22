#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import with_statement

from krpy import PYTHON_VERSION_STRING


def main():

    py_ver_msg = '\nPython version: {pv}\n'.format(pv=PYTHON_VERSION_STRING)
    print(py_ver_msg)

    from krpy import entrez

    return_value = entrez.esearch(db='taxonomy', term='Brassicaceae')
    id_list = return_value['IdList']

    term1 = str('txid' + id_list[0] + '[Organism]' +
                ' AND mitochondrion[filter]')

    term2 = str('txid' + id_list[0] + '[Organism]' +
                ' AND chloroplast[filter]')

    term3 = str('txid' + id_list[0] + '[Organism]' +
                ' AND refseq[filter]')

    term4 = str('txid' + id_list[0] + '[Organism]' +
                ' NOT mitochondrion[filter] NOT chloroplast[filter] NOT refseq[filter]')

    print(term1)
    return_value = entrez.esearch(db='nuccore', term=term1)
    id_list = return_value['IdList']
    print(len(id_list))

    print(term2)
    return_value = entrez.esearch(db='nuccore', term=term2)
    id_list = return_value['IdList']
    print(len(id_list))

    print(term3)
    return_value = entrez.esearch(db='nuccore', term=term3)
    id_list = return_value['IdList']
    print(len(id_list))

    print(term4)
    return_value = entrez.esearch(db='nuccore', term=term4)
    id_list = return_value['IdList']
    print(len(id_list))

if __name__ == '__main__':
    main()
