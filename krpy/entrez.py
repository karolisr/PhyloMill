# -*- coding: utf-8 -*-

"""

This module fascilitates interaction with NCBI's Entrez Programming
Utilities (E-utilities).

More information on E-utilities at:
    http://www.ncbi.nlm.nih.gov/books/NBK25497

Database names and unique identifiers returned can be found here:
    http://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import with_statement

import sys
import locale
from xml.etree import ElementTree

from krpy import Error
from krpy import internet

ENTREZ_BASE_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
ENCODING = locale.getdefaultlocale()[1]


def esearch(db, term):

    """

    Wraps ESearch E-utility.

    :param db: Name of the Entrez database to search.
    :type db: str

    :param term: Search terms.
    :type term: str

    :returns: A dictionary with these keys: Database, Count, IdList, QueryKey,
        WebEnv
    :rtype: dict

    """

    eutil = 'esearch.fcgi'

    term = internet.quote_plus(term)

    url = ENTREZ_BASE_URL + eutil + '?db=' + db + '&term=' + term \
        + '&retmax=100000' + '&usehistory=y'

    response = internet.urlopen(url)
    tree = ElementTree.parse(response)
    root = tree.getroot()
    response.close()

    count = None
    query_key = None
    web_env = None
    id_set = set()

    for child in root:
        if child.tag == 'IdList':
            for id in child:
                id_set.add(id.text)
        if child.tag == 'QueryKey':
            query_key = child.text
        if child.tag == 'WebEnv':
            web_env = child.text

    id_list = list(id_set)
    id_list.sort()

    count = len(id_list)
    if count >= 99999:
        message = (
            'There are more than 99,999 unique identifiers: {c}.')
        message = message.format(c=count)
        raise Error(message)

    return_value = {
        'Database': db,
        'Count': count,
        'IdList': id_list,
        'QueryKey': query_key,
        'WebEnv': web_env}

    return return_value


def epost(db, id_list):

    """

    Wraps EPost E-utility.

    :param db: Name of the Entrez database to search.
    :type db: str

    :param id_list: List of unique identifiers.
    :type id_list: list

    :returns: A dictionary with these keys: Database, Count, QueryKey, WebEnv
    :rtype: dict

    """

    eutil = 'epost.fcgi'

    id_list_string = ','.join(id_list)

    url = ENTREZ_BASE_URL + eutil + '?db=' + db + '&id=' + id_list_string

    response = internet.urlopen(url)
    tree = ElementTree.parse(response)
    root = tree.getroot()
    response.close()

    query_key = None
    web_env = None

    for child in root:
        if child.tag == 'QueryKey':
            query_key = child.text
        if child.tag == 'WebEnv':
            web_env = child.text

    count = len(id_list)

    return_value = {
        'Database': db,
        'Count': count,
        'QueryKey': query_key,
        'WebEnv': web_env}

    return return_value


def efetch(
    data,
    ret_type,
    ret_mode,
    parser,
    return_string=False
):

    """

    Wraps EFetch E-utility.

    :param data: A dictionary returned by :func:`esearch` or :func:`epost` with
        these keys: Database, Count, QueryKey, WebEnv
    :type data: dict

    :param ret_type: Retrieval type. This parameter specifies the record view
        returned, such as Abstract or MEDLINE from PubMed, or GenPept or FASTA
        from protein.
    :type ret_type: str

    :param ret_mode: Retrieval mode. This parameter specifies the data format
        of the records returned, such as plain text, HMTL or XML.
    :type ret_mode: str

    See the link below for the possible values of ret_type and ret_mode:
        http://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly

    :param parser: A function that will be called to interpret downloaded data.
        This function may be called several times as :func:`efetch` downloads
        downloads data in batches.
    :type parser: function

    :param return_string: If ``False`` (default), a file-like object produced
        by :func:`urllib.request.urlopen` will be passed to parser. If
        ``True``, a decoded string will be passed instead.
    :type return_string: bool

    :returns: A list of one or more items which will be of the type produced by
        the parser.
    :rtype: list

    """

    eutil = 'efetch.fcgi'

    db = data['Database']
    query_key = data['QueryKey']
    web_env = data['WebEnv']
    count = data['Count']

    retmax = 500
    retstart = 0

    return_value = list()

    for retstart in range(0, count, retmax):

        url = ENTREZ_BASE_URL + eutil + '?db=' + db \
            + '&query_key=' + query_key + '&WebEnv=' + web_env \
            + '&retstart=' + str(retstart) + '&retmax=' + str(retmax) \
            + '&rettype=' + ret_type + '&retmode=' + ret_mode

        response = internet.urlopen(url)

        parsed = None

        if return_string:
            response_contents = response.read().decode(ENCODING)
            parsed = parser(response_contents)
        else:
            parsed = parser(response)

        return_value.append(parsed)

        response.close()

    return return_value
