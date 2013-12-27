'''
    Created on Dec 26, 2013
    @author: Karolis Ramanauskas
    @copyright: 2013 Karolis Ramanauskas. All rights reserved.
    @license: GPLv3
'''

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

__all__ = []
__version__ = 0.1
__date__ = '2013-12-26'
__updated__ = '2013-12-26'


def attr(obj):
    for prop, value in vars(obj).iteritems():
        print(prop, " : ", value)
        print('---')


def random_id(length):
    import os
    return str(abs(hash(str(os.urandom(200)))))[0:length]


def parse_line(line, delimiter, quotechar, stripchar):

    '''
        Will parse a line separated by delimiter and will return a list. Parts
        surrounded by quotechar will be grouped even if delimiter is found
        within.
    '''

    if not isinstance(line, basestring):
        return ''

    if not isinstance(delimiter, basestring):
        return ''

    # if not isinstance(quotechar, basestring):
    #    return ''

    # print(repr(line))

    nl = '\n'
    l_rem_nl = line.split(nl)[0]

    if stripchar:
        l_rem_nl = l_rem_nl.replace(stripchar, '')

    # print(repr(l_rem_nl))
    l_spl_q = [l_rem_nl]
    if quotechar:
        l_spl_q = l_rem_nl.split(quotechar)

        # print(l_spl_q)

        if l_rem_nl.startswith(quotechar):
            l_spl_q.remove('')
        if l_rem_nl.endswith(quotechar):
            l_spl_q.pop()

    # print(l_spl_q)

    final_list = list()

    for x in l_spl_q:
        if x.startswith(delimiter) or len(l_spl_q) == 1:
            if x is delimiter:
                continue
            else:
                x_spl_d = x.split(delimiter)
                if x.startswith(delimiter):
                    x_spl_d.remove('')
                if x.endswith(delimiter):
                    x_spl_d.pop()
                for y in x_spl_d:
                    y = y.strip()
                    final_list.append(y)
        else:
            final_list.append(x)

    # print(final_list)
    return final_list


def timestamp():
    import datetime
    return str(datetime.datetime.now()).split('.')[0]
