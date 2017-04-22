# -*- coding: utf-8 -*-

"""

This module fascilitates interaction with the internet.

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import with_statement

import sys

# Make sure this works with Python 3.X and 2.7
if sys.hexversion < 0x03000000:
    from urllib import urlretrieve as urlretrieve_kr
    from urllib2 import Request as Request_kr
    from urllib2 import urlopen as urlopen_kr
    from urllib import quote_plus as quote_plus_kr
    # from urllib2 import URLError
    # from urllib2 import HTTPError
else:
    from urllib.request import urlretrieve as urlretrieve_kr
    from urllib.request import Request as Request_kr
    from urllib.request import urlopen as urlopen_kr
    from urllib.parse import quote_plus as quote_plus_kr
    # from urllib.error import URLError
    # from urllib.error import HTTPError

import socket

# On Python 2.7.9 on Mac OSX Yosemite 10.10.1 without this urlretrieve
# is very slow
# socket.setdefaulttimeout(5.0)
# print(socket.getdefaulttimeout())


def quote_plus(value):
    return quote_plus_kr(value)


def urlopen(url):

    request = Request_kr(url=url)
    # ToDo: deal with URLError/HTTPError
    response = urlopen_kr(request)

    return response


def urlretrieve(url, filename):
    # ToDo: deal with URLError/HTTPError
    return urlretrieve_kr(url=url, filename=filename)
