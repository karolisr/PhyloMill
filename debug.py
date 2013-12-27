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
__updated__ = '2013-12-27'

RUN_DEBUG_CODE = True
SILENCE_DEBUG_MESSAGES = False


def message(msg, sender):

    '''
        Used for simple and unified console logging.
    '''
    if not SILENCE_DEBUG_MESSAGES:
        print(str(sender) + ': ' + str(msg))

if __name__ == '__main__':
    pass
