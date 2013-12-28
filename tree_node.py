# -*- coding: utf-8 -*-

'''
Created on Dec 27, 2013
@author: Karolis Ramanauskas
@copyright: 2013 Karolis Ramanauskas. All rights reserved.
@license: GPLv3
'''

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

__all__ = []
__version__ = 0.1
__updated__ = '2013-12-27'


class TreeNode(object):
    '''
    documentation
    '''

    def __init__(self, name, parent=None):
        '''
        documentation
        '''

        self._name = name
        self._children = []
        self._parent = parent

        if parent is not None:
            parent.add_child(self)

#     def typeInfo(self):
#         return "NODE"

    def add_child(self, child):
        '''
        documentation
        '''
        self._children.append(child)

#     def insertChild(self, position, child):
#
#         if position < 0 or position > len(self._children):
#             return False
#
#         self._children.insert(position, child)
#         child._parent = self
#         return True

#     def removeChild(self, position):
#
#         if position < 0 or position > len(self._children):
#             return False
#
#         child = self._children.pop(position)
#         child._parent = None
#
#         return True

    def name(self):
        '''
        documentation
        '''
        return self._name

    def set_name(self, name):
        '''
        documentation
        '''
        self._name = name

    def child(self, row):
        '''
        documentation
        '''
        return self._children[row]

    def children_count(self):
        '''
        documentation
        '''
        return len(self._children)

    def parent(self):
        '''
        documentation
        '''
        return self._parent

    def row(self):
        '''
        documentation
        '''
        if self._parent is not None:
            return self._parent._children.index(self)

    def log(self, tab_level=-1):
        '''
        documentation
        '''

        output = ""
        tab_level += 1

        for current_tab_level in range(tab_level):
            output += "\t"

        output += "|------" + self._name + "\n"

        for child in self._children:
            output += child.log(tab_level)

        tab_level -= 1
        output += "\n"

        return output

    def __repr__(self):

        return self.log()

if __name__ == '__main__':
    pass
