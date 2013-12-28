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
__updated__ = '2013-12-28'


class Node(object):
    '''
    A representation of a node in a tree data structure.
    '''

    def __init__(self, name, data=None, parent=None):
        '''
        Constructor.
        '''

        self._name = name
        self._data = data
        self._parent = parent
        self._children = list()

        if parent is not None:
            parent.add_child(self)

    def add_child(self, child):
        '''
        documentation
        TODO: make sure all children names are unique.
        '''
        self._children.append(child)
        child._parent = self
        return child

    def insert_child(self, position, child):
        '''
        documentation
        TODO: make sure all children names are unique.
        '''
        if position < 0 or position > len(self._children):
            return False

        self._children.insert(position, child)
        child._parent = self
        return child

    def remove_child(self, position):
        '''
        documentation
        '''
        if position < 0 or position > len(self._children):
            return False

        child = self._children.pop(position)
        child._parent = None

        return True

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

    def data(self):
        '''
        documentation
        '''
        return self._data

    def set_data(self, data):
        '''
        documentation
        '''
        self._data = data

    def children_count(self):
        '''
        documentation
        '''
        return len(self._children)

    def children(self):
        '''
        documentation
        '''
        return self._children

    def children_names(self):
        '''
        documentation
        '''
        return [x.name() for x in self._children]

    def child_by_index(self, index):
        '''
        documentation
        '''
        return self._children[index]

    def child_by_name(self, name):
        '''
        documentation
        '''
        return [x for x in self._children if x.name() == name][0]

    def child(self, identifier):
        '''
        documentation
        '''
        if isinstance(identifier, basestring):
            return self.child_by_name(identifier)
        elif isinstance(identifier, int):
            return self.child_by_index(identifier)

    def parent(self):
        '''
        documentation
        '''
        return self._parent

    def index(self):
        '''
        documentation
        '''
        if self._parent is not None:
            return self._parent._children.index(self)

    def _representation(self, level=0):
        '''
        documentation
        '''

        output = ''
        output += '    ' * level

        if self._parent is None:
            output += '─'
        else:
            output += '└'

        output += self._name + ', data: ' + str(self._data) + '\n'

        level += 1

        for child in self._children:
            output += child._representation(level)

        level -= 1

        return output

    def __repr__(self):
        return self._representation()

if __name__ == '__main__':
    import krpy
    if krpy.debug.RUN_DEBUG_CODE:
        ROOT = Node('root')
        ROOT.add_child(Node('child1'))
        CHILD2 = ROOT.add_child(Node('child2'))
        CHILD2.add_child(Node('child3'))
        CHILD4 = ROOT.add_child(Node('child4'))
        CHILD4.add_child(Node('child5'))
        CHILD6 = CHILD4.add_child(Node('child6'))
        CHILD6.add_child(Node('child7'))
        print(ROOT)
