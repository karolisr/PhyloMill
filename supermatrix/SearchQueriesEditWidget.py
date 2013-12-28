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

__updated__ = '2013-12-28'

from PyQt4 import QtGui
from PyQt4 import QtCore


class SearchQueriesTreeModel(QtCore.QAbstractItemModel):
    '''
    SearchQueriesTreeModel
    '''

    def __init__(self, root_node, parent=None):
        '''
        Constructor
        '''
        super(SearchQueriesTreeModel, self).__init__(parent)
        self._root_node = root_node

    def rowCount(self, parent_q_model_index):
        '''
        Returns number of children for a given QModelIndex
        '''
        if not parent_q_model_index.isValid():
            parent_node = self._root_node
        else:
            parent_node = parent_q_model_index.internalPointer()
        return parent_node.children_count()

    def columnCount(self, parent_q_model_index):
        '''
        columnCount
        '''
        return 1

    def data(self, q_model_index, role):
        '''
        Returns display data to the view based on QModelIndex and role.
        '''
        if not q_model_index.isValid():
            return None
        node = q_model_index.internalPointer()
        if role == QtCore.Qt.DisplayRole or role == QtCore.Qt.EditRole:
            if q_model_index.column() == 0:
                return node.name()

#         if role == QtCore.Qt.DecorationRole:
#             if index.column() == 0:
#                 typeInfo = node.typeInfo()
#
#                 if typeInfo == "LIGHT":
#                     return QtGui.QIcon(QtGui.QPixmap(":/Light.png"))
#
#                 if typeInfo == "TRANSFORM":
#                     return QtGui.QIcon(QtGui.QPixmap(":/Transform.png"))
#
#                 if typeInfo == "CAMERA":
#                     return QtGui.QIcon(QtGui.QPixmap(":/Camera.png"))

    def setData(self, index, value, role=QtCore.Qt.EditRole):
        '''
        setData
        '''
        if index.isValid():
            if role == QtCore.Qt.EditRole:
                node = index.internalPointer()
                node.set_name(value)
                return True
        return False

    def headerData(self, section, orientation, role):
        '''
        headerData
        '''
        if role == QtCore.Qt.DisplayRole:
            return 'Search Queries'

    def flags(self, q_model_index):
        '''
        flags
        '''
        return QtCore.Qt.ItemIsEnabled | \
            QtCore.Qt.ItemIsSelectable | \
            QtCore.Qt.ItemIsEditable

    def parent(self, q_model_index):
        '''
        Returns the parent of the node at the given QModelIndex.
        '''
        node = self.get_node(q_model_index)
        parent_node = node.parent()
        if parent_node == self._root_node:
            return QtCore.QModelIndex()
        return self.createIndex(parent_node.index(), 0, parent_node)

    def index(self, row, column, parent_q_model_index):
        '''
        Returns a QModelIndex that corresponds to the given row, column and
        parent QModelIndex.
        '''
        parent_node = self.get_node(parent_q_model_index)
        child = parent_node.child(row)
        if child:
            return self.createIndex(row, column, child)
        else:
            return QtCore.QModelIndex()

    def get_node(self, q_model_index):
        '''
        Returns a Node that corresponds to the given QModelIndex.
        '''
        if q_model_index.isValid():
            node = q_model_index.internalPointer()
            if node:
                return node
        return self._root_node


#     """INPUTS: int, int, QModelIndex"""
#     def insertRows(self, position, rows, parent=QtCore.QModelIndex()):
#
#         parentNode = self.getNode(parent)
#
#         self.beginInsertRows(parent, position, position + rows - 1)
#
#         for row in range(rows):
#
#             childCount = parentNode.childCount()
#             childNode = Node("untitled" + str(childCount))
#             success = parentNode.insertChild(position, childNode)
#
#         self.endInsertRows()
#
#         return success
#
#     def insertLights(self, position, rows, parent=QtCore.QModelIndex()):
#
#         parentNode = self.getNode(parent)
#
#         self.beginInsertRows(parent, position, position + rows - 1)
#
#         for row in range(rows):
#
#             childCount = parentNode.childCount()
#             childNode = LightNode("light" + str(childCount))
#             success = parentNode.insertChild(position, childNode)
#
#         self.endInsertRows()
#
#         return success
#
#     """INPUTS: int, int, QModelIndex"""
#     def removeRows(self, position, rows, parent=QtCore.QModelIndex()):
#
#         parentNode = self.getNode(parent)
#         self.beginRemoveRows(parent, position, position + rows - 1)
#
#         for row in range(rows):
#             success = parentNode.removeChild(position)
#
#         self.endRemoveRows()
#
#         return success


class SearchQueriesTreeView(QtGui.QTreeView):
    '''
    SearchQueriesTreeView
    '''

    def __init__(self):
        '''
        Constructor
        '''
        super(SearchQueriesTreeView, self).__init__()


class SearchQueriesEditWidget(QtGui.QWidget):
    '''
    SearchQueriesEditWidget
    '''

    def __init__(self):
        '''
        Constructor
        '''
        super(SearchQueriesEditWidget, self).__init__()

if __name__ == '__main__':

    import sys

    import krpy

    from krpy import supermatrix

    if krpy.debug.RUN_DEBUG_CODE:

        with open('test_data/search_queries.tsv', 'r') as\
            SEARCH_QUERIES_HANDLE:
            SEARCH_QUERIES = supermatrix.io.parse_search_queries_file(
                                                        SEARCH_QUERIES_HANDLE)

        APP = QtGui.QApplication(sys.argv)

        MODEL = SearchQueriesTreeModel(SEARCH_QUERIES)

        TREE_VIEW = QtGui.QTreeView()
        TREE_VIEW.setModel(MODEL)
        TREE_VIEW.show()
        TREE_VIEW.header().show()

        sys.exit(APP.exec_())
