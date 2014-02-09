#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
DOCSTRING
'''

from __future__ import print_function, unicode_literals
from PyQt4 import QtGui


class KRWorkflowStepWidget(QtGui.QGroupBox):

    '''
    DOCSTRING
    '''

    def __init__(self, title, configure_slot, run_slot):

        super(KRWorkflowStepWidget, self).__init__(title)
        self.init_ui(configure_slot, run_slot)

    def init_ui(self, configure_slot, run_slot):

        '''
        DOCSTRING
        '''

        status_label = QtGui.QLabel('Status')

        if configure_slot:
            configure_button = QtGui.QPushButton("Configure")
            configure_button.setFixedSize(100, 32)
            configure_button.clicked.connect(configure_slot)

        if run_slot:
            run_button = QtGui.QPushButton("Run")
            run_button.setFixedSize(100, 32)
            configure_button.clicked.connect(run_slot)

        h_box = QtGui.QHBoxLayout()

        h_box.addWidget(status_label)

        h_box.addSpacing(10)
        h_box.addStretch(0)

        if configure_slot:
            h_box.addWidget(configure_button)

        if run_slot:
            h_box.addWidget(run_button)

        self.setLayout(h_box)


def test_widget():

    '''
    DOCSTRING
    '''

    import sys

    application = QtGui.QApplication(sys.argv)
    kr_workflow_step_widget = KRWorkflowStepWidget(
                                                title='KRWorkflowStepWidget',
                                                configure_slot=None,
                                                run_slot=None)
    kr_workflow_step_widget.show()
    sys.exit(application.exec_())

if __name__ == '__main__':

    test_widget()
