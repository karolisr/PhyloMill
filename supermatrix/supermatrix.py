# !/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import unicode_literals

import sys
import os

from PyQt4 import QtCore
from PyQt4 import QtGui

PS = os.path.sep
RESOURCES_DIR = 'resources'
ICONS_DIR = RESOURCES_DIR + PS + 'icons'

class KRWorkflowStep(QtGui.QGroupBox):

    def __init__(self, title='title'):

        super(KRWorkflowStep, self).__init__(title)

        self.init_ui()

    def init_ui(self):

        # name_label = QtGui.QLabel(name)

        configure_button = QtGui.QPushButton("Configure")
        run_button = QtGui.QPushButton("Run")

        hbox = QtGui.QHBoxLayout()
        # hbox.addWidget(name_label)
        # hbox.addSpacing(0)
        hbox.addStretch(0)
        hbox.addWidget(configure_button)
        hbox.addWidget(run_button)

        self.setLayout(hbox)

class KRMainWindow(QtGui.QMainWindow):

    def __init__(self, name='KRMainWindow', width=500, height=300):

        super(KRMainWindow, self).__init__()

        self.name = name
        self.init_ui(width, height)

    def init_ui(self, width, height):

        step_1 = KRWorkflowStep(title='1. Search && Download')
        step_2 = KRWorkflowStep(title='2. Filter')

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(step_1)
        vbox.addWidget(step_2)
        vbox.addStretch(0)

        self.widget = QtGui.QWidget()
        self.widget.setLayout(vbox)
        self.setCentralWidget(self.widget)

        self.create_actions()
        self.create_menus()

        self.resize(width, height)

        self.setWindowTitle(self.name)

        self.show()

        self.center_window()

    def open(self):

        chosen = QtGui.QFileDialog(self).getExistingDirectory(self)
        print(chosen)
        return(chosen)

    def save(self):

        chosen = QtGui.QFileDialog(self).getExistingDirectory(self)
        print(chosen)
        return(chosen)

    def about(self):

        QtGui.QMessageBox.about(self, "About",
            "Supermatrix by Karolis Ramanauskas")

    def create_actions(self):

        self.open_action = QtGui.QAction("&Open...", self,
                shortcut=QtGui.QKeySequence.Open,
                statusTip="Open an existing configuration directory",
                triggered=self.open)

        self.save_action = QtGui.QAction("&Save", self,
                shortcut=QtGui.QKeySequence.Save,
                statusTip="Save the configuration to a directory",
                triggered=self.save)

        self.exit_action = QtGui.QAction("E&xit", self,
                shortcut=QtGui.QKeySequence.Quit,
                statusTip="Exit the application",
                triggered=self.close)

        self.about_action = QtGui.QAction("&About", self,
                statusTip="Show the application's About box",
                triggered=self.about)

        self.about_qt_action = QtGui.QAction("About &Qt", self,
                statusTip="Show the Qt library's About box",
                triggered=QtGui.qApp.aboutQt)

    def create_menus(self):

        self.fileMenu = self.menuBar().addMenu("&File")
        self.fileMenu.addAction(self.open_action)
        self.fileMenu.addAction(self.save_action)
        self.fileMenu.addSeparator();
        self.fileMenu.addAction(self.exit_action)

        self.helpMenu = self.menuBar().addMenu("&Help")
        self.helpMenu.addAction(self.about_action)
        self.helpMenu.addAction(self.about_qt_action)

    def center_window(self):

        self_geometry = self.frameGeometry()
        desktop_center = QtGui.QDesktopWidget().availableGeometry().center()
        self_geometry.moveCenter(desktop_center)
        self.move(self_geometry.topLeft())

def main():

    application = QtGui.QApplication(sys.argv)
    main_window = KRMainWindow(name='Supermatrix', width=0, height=0)
    sys.exit(application.exec_())

if __name__ == '__main__':

    main()


# main_tab_widget = QtGui.QToolBox()
# search_tab_widget = QtGui.QWidget()
# main_tab_widget.addItem(search_tab_widget, 'Search')
# settings_tab_widget = QtGui.QWidget()
# main_tab_widget.addItem(settings_tab_widget, 'Settings')
# main_tab_widget = QtGui.QTabWidget()
# self.setCentralWidget(main_tab_widget)
# main_tab_widget.setTabShape(QtGui.QTabWidget.Rounded)
# main_tab_widget.setTabPosition(QtGui.QTabWidget.North)
# search_tab_widget = QtGui.QWidget()
# main_tab_widget.addTab(search_tab_widget, 'Search')
# settings_tab_widget = QtGui.QWidget()
# main_tab_widget.addTab(settings_tab_widget, 'Settings')
# status_bar = self.statusBar()
# exit_action_icon = ICONS_DIR + PS + 'Switch_32.png'
# exit_action = QtGui.QAction(QtGui.QIcon(exit_action_icon), 'Exit', self)
# exit_action.triggered.connect(self.close)
# menu_bar = self.menuBar()
# file_menu = menu_bar.addMenu('&File')
# file_menu.addAction(exit_action)
# self.setUnifiedTitleAndToolBarOnMac(True)
# print(self.unifiedTitleAndToolBarOnMac())
# tool_bar = QtGui.QToolBar()
# self.addToolBar(tool_bar)
# tool_bar = self.addToolBar('Main Toolbar')
# tool_bar_icon_size = QtCore.QSize(32, 32)
# tool_bar.setIconSize(tool_bar_icon_size)
# tool_bar.setMovable(False)
# tool_bar.addAction(exit_action)
# self.widget = QtGui.QWidget()
# self.setCentralWidget(self.widget)
