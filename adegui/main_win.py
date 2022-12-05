from PyQt5.QtWidgets import (QMainWindow, QVBoxLayout, QWidget,
                             qApp, QStyle)
from PyQt5.QtCore import pyqtSlot
from adegui.work_area.work_area import WorkAreaTabs
from adegui.main_buttons import MainButtons
from adegui.generate_script import write_ade_script_from_config


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("autodE GUI")

        # The layout for the main window => buttons at bottom,
        # tabs and work area at top
        main_layout = QVBoxLayout()  # create a vertical stacked layout
        self.work_area_tabs = WorkAreaTabs()
        self.main_btns = MainButtons()
        main_layout.addWidget(self.work_area_tabs)  # first part of vertical layout <= main work area
        main_layout.addWidget(self.main_btns)  # second part of vertical layout <= buttons
        self.main_btns.btn2.clicked.connect(self.prev_tab)
        self.main_btns.btn3.clicked.connect(self.next_tab)
        self.work_area_tabs.currentChanged.connect(self.tab_changed)
        self._third_btn_at_gen = False  # to reduce unnecessary ui changes

        # QMainWindow already has a fixed layout, so create a widget then
        # apply layout on it, and set it as the central widget
        main_widget = QWidget()
        main_widget.setLayout(main_layout)
        self.setCentralWidget(main_widget)

    @pyqtSlot()
    def tab_changed(self):
        # when last tab is on, change next button to generate button
        current_tab_idx = self.work_area_tabs.currentIndex()
        if current_tab_idx == (self.work_area_tabs.count()-1):  # if on last tab
            self.main_btns.btn3.clicked.disconnect()  # remove all triggers
            self.main_btns.btn3.clicked.connect(self.generate_script)  # can write script now
            self.main_btns.btn3.setText("Generate")
            self.main_btns.btn3.setIcon(qApp.style().standardIcon(QStyle.SP_DialogSaveButton))
            self._third_btn_at_gen = True
        else:  # not on last tab
            if not self._third_btn_at_gen:  # if button not changed, do nothing
                return None
            # otherwise restore original button
            self.main_btns.btn3.clicked.disconnect()
            self.main_btns.btn3.clicked.connect(self.next_tab)
            self.main_btns.btn3.setText("Next")
            self.main_btns.btn3.setIcon(qApp.style().standardIcon(QStyle.SP_ArrowRight))
            self._third_btn_at_gen = False

    @pyqtSlot()
    def next_tab(self):
        current_tab_idx = self.work_area_tabs.currentIndex()
        if current_tab_idx < (self.work_area_tabs.count()-1):
            self.work_area_tabs.setCurrentIndex(current_tab_idx+1)  # go to next tab

    @pyqtSlot()
    def prev_tab(self):
        current_tab_idx = self.work_area_tabs.currentIndex()
        if current_tab_idx > 0:
            self.work_area_tabs.setCurrentIndex(current_tab_idx-1)  # go to previous tab

    @pyqtSlot()
    def generate_script(self):
        write_ade_script_from_config(self)
