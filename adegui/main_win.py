from PyQt5.QtWidgets import QMainWindow, QVBoxLayout, QWidget
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
        self.main_btns.next_btn.clicked.connect(self.next_tab)
        self.main_btns.generate_btn.clicked.connect(self.generate_script)

        # QMainWindow already has a fixed layout, so create a widget then
        # apply layout on it, and set it as the central widget
        main_widget = QWidget()
        main_widget.setLayout(main_layout)
        self.setCentralWidget(main_widget)

    @pyqtSlot()
    def next_tab(self):
        current_tab_idx = self.work_area_tabs.currentIndex()
        if current_tab_idx < (self.work_area_tabs.count()-1):
            self.work_area_tabs.setCurrentIndex(current_tab_idx+1)  # go to next tab
        if (self.work_area_tabs.count()-1) == (current_tab_idx+1):  # if at last tab
            self.main_btns.btn3.setCurrentIndex(1)  # show generate button

    @pyqtSlot()
    def prev_tab(self):
        current_tab_idx = self.work_area_tabs.currentIndex()
        if current_tab_idx > 0:
            self.work_area_tabs.setCurrentIndex(current_tab_idx-1)  # go to previous tab
        if current_tab_idx <= (self.work_area_tabs.count()-1):  # if not at last tab
            self.main_btns.btn3.setCurrentIndex(0)

    @pyqtSlot()
    def generate_script(self):
        write_ade_script_from_config(self)
