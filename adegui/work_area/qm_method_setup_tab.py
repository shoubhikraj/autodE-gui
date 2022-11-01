from PyQt5.QtWidgets import (QWidget, QGroupBox, QComboBox,
                             QVBoxLayout, QHBoxLayout)
from PyQt5.QtCore import pyqtSlot
from adegui import Config


class QMMethodSetupTab(QWidget):
    """
    This part sets up the QM Method tab in the work area.
    QM Method tab sets up the underlying softwares that autodE uses
    (e.g. xTB, Gaussian, ORCA etc.)
    """
    def __init__(self):
        super().__init__()

        lmethod_setup = QGroupBox("Low-accuracy method")
        self.lmethod_drop_menu = QComboBox()
        for item in Config.ade_avail_lmethods:
            self.lmethod_drop_menu.addItem(item)
        Config.ade_lmethod = self.lmethod_drop_menu.currentText()  # initialize lmethod
        self.lmethod_drop_menu.currentIndexChanged.connect(self.lmethod_changed)

        lmethod_layout = QVBoxLayout()
        lmethod_layout.addWidget(self.lmethod_drop_menu)
        lmethod_layout.addStretch()
        lmethod_setup.setLayout(lmethod_layout)

        hmethod_setup = QGroupBox("High-accuracy method")
        self.hmethod_drop_menu = QComboBox()
        for item in Config.ade_avail_hmethods:
            self.hmethod_drop_menu.addItem(item)
        Config.ade_hmethod = self.hmethod_drop_menu.currentText()  # initialize hmethod

        hmethod_layout = QVBoxLayout()
        hmethod_layout.addWidget(self.hmethod_drop_menu)
        hmethod_layout.addStretch()
        hmethod_setup.setLayout(hmethod_layout)

        large_layout = QHBoxLayout() # lmethod on left, hmethod on right
        large_layout.addWidget(lmethod_setup, stretch=2)
        large_layout.addWidget(hmethod_setup, stretch=3)
        self.setLayout(large_layout)

    @pyqtSlot()
    def lmethod_changed(self):
        Config.ade_lmethod = self.lmethod_drop_menu.currentText()


    @pyqtSlot()
    def hmethod_change(self):
        current_hmethod = self.hmethod_drop_menu.currentText()
        if current_hmethod == 'Gaussian09': # fix the names of Gaussian
            Config.ade_hmethod = 'G09'
        elif current_hmethod == 'Gaussian16':
            Config.ade_hmethod = 'G16'
        else:
            Config.ade_hmethod = current_hmethod

