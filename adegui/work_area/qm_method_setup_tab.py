from PyQt5.QtWidgets import (QWidget, QGroupBox, QComboBox,
                             QVBoxLayout, QHBoxLayout, QSizePolicy,
                             QLabel)
from PyQt5.QtCore import pyqtSlot
from adegui import Config


avail_lmethods = ['XTB', 'MOPAC']
avail_hmethods = ['ORCA', 'Gaussian09', 'Gaussian16', 'NWChem', 'QChem']

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
        for item in avail_lmethods:
            self.lmethod_drop_menu.addItem(item)
        Config.ade_lmethod = self.lmethod_drop_menu.currentText()  # initialize lmethod
        self.lmethod_drop_menu.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)  # no horizontal stretch
        self.lmethod_drop_menu.currentIndexChanged.connect(self.lmethod_changed)  # signal connect

        lmethod_layout = QVBoxLayout()
        lmethod_layout.addWidget(self.lmethod_drop_menu)
        lmethod_layout.addStretch()
        lmethod_setup.setLayout(lmethod_layout)

        hmethod_setup = QGroupBox("High-accuracy method")
        self.hmethod_drop_menu = QComboBox()
        for item in avail_hmethods:
            self.hmethod_drop_menu.addItem(item)
        Config.ade_hmethod = self.hmethod_drop_menu.currentText()  # initialize hmethod
        self.hmethod_drop_menu.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        self.hmethod_drop_menu.currentIndexChanged.connect(self.hmethod_changed)

        qm_hmethod_mod_widget = QWidget()  # modify hmethod (basis and functionals)
        qm_mod_layout = QHBoxLayout()
        qm_mod_layout.addWidget(QMModifyWidget('Single-point', 'sp'))
        qm_mod_layout.addWidget(QMModifyWidget('Geometries (gradient+hessian)', 'geom'))
        qm_hmethod_mod_widget.setLayout(qm_mod_layout)

        hmethod_layout = QVBoxLayout()
        hmethod_layout.addWidget(self.hmethod_drop_menu)
        hmethod_layout.addWidget(qm_hmethod_mod_widget)
        hmethod_layout.addStretch()
        hmethod_setup.setLayout(hmethod_layout)

        large_layout = QHBoxLayout()  # lmethod on left, hmethod on right
        large_layout.addWidget(lmethod_setup, stretch=2)
        large_layout.addWidget(hmethod_setup, stretch=3)
        self.setLayout(large_layout)

    @pyqtSlot()
    def lmethod_changed(self):
        Config.ade_lmethod = self.lmethod_drop_menu.currentText()

    @pyqtSlot()
    def hmethod_changed(self):
        current_hmethod = self.hmethod_drop_menu.currentText()
        if current_hmethod == 'Gaussian09': # fix the names of Gaussian
            Config.ade_hmethod = 'G09'
        elif current_hmethod == 'Gaussian16':
            Config.ade_hmethod = 'G16'
        else:
            Config.ade_hmethod = current_hmethod


avail_basis_sets = ['(default)', 'def2-SVP', 'def2-TZVP', 'ma-def2-SVP',
                    'ma-def2-TZVP', 'cc-pVDZ', 'aug-cc-pVDZ', 'cc-pVTZ',
                    'aug-cc-pVTZ']
avail_functionals = ['(default)', 'B3LYP', 'M06-2X', 'X3LYP',
                     'PBE0', 'TPSS']


class QMModifyWidget(QGroupBox):
    """
    This widget modifies the basis set and the functional for the current QM Method,
    either for the single point calculations or geometry calculations
    """
    def __init__(self, name: str, sp_or_geom: str):
        """
        Creates a new QM Method Modify widget instance

        :param name: Heading of widget
        :param sp_or_geom: whether to modify single point(SP) or geometry config
        """
        super().__init__(name)
        self.sp_or_geom_mod = sp_or_geom  # which to modify -> used later

        self.basis_menu = QComboBox()
        for basis_set in avail_basis_sets:
            self.basis_menu.addItem(basis_set)
        self.basis_menu.currentIndexChanged.connect(self.basis_changed)
        basis_layout = QHBoxLayout()
        basis_layout.addWidget(QLabel('Basis Set:'))
        basis_layout.addWidget(self.basis_menu)
        basis_widget = QWidget()  # use a dummy widget to push the layout onto
        basis_widget.setLayout(basis_layout)

        self.func_menu = QComboBox()
        for functional in avail_functionals:
            self.func_menu.addItem(functional)
        self.func_menu.currentIndexChanged.connect(self.func_changed)
        func_layout = QHBoxLayout()
        func_layout.addWidget(QLabel('Functional:'))
        func_layout.addWidget(self.func_menu)
        func_widget = QWidget()
        func_widget.setLayout(func_layout)

        large_layout = QVBoxLayout()
        large_layout.addWidget(basis_widget)
        large_layout.addWidget(func_widget)
        self.setLayout(large_layout)

    @pyqtSlot()
    def basis_changed(self):
        basis_selected = self.basis_menu.currentText()
        if basis_selected == '(default)': basis_selected = ''
        if self.sp_or_geom_mod == 'sp':
            Config.ade_hmethod_sp_basis = basis_selected
        elif self.sp_or_geom_mod == 'geom':
            Config.ade_hmethod_geom_basis = basis_selected
        else:
            raise Exception('unknown option')

    @pyqtSlot()
    def func_changed(self):
        func_selected = self.func_menu.currentText()
        if func_selected == '(default)': func_selected = ''
        if self.sp_or_geom_mod == 'sp':
            Config.ade_hmethod_sp_func = func_selected
        elif self.sp_or_geom_mod == 'geom':
            Config.ade_hmethod_geom_func = func_selected
        else:
            raise Exception('unknown option')
