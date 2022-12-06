import os
from typing import List

import rdkit.Chem
from PyQt5.QtCore import pyqtSlot, QProcess
from PyQt5.QtWidgets import (QWidget, QHBoxLayout, QGroupBox,
                             QVBoxLayout, QLineEdit, QPushButton,
                             QLabel, QMessageBox, QSpinBox,
                             QFrame, qApp, QStyle)
from adegui import Config
from adegui.config import AdeGuiMolecule
from adegui.common import smiles_to_3d_rdkmol

# file operations in this script so need scratchdir
scrdir = Config.adegui_scratchdir


class MoleculeSelectTab(QWidget):
    """
    This area includes code for inputting SMILES or drawing
    a molecule for reactant and product
    """
    def __init__(self):
        super().__init__()
        # setup grid
        large_layout = QHBoxLayout()  # reactant on left, product on right

        # SMILES input or draw reactants
        rct_input = QGroupBox("Reactant(s)")
        rct_add_btn = QPushButton("+ Add reactant")
        rct_add_btn.clicked.connect(self.add_reactant_field)
        rct_lower_bar = QHBoxLayout()
        rct_lower_bar.addWidget(rct_add_btn)
        rct_lower_bar.addStretch()

        self.rct_dynamic_area = QVBoxLayout()
        self.rct_widgets = []  # holds references to widgets
        self.rct_max_len = 0

        rct_layout = QVBoxLayout()
        rct_layout.addLayout(self.rct_dynamic_area)
        rct_layout.addLayout(rct_lower_bar)
        rct_layout.addStretch()  # to prevent unnecessary whitespace between widgets
        # two reactants right now can add more later
        rct_input.setLayout(rct_layout)

        # SMILES input or draw products
        prod_input = QGroupBox("Product(s)")
        prod_add_btn = QPushButton("+ Add product")
        prod_add_btn.clicked.connect(self.add_product_field)
        prod_lower_bar = QHBoxLayout()
        prod_lower_bar.addWidget(prod_add_btn)
        prod_lower_bar.addStretch()

        prod_layout = QVBoxLayout()
        prod_layout.addWidget(MoleculeDrawOrType(Config.ade_prod_mols, 0, 'prod'))
        prod_layout.addLayout(prod_lower_bar)
        prod_layout.addStretch()
        prod_input.setLayout(prod_layout)

        large_layout.addWidget(rct_input)
        large_layout.addWidget(prod_input)
        self.setLayout(large_layout)

    @pyqtSlot()
    def add_reactant_field(self):
        Config.ade_rct_mols.append(AdeGuiMolecule(''))
        rct_widget = MoleculeDrawOrType(Config.ade_rct_mols, self.rct_max_len, 'rct')
        self.rct_widgets.append(rct_widget)
        self.rct_dynamic_area.addWidget(rct_widget)
        self.rct_max_len += 1

    @pyqtSlot()
    def add_product_field(self):
        pass


class MoleculeDrawOrType(QFrame):
    """
    Each widget for typing or drawing molecule (drawn molecules are converted to SMILES)
    Initialize by passing a list of SMILES strings, and integer denoting position of the
    list that is modified by widget
    """
    def __init__(self, mols_list: List[AdeGuiMolecule], index: int, rct_or_prod: str):
        """
        Initialize a molecule draw/type widget

        :param mols_list: List of molecules, each molecule has filename/smiles, charge, mult
        :param index: index to modify in list
        :param rct_or_prod: needed to get unique filenames on disk
        """
        super().__init__()
        self.mol_list = mols_list
        self.mol_index = index
        self.mol_fname = 'molecule-'+str(rct_or_prod)+str(index)+'.mol'
        self.xyz_fname = self.mol_fname[:-4] + '.xyz'
        self.editor_proc = None  # molecular editor process

        self.smi_textbox = QLineEdit()
        self.smi_textbox.setPlaceholderText("Type SMILES or draw")
        self.smi_textbox.textEdited.connect(self.molecule_written)
        self.draw_btn = QPushButton("Draw")
        self.draw_btn.setIcon(qApp.style().standardIcon(QStyle.SP_DialogResetButton))
        self.draw_btn.clicked.connect(self.molecule_draw_start)

        upper_bar_layout = QHBoxLayout()
        upper_bar_layout.addWidget(QLabel("SMILES:"), stretch=0)
        upper_bar_layout.addWidget(self.smi_textbox, stretch=2)
        upper_bar_layout.addWidget(self.draw_btn, stretch=1)
        upper_bar = QWidget()
        upper_bar.setLayout(upper_bar_layout)

        self.charge_dial = QSpinBox()  # set the charge
        self.charge_dial.setRange(-20, 20)
        self.charge_dial.valueChanged.connect(self.charge_changed)
        self.mult_dial = QSpinBox()  # set the multiplicity
        self.mult_dial.setRange(1, 20)
        self.mult_dial.valueChanged.connect(self.mult_changed)

        lower_bar_layout = QHBoxLayout()
        lower_bar_layout.addWidget(QLabel("Charge:"), stretch=0)
        lower_bar_layout.addWidget(self.charge_dial, stretch=1)
        lower_bar_layout.addStretch(1)
        lower_bar_layout.addWidget(QLabel("Multiplicity:"), stretch=0)
        lower_bar_layout.addWidget(self.mult_dial, stretch=1)
        lower_bar = QWidget()
        lower_bar.setLayout(lower_bar_layout)

        large_layout = QVBoxLayout()
        large_layout.addWidget(upper_bar)
        large_layout.addWidget(lower_bar)
        large_layout.setContentsMargins(0, 0, 0, 0)  # remove padding

        self.setFrameStyle(QFrame.StyledPanel | QFrame.Plain)
        self.setLayout(large_layout)

    @pyqtSlot()
    def molecule_written(self):
        """ Triggered when molecule is typed in textbox """
        self.mol_list[self.mol_index].molecule = self.smi_textbox.text()  # assign SMILES
        # warning! This does not check if SMILES is sane
        if os.path.isfile(scrdir/self.mol_fname): os.remove(scrdir/self.mol_fname)  # remove old molfile
        return None

    @pyqtSlot()
    def molecule_draw_start(self):
        """ Triggered when molecule draw button is pressed """
        # if there is no editor then warn
        if Config.adegui_moleditor is None:
            QMessageBox.warning(self,
                                "autodE-GUI",
                                "No molecule drawing software available!")
            return None

        if not os.path.isfile(scrdir/self.mol_fname):  # if previous file does not exist make new mol
            mol = smiles_to_3d_rdkmol(self.smi_textbox.text())
            # if there is no SMILES or illegal SMILES, it will go to default CH4
            if mol is None:
                self.smi_textbox.setText('')  # empty textbox
                mol = smiles_to_3d_rdkmol('C')  # use CH4 as starting point
            rdkit.Chem.MolToMolFile(mol, str(scrdir/self.mol_fname))

        if self.editor_proc is None:  # protect against multiple launch
            self.draw_btn.setEnabled(False)  # disable button
            self.draw_btn.setToolTip("Avogadro is still running!")  # explain why it is disabled
            self.editor_proc = QProcess()  # async process
            self.editor_proc.finished.connect(self.molecule_draw_finished)  # when finished
            # use editor to edit the molfile
            self.editor_proc.start(Config.adegui_moleditor, [str(scrdir/self.mol_fname)])

    @pyqtSlot()
    def molecule_draw_finished(self):
        """ Triggered when molecule editor exits, i.e. drawing is finished """
        self.draw_btn.setEnabled(True)  # enable button again
        self.draw_btn.setToolTip('')  # clear tooltip
        self.editor_proc = None
        # now handle the molfile
        if os.path.isfile(scrdir/self.mol_fname):
            mol = rdkit.Chem.MolFromMolFile(str(scrdir/self.mol_fname), sanitize=False, removeHs=False)
        else:
            mol = None
        # if molfile is deleted or corrupted somehow
        if mol is None:
            QMessageBox.warning(self, "autodE-GUI", "Molfile unavailable or corrupted!")
            return None
        # try to display it as sanitized SMILES
        mol_copy = rdkit.Chem.Mol(mol)  # get a copy
        istat = rdkit.Chem.SanitizeMol(mol_copy, catchErrors=True)
        if istat == 0:  # if sanitize possible
            mol_copy = rdkit.Chem.RemoveHs(mol_copy)
            smi = rdkit.Chem.MolToSmiles(mol_copy)
            self.smi_textbox.setText(smi)
        else:  # just display unsanitized version
            smi = rdkit.Chem.MolToSmiles(mol)
            self.smi_textbox.setText(smi)
        # then convert to xyz and save its path
        rdkit.Chem.MolToXYZFile(mol, self.xyz_fname)
        self.mol_list[self.mol_index].molecule = scrdir/self.xyz_fname
        return None

    @pyqtSlot()
    def charge_changed(self):
        charge = self.charge_dial.value()
        self.mol_list[self.mol_index].charge = charge

    @pyqtSlot()
    def mult_changed(self):
        mult = self.mult_dial.value()
        if mult >= 1:
            self.mol_list[self.mol_index].mult = mult
        else:
            self.mult_dial.setValue(1)

