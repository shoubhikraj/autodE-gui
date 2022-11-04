import os
import subprocess
from typing import List

import rdkit.Chem
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import (QWidget, QHBoxLayout, QGroupBox,
                             QVBoxLayout, QLineEdit, QPushButton,
                             QLabel)
from adegui import Config
from adegui.common import smiles_to_3d_rdkmol

# file operations in this script so need workdir
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
        rct_layout = QVBoxLayout()
        rct_layout.addWidget(MoleculeDrawOrType(Config.ade_rct_mols, 0, 'rct'))
        rct_layout.addWidget(MoleculeDrawOrType(Config.ade_rct_mols, 1, 'rct'))
        rct_layout.addStretch()  # to prevent unnecessary whitespace between widgets
        # two reactants right now can add more later
        rct_input.setLayout(rct_layout)

        # SMILES input or draw products
        prod_input = QGroupBox("Product(s)")
        prod_layout = QVBoxLayout()
        prod_layout.addWidget(MoleculeDrawOrType(Config.ade_prod_mols, 0, 'prod'))
        prod_layout.addWidget(MoleculeDrawOrType(Config.ade_prod_mols, 1, 'prod'))
        prod_layout.addStretch()
        prod_input.setLayout(prod_layout)

        large_layout.addWidget(rct_input)
        large_layout.addWidget(prod_input)
        self.setLayout(large_layout)


class MoleculeDrawOrType(QWidget):
    """
    Each widget for typing or drawing molecule (drawn molecules are converted to SMILES)
    Initialize by passing a list of SMILES strings, and integer denoting position of the
    list that is modified by widget
    """
    def __init__(self, mols_list: List[dict], index: int, rct_or_prod: str):
        """
        Initialize a molecule draw/type widget
        :param mols_list: List of molecules, each molecule is a dict of filename, charge, mult
        :param index: index to modify in list
        :param rct_or_prod: needed to get unique filenames on disk
        """
        super().__init__()
        self.mol_list = mols_list
        self.smi_index = index
        self.mol_fname = 'molecule-'+str(rct_or_prod)+str(index)+'.mol'

        self.smi_textbox = QLineEdit()
        self.smi_textbox.textEdited.connect(self.molecule_written)
        self.smi_textbox.textChanged.connect(self.molecule_changed)
        self.mol_list[self.smi_index] = self.smi_textbox.text()  # initialize
        draw_btn = QPushButton("Draw")
        draw_btn.clicked.connect(self.molecule_drawn)

        layout = QHBoxLayout()
        layout.addWidget(QLabel("SMILES:"), 0)
        layout.addWidget(self.smi_textbox, 2)
        layout.addWidget(draw_btn, 1)
        self.setLayout(layout)

    @pyqtSlot()
    def molecule_written(self):
        """ Triggered when molecule is typed in textbox """
        os.remove(scrdir/self.mol_fname) if os.path.isfile(scrdir/self.mol_fname) else None
        return None

    @pyqtSlot()
    def molecule_drawn(self):
        """ Triggered when molecule is drawn """
        if not os.path.isfile(scrdir/self.mol_fname):
            mol = smiles_to_3d_rdkmol(self.smi_textbox.text())
            # if there is no SMILES or illegal SMILES, it will go to default CH4
            if mol is None:
                mol = smiles_to_3d_rdkmol('C')
            rdkit.Chem.MolToMolFile(mol, str(scrdir/self.mol_fname))
        # else carry on to display
        subprocess.run([Config.adegui_moleditor, str(scrdir/self.mol_fname)])
        mol = rdkit.Chem.MolFromMolFile(str(scrdir/self.mol_fname))
        smi = rdkit.Chem.MolToSmiles(mol)
        self.smi_textbox.setText(smi)
        return None

    @pyqtSlot()
    def molecule_changed(self):
        """ Any change to molecule (typing or drawing) """
        self.mol_list[self.smi_index] = self.smi_textbox.text()
