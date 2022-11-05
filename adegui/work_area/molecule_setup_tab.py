import os
import subprocess
from typing import List

import rdkit.Chem
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import (QWidget, QHBoxLayout, QGroupBox,
                             QVBoxLayout, QLineEdit, QPushButton,
                             QLabel, QMessageBox)
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

        self.smi_textbox = QLineEdit()
        self.smi_textbox.textEdited.connect(self.molecule_written)
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
        self.mol_list[self.mol_index].molecule = self.smi_textbox.text()  # assign SMILES
        # warning! This does not check if SMILES is sane
        return None

    @pyqtSlot()
    def molecule_drawn(self):
        """ Triggered when molecule is drawn """
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
                mol = smiles_to_3d_rdkmol('C')
            rdkit.Chem.MolToMolFile(mol, str(scrdir/self.mol_fname))
        # use editor to edit the molfile
        subprocess.run([Config.adegui_moleditor, str(scrdir/self.mol_fname)])
        mol = rdkit.Chem.MolFromMolFile(str(scrdir/self.mol_fname), sanitize=False, removeHs=False)
        # try to display it as sanitized SMILES
        mol_copy = rdkit.Chem.Mol(mol)  # get a copy
        istat = rdkit.Chem.SanitizeMol(mol_copy, catchErrors=True)
        if istat == 0: # if sanitize possible
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
