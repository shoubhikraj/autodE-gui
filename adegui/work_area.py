import subprocess
from PyQt5.QtWidgets import (QTabWidget, QWidget, QGridLayout,
                             QVBoxLayout, QLineEdit, QLabel,
                             QPushButton, QComboBox)
from PyQt5.QtCore import pyqtSlot
import rdkit.Chem
from adegui.common import smiles_to_3d_rdkmol


class WorkAreaTabs(QTabWidget):
    """ Sets up the working area, with multiple tabs """

    def __init__(self):
        super().__init__()
        # set up all the tabs one by one
        self.addTab(MoleculeSelectTab(), "Molecule details")
        self.addTab(JobSetupTab(), "Job Setup")
        self.addTab(QWidget(), "QM Method")


class MoleculeSelectTab(QWidget):
    """
    This area includes code for inputting SMILES or drawing
    a molecule for reactant and product
    """

    def __init__(self):
        super().__init__()
        # setup grid
        layout = QGridLayout()  # 2 rows, 6 columns
        layout.setColumnStretch(0, 1)
        layout.setColumnStretch(1, 4)  # =>
        layout.setColumnStretch(2, 1)  #
        layout.setColumnStretch(3, 1)  #
        layout.setColumnStretch(4, 4)  # => these are the SMILES textbox columns
        layout.setColumnStretch(5, 1)

        # SMILES input or draw reactants
        self.smi_textbox_rct = QLineEdit()  # for reactant SMILES input
        draw_ext_rct = QPushButton("Draw")  # to draw reactant using an external sketcher
        draw_ext_rct.clicked.connect(self.reactant_drawn)
        # SMILES input or draw products
        self.smi_textbox_prod = QLineEdit()
        draw_ext_prod = QPushButton("Draw")
        draw_ext_prod.clicked.connect(self.product_drawn)

        layout.addWidget(QLabel("Reactants"), 0, 0)
        layout.addWidget(QLabel("SMILES: "), 1, 0)
        layout.addWidget(self.smi_textbox_rct, 1, 1)
        layout.addWidget(draw_ext_rct, 1, 2)

        layout.addWidget(QLabel("Products"), 0, 3)
        layout.addWidget(QLabel("SMILES: "), 1, 3)
        layout.addWidget(self.smi_textbox_prod, 1, 4)
        layout.addWidget(draw_ext_prod, 1, 5)
        self.setLayout(layout)

    @pyqtSlot()
    def reactant_written(self):
        pass

    @pyqtSlot()
    def reactant_drawn(self):
        # TODO : take SMILES from the textbox and then convert to mol
        rct_mol = smiles_to_3d_rdkmol('C')
        rdkit.Chem.MolToMolFile(rct_mol, 'rct.mol')
        subprocess.run(['avogadro', 'rct.mol'])
        rct_mol = rdkit.Chem.MolFromMolFile('rct.mol')
        rct_smi = rdkit.Chem.MolToSmiles(rct_mol)
        self.smi_textbox_rct.setText(rct_smi)

    @pyqtSlot()
    def product_written(self):
        pass

    @pyqtSlot()
    def product_drawn(self):
        pass


class JobSetupTab(QWidget):
    """
    This part sets up the Job Setup tab in the work area.
    Job setup means the type of calculation for autodE to run.
    """

    def __init__(self):
        super().__init__()
        # setup the drop down menu
        drop_menu = QComboBox()
        drop_menu.addItem("Reaction Profile")
        drop_menu.addItem("Transition State: adaptive")
        drop_menu.addItem("Transition State: CI-NEB")
        drop_menu.addItem("Reaction Path: NEB")
        drop_menu.addItem("Reaction Path: CI-NEB")
        drop_menu.addItem("1D PES")
        drop_menu.addItem("2D PES")

        # put label and the dropdown menu below it
        layout = QVBoxLayout()
        layout.addWidget(QLabel("autodE Job Type:"))
        layout.addWidget(drop_menu)
        self.setLayout(layout)
