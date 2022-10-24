import subprocess
import os
from typing import List, Tuple
from PyQt5.QtWidgets import (QTabWidget, QWidget, QGroupBox,
                             QVBoxLayout, QLineEdit, QStackedWidget,
                             QPushButton, QComboBox, QHBoxLayout,
                             QPlainTextEdit,QLabel)
from PyQt5.QtCore import pyqtSlot
import rdkit.Chem
from adegui.common import smiles_to_3d_rdkmol, _read_priv_rsrc_txt
from adegui import Config

# file operations in this script so need workdir
cwd = Config.adegui_workdir


class WorkAreaTabs(QTabWidget):
    """ Sets up the working area, with multiple tabs """

    def __init__(self):
        super().__init__()
        # set up all the tabs one by one
        self.addTab(MoleculeSelectTab(), "Molecule details")
        self.addTab(JobSetupTab(), "Job Setup")
        self.addTab(QMMethodSetupTab(), "QM Method")


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
        rct_layout.addWidget(MoleculeDrawOrType(Config.ade_rct_smis, 0, 'rct'))
        rct_layout.addWidget(MoleculeDrawOrType(Config.ade_rct_smis, 1, 'rct'))
        rct_layout.addStretch()  # to prevent unnecessary whitespace between widgets
        # two reactants right now can add more later
        rct_input.setLayout(rct_layout)

        # SMILES input or draw products
        prod_input = QGroupBox("Product(s)")
        prod_layout = QVBoxLayout()
        prod_layout.addWidget(MoleculeDrawOrType(Config.ade_prod_smis, 0, 'prod'))
        prod_layout.addWidget(MoleculeDrawOrType(Config.ade_prod_smis, 1, 'prod'))
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
    def __init__(self, smi_list: list, index: int, rct_or_prod: str):
        """
        Initialize a molecule draw/type widget
        :param smi_list: List of SMILES strings
        :param index: index to modify in list
        :param rct_or_prod: needed to get unique filenames on disk
        """
        super().__init__()
        self.smi_list = smi_list
        self.smi_index = index
        self.mol_fname = 'molecule-'+str(rct_or_prod)+str(index)+'.mol'

        self.smi_textbox = QLineEdit()
        self.smi_textbox.textEdited.connect(self.molecule_written)
        self.smi_textbox.textChanged.connect(self.molecule_changed)
        self.smi_list[self.smi_index] = self.smi_textbox.text()  # initialize
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
        os.remove(cwd+self.mol_fname) if os.path.isfile(cwd+self.mol_fname) else None
        return None

    @pyqtSlot()
    def molecule_drawn(self):
        """ Triggered when molecule is drawn """
        if not os.path.isfile(cwd+self.mol_fname):
            mol = smiles_to_3d_rdkmol(self.smi_textbox.text())
            # if there is no SMILES or illegal SMILES, it will go to default CH4
            if mol is None:
                mol = smiles_to_3d_rdkmol('C')
            rdkit.Chem.MolToMolFile(mol, cwd+self.mol_fname)
        # else carry on to display
        subprocess.run(['avogadro', cwd+self.mol_fname])
        mol = rdkit.Chem.MolFromMolFile(cwd+self.mol_fname)
        smi = rdkit.Chem.MolToSmiles(mol)
        self.smi_textbox.setText(smi)
        return None

    @pyqtSlot()
    def molecule_changed(self):
        """ Any change to molecule (typing or drawing) """
        self.smi_list[self.smi_index] = self.smi_textbox.text()

# TODO replace the second item with a widget
gui_avail_job_typs: List[Tuple[str, str]] = [
    ("Reaction Profile", "reaction_profile_info.txt"),
    ("Transition State: adaptive", "adaptive_ts.txt"),
    ("Transition State: CI-NEB", "ci_neb_ts.txt"),
    ("Reaction path: NEB", "")
]


class JobSetupTab(QWidget):
    """
    This part sets up the Job Setup tab in the work area.
    Job setup means the type of calculation for autodE to run.
    """

    def __init__(self):
        super().__init__()
        # setup the drop down menu
        self.drop_menu = QComboBox()
        for item in gui_avail_job_typs:
            self.drop_menu.addItem(item[0])
            # adds job type names i.e. Reaction Profile etc.
        Config.ade_job_typ = self.drop_menu.currentText()  # initialize
        self.drop_menu.currentIndexChanged.connect(self.jobtyp_changed)

        # put label and the dropdown menu below it
        job_setup_box = QGroupBox("autodE Job Type:")
        job_setup_layout = QVBoxLayout()
        job_setup_layout.addWidget(self.drop_menu)
        job_setup_box.setLayout(job_setup_layout)

        # now add other options and info-box
        self.opt_info_box = OptInfoBox()

        # finally, set the layout
        large_layout = QVBoxLayout()
        large_layout.addWidget(job_setup_box)
        large_layout.addWidget(self.opt_info_box)
        self.setLayout(large_layout)

    @pyqtSlot()
    def jobtyp_changed(self):
        Config.ade_job_typ = str(self.drop_menu.currentText())
        current_job_idx = self.drop_menu.currentIndex()
        self.opt_info_box.setCurrentIndex(current_job_idx)


class OptInfoBox(QStackedWidget):
    """
    Creates the Information box in Job Type tab,
    with additional options for each autodE Job Type
    """
    def __init__(self):
        super().__init__()
        for item in gui_avail_job_typs:
            self.addWidget(QPlainTextEdit())
        #rct_profile_box = QPlainTextEdit()
        #rct_profile_box.setPlainText(_read_priv_rsrc_txt('reaction_profile_info.txt'))
        #rct_profile_box.setReadOnly(True)
        #self.addWidget(rct_profile_box)



class QMMethodSetupTab(QWidget):
    """
    This part sets up the QM Method tab in the work area.
    QM Method tab sets up the underlying softwares that autodE uses
    (e.g. xTB, Gaussian, ORCA etc.)
    """
    def __init__(self):
        super().__init__()
        large_layout = QHBoxLayout() # lmethod on left, hmethod on right

        lmethod_setup = QGroupBox("Low-accuracy method")
        self.lmethod_drop_menu = QComboBox()
        for item in Config.ade_avail_lmethods:
            self.lmethod_drop_menu.addItem(item)

        lmethod_layout = QVBoxLayout()
        lmethod_layout.addWidget(self.lmethod_drop_menu)
        lmethod_layout.addStretch()
        lmethod_setup.setLayout(lmethod_layout)

        hmethod_setup = QGroupBox("High-accuracy method")
        self.hmethod_drop_menu = QComboBox()
        for item in Config.ade_avail_hmethods:
            self.hmethod_drop_menu.addItem(item)

        hmethod_layout = QVBoxLayout()
        hmethod_layout.addWidget(self.hmethod_drop_menu)
        hmethod_layout.addStretch()
        hmethod_setup.setLayout(hmethod_layout)

        large_layout.addWidget(lmethod_setup)
        large_layout.addWidget(hmethod_setup)
        self.setLayout(large_layout)


