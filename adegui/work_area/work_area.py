from PyQt5.QtWidgets import (QTabWidget)
from adegui.work_area.job_setup_tab import JobSetupTab
from adegui.work_area.molecule_setup_tab import MoleculeSelectTab
from adegui.work_area.qm_method_setup_tab import QMMethodSetupTab


class WorkAreaTabs(QTabWidget):
    """ Sets up the working area, with multiple tabs """

    def __init__(self):
        super().__init__()
        # set up all the tabs one by one
        self.addTab(MoleculeSelectTab(), "Molecule details")
        self.addTab(JobSetupTab(), "Job Setup")
        self.addTab(QMMethodSetupTab(), "QM Method")
