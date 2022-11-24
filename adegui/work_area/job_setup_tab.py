from typing import List
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import (QWidget, QComboBox, QGroupBox,
                             QHBoxLayout, QVBoxLayout, QStackedWidget,
                             QPlainTextEdit)
from adegui import Config
from adegui.common import _read_priv_rsrc_txt

gui_avail_job_typs: List[str] = [
    "Reaction Profile",
    "Transition State",
    "Reaction path: CI-NEB",
    "Reaction path: NEB"
]  # we need the list here to make the order right


class JobSetupTab(QWidget):
    """
    This part sets up the Job Setup tab in the work area.
    Job setup means the type of calculation for autodE to run.
    """

    def __init__(self):
        super().__init__()
        # setup the drop-down menu
        self.drop_menu = QComboBox()
        for item in gui_avail_job_typs:
            self.drop_menu.addItem(item)
            # adds job type names i.e. Reaction Profile etc.
        Config.ade_job_typ = self.drop_menu.currentText()  # initialize
        self.drop_menu.currentIndexChanged.connect(self.jobtyp_changed)

        # put label and the dropdown menu below it
        job_setup_box = QGroupBox("autodE Job Type:")
        job_setup_layout = QHBoxLayout()
        job_setup_layout.addWidget(self.drop_menu)
        job_setup_layout.addStretch()
        job_setup_box.setLayout(job_setup_layout)

        # now add other Options and info-box
        self.opt_info_box = OptInfoBox()

        # finally, set the larger layout
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
    # Different job types need different widgets, so it has to be hand coded
    def __init__(self):
        super().__init__()
        rct_profile_box = QPlainTextEdit()
        rct_profile_box.setPlainText(_read_priv_rsrc_txt('reaction_profile_info.txt'))
        rct_profile_box.setReadOnly(True)
        self.addWidget(rct_profile_box)
        # TODO fix this part
        ts_box = QPlainTextEdit()
        ts_box.setPlainText(_read_priv_rsrc_txt('ts_info.txt'))
        ts_box.setReadOnly(True)
        self.addWidget(ts_box)
        ci_neb_box = QPlainTextEdit()
        ci_neb_box.setPlainText(_read_priv_rsrc_txt('ci_neb_info.txt'))
        ci_neb_box.setReadOnly(True)
        self.addWidget(ci_neb_box)
        neb_box = QPlainTextEdit()
        neb_box.setPlainText(_read_priv_rsrc_txt('neb_info.txt'))
        neb_box.setReadOnly(True)
        self.addWidget(neb_box)

