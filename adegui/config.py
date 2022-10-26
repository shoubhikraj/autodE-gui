# Contains the config class which is instantiated once at the beginning
# of the program
from typing import List

class _AdeGuiConfig:
    def __init__(self):
        self.adegui_workdir: str = ''  # empty means current working dir for python

        self.adegui_moleditor: str = '/Applications/Avogadro.app/Contents/MacOS/Avogadro'  # path of the molecule editor

        self.ade_n_cores: int = 8  # Number of cores in autodE
        self.ade_job_name: str = ''

        self.ade_rct_smis: list = ['', '']  # Reactant(s)
        self.ade_prod_smis: list = ['', '']  # Product(s)

        self.ade_job_typ: str = ''  # Job Type

        self.ade_avail_lmethods: List[str] = ['XTB', 'MOPAC']
        self.ade_avail_hmethods: List[str] = ['ORCA', 'Gaussian09', 'Gaussian16', 'NWChem', 'QChem']
        self.ade_lmethod: str = ''  # chosen lmethod
        self.ade_hmethod: str = ''  # chosen hmethod



Config = _AdeGuiConfig()
