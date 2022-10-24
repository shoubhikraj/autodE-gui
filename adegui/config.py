# Contains the config class which is instantiated once at the beginning
# of the program

class _AdeGuiConfig:
    def __init__(self):
        self.adegui_workdir: str = ''  # empty means current working dir for python

        self.ade_n_cores: int = 8  # Number of cores in autodE
        self.ade_job_name: str = ''

        self.ade_rct_smis: list = ['', '']  # Reactant(s)
        self.ade_prod_smis: list = ['', '']  # Product(s)

        self.ade_job_typ: str = None  # Job Type


        self.ade_avail_lmethods: list[str] = ['XTB', 'MOPAC']
        self.ade_avail_hmethods: list[str] = ['ORCA', 'Gaussian09', 'Gaussian16', 'NWChem', 'QChem']



Config = _AdeGuiConfig()
