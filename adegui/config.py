# Contains the config class which is instantiated once at the beginning
# of the program
import platform, subprocess
from typing import List, Union
import tempfile, pathlib, atexit, shutil
import errno


class AdeGuiMolecule:
    """
    Class for a molecule from the molecule setup tab
    """
    def __init__(self, molecule: Union[pathlib.Path, str], charge: int = 0, mult: int = 1):
        """
        Creates a new molecule instance. The first argument (molecule) has to be
        a path or a string. The path is always of a xyzfile, but if it's a string
        it should be a SMILES.

        :param molecule: The molecule path (xyzfile) or SMILES string (can be empty string)
        :param charge: charge of molecule, default 0
        :param mult: multiplicity of molecule, default 1
        """
        self.molecule = molecule  # can be path or SMILES string
        self.charge = charge
        self.mult = mult


class _AdeGuiConfig:
    def __init__(self):

        self.adegui_workdir: str = ''  # empty means current working dir for python
        temp_dir = tempfile.mkdtemp() # temp directory
        #self.adegui_scratchdir: pathlib.Path = pathlib.Path(temp_dir)
        self.adegui_scratchdir = pathlib.Path.cwd()

        # search for molecule editor
        self.adegui_moleditor = None
        for editor in ['avogadro', 'Avogadro']:  # all editor names possible
            if shutil.which(editor) is not None:
                self.adegui_moleditor = shutil.which(editor)
                break

        self.ade_n_cores: int = 8  # Number of cores in autodE
        self.ade_job_name: str = ''

        # initialize reactants and products
        self.ade_rct_mols: List[AdeGuiMolecule] = [AdeGuiMolecule(''), AdeGuiMolecule('')]  # Reactant(s)
        self.ade_prod_mols: List[AdeGuiMolecule] = [AdeGuiMolecule(''), AdeGuiMolecule('')]  # Product(s)

        self.ade_job_typ: str = ''  # Job Type

        self.ade_avail_lmethods: List[str] = ['XTB', 'MOPAC']
        self.ade_avail_hmethods: List[str] = ['ORCA', 'Gaussian09', 'Gaussian16', 'NWChem', 'QChem']
        self.ade_lmethod: str = ''  # chosen lmethod
        self.ade_hmethod: str = ''  # chosen hmethod

        #atexit.register(self._cleanup)  # TODO activate this once debugging done

    def _cleanup(self):
        try:
            shutil.rmtree(self.adegui_scratchdir)  # remove the temp directory
        except OSError as exc:
            if exc.errno != errno.ENOENT:
                raise
        return None
        # https://stackoverflow.com/questions/6884991/how-to-delete-a-directory-created-with-tempfile-mkdtemp


Config = _AdeGuiConfig()
