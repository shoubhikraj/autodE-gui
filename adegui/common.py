# This file contains common functions, constants etc. which area used
# in the main code
import importlib.resources
from typing import Type, Optional
from rdkit import Chem
from rdkit.Chem import AllChem


def smiles_to_3d_rdkmol(smi: str) -> Optional[Type[Chem.rdchem.Mol]]:
    """
    Takes a SMILES as a string and returns an RDKit mol object
    with a generated 3D geometry (with hydrogens added)
    :param smi: SMILES string
    :return: RDKit mol object
    """
    # TODO should I try to have a timeout in case EmbedMolecule hangs??
    try:
        smi = str(smi)
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)
        istat = AllChem.EmbedMolecule(mol, maxAttempts=5)
        if istat == 0:
            return mol
        else:
            return None
    except:
        return None

def _read_priv_rsrc_txt(fname: str) -> str:
    """
    Take a resource file name in adegui
    :param fname: Name of text file
    :return: String containing the contents of text file
    """
    with importlib.resources.path('adegui.resources',str(fname)) as fh:
        text = fh.read_text()
    return text

