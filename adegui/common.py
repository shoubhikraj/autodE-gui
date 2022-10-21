# This file contains common functions, constants etc. which is used
# repeatedly in the main code
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
