# This file contains common functions, constants etc. which are used
# in the main code
import shutil

import importlib.resources
import pathlib, os
from typing import Type, Optional, Union
from rdkit import Chem
from rdkit.Chem import AllChem
from PyQt5.QtWidgets import QMessageBox


def smiles_to_3d_rdkmol(smi: str) -> Optional[Type[Chem.rdchem.Mol]]:
    """
    Takes a SMILES as a string and returns an RDKit mol object
    with a generated 3D geometry (with hydrogens added)

    :param smi: SMILES string
    :return: RDKit mol object
    """
    # TODO should I try to have a timeout in case EmbedMolecule hangs??
    if str(smi) == '':
        return None
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
    Read a resource file in adegui

    :param fname: Name of text file
    :return: String containing the contents of text file
    """
    with importlib.resources.path('adegui.resources', str(fname)) as fh:
        text = fh.read_text()
    return text


def _safe_copy_file(source: pathlib.Path, dest: pathlib.Path, obj) -> int:
    """
    Copy a file from source to destination with multiple checks for
    extra safety

    :param source: The source file (pathlib Path)
    :param dest: The destination file (pathlib Path)
    :param obj: The parent QWidget instance
    :return: 0 if successfull, -1 if unsuccessful
    """
    if not os.path.isfile(source):
        QMessageBox.critical(obj,
                             "autodE-GUI",
                             f"Error, file {source.name} has been deleted or file system "
                             "permission error! Unable to write file.")
        return -1
    else:
        try:
            shutil.copyfile(source, dest)
        except shutil.SameFileError:  # if same file, no need to copy
            return 0
        return 0

