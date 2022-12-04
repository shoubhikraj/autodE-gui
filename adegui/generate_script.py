import pathlib, os, shutil

from PyQt5.QtWidgets import QMessageBox, QFileDialog
from adegui import Config
from adegui.common import _safe_copy_file
from adegui.work_area.job_setup_tab import gui_avail_job_typs


def write_ade_script_from_config(obj) -> None:
    """
    Writes an autodE run script (Python) based on current config vars (adegui.Config)

    :param obj: The parent QWidget instance (required for the warning windows)
    :return: None
    """

    # stop if there are no reactants or products
    if (Config.ade_rct_mols[0].molecule, Config.ade_rct_mols[1].molecule) == ('', ''):
        QMessageBox.warning(obj,
                            "autodE-GUI",
                            "There are no reactants! Unable to write script.")
        return None
    if (Config.ade_prod_mols[0].molecule, Config.ade_prod_mols[1].molecule) == ('', ''):
        QMessageBox.warning(obj,
                            "autodE-GUI",
                            "There are no products! Unable to write script.")
        return None

    save_fname, _ = QFileDialog.getSaveFileName(obj,
                                                caption="Save autodE script",
                                                filter="Python Files (*.py)")
    if save_fname == '':
        return None
    else:
        save_fpath = pathlib.Path(save_fname)
    cwd_path = save_fpath.parent  # get path for the working directory (from save file dialog)


    gen_script = []  # list of lines of text for the script to be generated

    # first import and set the autodE config variables
    gen_script.append("import autode as ade\n\n")
    gen_script.append(f"ade.Config.n_cores = {Config.ade_n_cores}\n")
    gen_script.append(f"ade.Config.max_core = {Config.ade_max_core_mem:.2f}\n")
    gen_script.append("\n")  # line breaks to make it look better

    # setup the lmethod and hmethod
    gen_script.append(f"ade.Config.lcode = '{Config.ade_lmethod}'\n")
    gen_script.append(f"ade.Config.hcode = '{Config.ade_hmethod}'\n")

    # options for hmethod (basis and functional)
    # TODO what is low_sp and how should it be handled?
    if Config.ade_hmethod_sp_basis != '':
        gen_script.append(f"ade.Config.{Config.ade_hmethod}.keywords."
                          f"sp.basis_set = '{Config.ade_hmethod_sp_basis}'\n")
        gen_script.append(f"ade.Config.{Config.ade_hmethod}.keywords."
                          f"low_sp.basis_set = '{Config.ade_hmethod_sp_basis}'\n")
    if Config.ade_hmethod_sp_func != '':
        gen_script.append(f"ade.Config.{Config.ade_hmethod}.keywords."
                          f"sp.functional = '{Config.ade_hmethod_sp_func}'\n")
        gen_script.append(f"ade.Config.{Config.ade_hmethod}.keywords."
                          f"low_sp.basis_set = '{Config.ade_hmethod_sp_basis}'\n")
    if Config.ade_hmethod_geom_basis != '':
        gen_script.append(f"ade.Config.{Config.ade_hmethod}.keywords."
                          f"set_opt_basis('{Config.ade_hmethod_geom_basis}')\n")
    if Config.ade_hmethod_geom_func != '':
        gen_script.append(f"ade.Config.{Config.ade_hmethod}.keywords."
                          f"set_opt_functional('{Config.ade_hmethod_geom_func}')\n")
    gen_script.append("\n")

    # then set reactants and products
    rct_and_prod = ''  # need for the final reaction setup line

    for index, rct_molecule in enumerate(Config.ade_rct_mols):
        if not rct_molecule.molecule == '':
            if isinstance(rct_molecule.molecule, pathlib.Path):
                xyz_fname = rct_molecule.molecule.name
                gen_script.append(f"rct{index} = ade.Reactant('{xyz_fname}', ")
                if _safe_copy_file(rct_molecule.molecule, cwd_path/xyz_fname, obj) != 0:
                    return None
            elif isinstance(rct_molecule.molecule, str):
                gen_script.append(f"rct{index} = ade.Reactant(smiles='{rct_molecule.molecule}', ")
            gen_script.append(f"charge={rct_molecule.charge}, "
                              f"mult={rct_molecule.mult})\n")
            rct_and_prod += f"rct{index},"

    for index, prod_molecule in enumerate(Config.ade_prod_mols):
        if not prod_molecule.molecule == '':
            if isinstance(prod_molecule.molecule, pathlib.Path):
                xyz_fname = prod_molecule.molecule.name
                gen_script.append(f"prod{index} = ade.Product('{xyz_fname}', ")
                if _safe_copy_file(prod_molecule.molecule, cwd_path/xyz_fname, obj) != 0:
                    return None
            elif isinstance(prod_molecule.molecule, str):
                gen_script.append(f"prod{index} = ade.Product(smiles='{prod_molecule.molecule}', ")
            gen_script.append(f"charge={prod_molecule.charge}, "
                              f"mult={prod_molecule.mult})\n")
            rct_and_prod += f"prod{index},"
    gen_script.append("\n")

    # calculation setup
    if Config.ade_job_typ not in gui_avail_job_typs:
        raise Exception("Job type is not available!")

    # Reaction profile =>
    if Config.ade_job_typ == 'Reaction Profile':
        gen_script.append(f"rxn = ade.Reaction({rct_and_prod})\n")
        gen_script.append(f"rxn.calculate_reaction_profile()")
    # <=

    # Transition state (adaptive) =>
    if Config.ade_job_typ == "Transition State":
        gen_script.append(f"rxn = ade.Reaction({rct_and_prod})\n")
        gen_script.append(f"rxn.locate_transition_state()\n")
        gen_script.append("rxn.ts.print_xyz_file(filename='ts.xyz')\n")
    # <=

    # (over)write script
    try:
        with open(save_fpath, 'w') as fh:
            fh.writelines(gen_script)
        # if success
        QMessageBox.information(obj, "autodE-GUI", "Finished generating script file")
        return None
    except OSError:
        QMessageBox.critical(obj,
                             "autodE-GUI",
                             "OSError: unable to writing file (file system error)")
        return None
