from PyQt5.QtWidgets import QMessageBox
from adegui import Config

cwd = Config.adegui_workdir


def write_ade_script_from_config(obj) -> None:
    """
    Writes an autodE run script (Python) based on current config vars (adegui.Config)
    Name of generated script is "aderun.py"
    :param obj: The parent QWidget instance (required for the warning windows)
    :return: None
    """
    # TODO: make the name of script editable (??)
    fname = 'aderun.py'  # name of file <= should this be changeable?
    gen_script = []  # list of lines of text for the script to be generated

    # first import and set the autodE config variables
    gen_script.append("import autode as ade\n\n")
    gen_script.append(f"ade.Config.n_cores = {Config.ade_n_cores}\n")
    gen_script.append("\n")

    # setup the lmethod and hmethod
    gen_script.append(f"ade.Config.lcode = '{Config.ade_lmethod}'\n")
    gen_script.append(f"ade.Config.hcode = '{Config.ade_hmethod}'\n")
    gen_script.append("\n")  # line break to make it look better

    # then set reactants and products
    rct_and_prod = ''  # need for the final reaction setup line

    if Config.ade_rct_smis == ['', '']:
        QMessageBox.warning(obj,
                            "autodE-GUI",
                            "There are no reactants! Unable to write script.")
        return None
    for index, rct_smi in enumerate(Config.ade_rct_smis):
        if not rct_smi == '':
            gen_script.append(f"rct{index} = ade.Reactant(smiles='{rct_smi}')\n")
            rct_and_prod += f"rct{index},"

    if Config.ade_prod_smis == ['', '']:
        QMessageBox.warning(obj,
                            "autodE-GUI",
                            "There are no products! Unable to write script.")
        return None
    for index, prod_smi in enumerate(Config.ade_prod_smis):
        if not prod_smi == '':
            gen_script.append(f"prod{index} = ade.Product(smiles='{prod_smi}')\n")
            rct_and_prod += f"prod{index},"

    # calculation setup
    if Config.ade_job_typ == 'Reaction Profile':
        gen_script.append(f"rxn = ade.Reaction({rct_and_prod})\n")
        gen_script.append(f"rxn.calculate_reaction_profile()")

    # (over)write script
    try:
        with open(cwd + fname, 'w') as fh:
            fh.writelines(gen_script)
        # if success
        QMessageBox.information(obj, "autodE-GUI", "Finished generating script file")
        return None
    except OSError:
        QMessageBox.critical(obj,
                             "autodE-GUI",
                             "OSError: unable to writing file (file system error)")
        return None
