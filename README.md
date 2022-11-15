## autodE GUI

This is a Python-based GUI for the software autodE (which can model reaction pathways, deal with conformers of stationary points etc.
with quantum chemical methods)

At the moment, it is in development and bugs/hangups are expected.

### How to install:

The code is written to be usable on all platforms. It has been tested on MacOS 12.6, Linux Ubuntu 22.0.4, and Windows 10.

The autodE GUI will work with Python >3.7, and needs the following packages:

1. PyQt5

If you are in a conda environment, install with `conda install pyqt` (This is PyQt5 as of October 2021, but the version
may change in future)

Otherwise, use `pip` to install it: `pip install pyqt5`.

2. RDKit

In a conda environment, use `conda install -c conda-forge rdkit`. You can also use pip (as of 2021) to install it: `pip 
install rdkit`.

3. Avogadro (or alternative):

It also needs a molecule editor for the "Draw" functionalities to work. adEGUI assumes that the editor is Avogadro, and 
will try to find it. The Avogadro executable must be available on PATH. (Any other editor should also work but have not
been tested so far. The editor needs to have a command-line interface that can open and edit a mol file.) The GUI
will work even if no editor exists, but then molecules can be input as SMILES strings only.

Avogadro or other similar softwares are usually installed through package managers or installers in the normal way (i.e.
it is not related to python or conda directly).

### How to use
If PyQt5 is installed, run `python main.py`. This will start the main GUI window. Draw or type the reactants and products,
and once done, press the "Generate" button, which will write a script file for autodE. You can then run the script file 
on a local system or submit to a cluster. In future, a feature to run autodE from the GUI will be implemented.

### Known issues
The molecule editor executable (default Avogadro) must be closed before the GUI interface will be responsive. This is
because it waits for Avogadro to finish editing the file, and prevents multiple instances of Avogadro attempting to edit
the same file.

### Development
autodE GUI is currently in development. Any feedback and advice is greatly appreciated.