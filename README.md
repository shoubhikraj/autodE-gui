## autodE GUI

This is a Python-based GUI for the software autodE (which can model reaction pathways, deal with conformers of stationary points etc.
with quantum chemical methods)

At the moment, it is in development.

### How to install:

The code is written to be usable on all platforms. It has been tested on MacOS 12.6, Linux Ubuntu 22.0.4, and Windows 10.

The autodE GUI will work with Python >3.7, and needs the following packages:

1. PyQt5

If you are in a conda environment, install with `conda install pyqt` (This is PyQt5 as of October 2021, but the version
may change in future)

Otherwise, use `pip` to install it: `pip install pyqt5`.

2. It also needs a molecule editor for the "Draw" functionalities to work. adEGUI assumes that the editor is Avogadro, and 
will try to find it. The Avogadro executable must be available on PATH. (Any other editor will work. The editor needs to
have a command-line interface that can open and edit a mol file.)

### How to use
If PyQt5 is installed, run `python main.py`. This is start the main GUI window. Draw or type the reactants and products,
and once done, press the "Generate" button, which will write a script file for autodE ("aderun.py") in the current working
directory. You can then run the script file 
on a local system or submit to a cluster. In future, a facility to run autodE from the GUI will be implemented.

### Known issues
The molecule editor executable (default Avogadro) must be closed before the GUI interface will be responsive. This is
because it waits for Avogadro to finish editing the file, and prevents multiple instances of Avogadro attempting to edit
the same file.

### Development
autodE GUI is in development. Any feedback and advice is greatly appreciated.