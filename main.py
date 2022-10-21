# Main script for application, run it to get the GUI

import sys
from PyQt5.QtWidgets import QApplication
from adegui.main_win import MainWindow

if __name__ == '__main__':
    app = QApplication([])
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
