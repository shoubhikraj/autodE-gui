# Main script for application, run it to get the GUI

import sys
from PyQt5.QtWidgets import QApplication
from adegui.main_win import MainWindow
from adegui.stylesheet import adegui_style_sheet

if __name__ == '__main__':
    app = QApplication([])
    app.setStyle("Fusion")
    #app.setStyleSheet(adegui_style_sheet)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
