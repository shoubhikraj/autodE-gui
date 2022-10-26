# defines the buttons in the main window
from PyQt5.QtWidgets import QHBoxLayout, QPushButton, QWidget, qApp
from PyQt5.QtCore import pyqtSlot
from adegui.generate_script import write_ade_script_from_config


class MainButtons(QWidget):
    def __init__(self):
        super().__init__()

        # create the buttons
        btn1 = QPushButton("Run")
        btn1.setEnabled(False)  # implement run function later
        btn2 = QPushButton("Generate")
        btn3 = QPushButton("Exit")

        # connect buttons to trigger functions
        btn2.clicked.connect(self.generate_script)
        btn3.clicked.connect(qApp.quit)

        # the buttons have to be stacked horizontally
        btn_area_layout = QHBoxLayout()
        btn_area_layout.addWidget(btn1)  # first button
        btn_area_layout.addWidget(btn2)  # second button
        btn_area_layout.addWidget(btn3)  # third button

        self.setLayout(btn_area_layout)

    @pyqtSlot()
    def generate_script(self):
        write_ade_script_from_config(self)


