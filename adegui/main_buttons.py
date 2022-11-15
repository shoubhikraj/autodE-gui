# defines the buttons in the main window
from PyQt5.QtWidgets import (QHBoxLayout, QPushButton, QWidget,
                             qApp, QStyle)
from PyQt5.QtCore import pyqtSlot
from adegui.generate_script import write_ade_script_from_config


class MainButtons(QWidget):
    def __init__(self):
        super().__init__()

        # create the buttons
        btn1 = QPushButton("Run")
        btn1.setEnabled(False)  # implement run function later
        btn2 = QPushButton("Generate")  # TODO put back and next buttons
        btn3 = QPushButton("Exit")

        # add icons to buttons
        btn1.setIcon(qApp.style().standardIcon(QStyle.SP_MediaPlay))
        btn2.setIcon(qApp.style().standardIcon(QStyle.SP_DialogSaveButton))
        btn3.setIcon(qApp.style().standardIcon(QStyle.SP_BrowserStop))

        # connect buttons to trigger functions
        btn2.clicked.connect(self.generate_script)
        btn3.clicked.connect(qApp.quit)

        # the buttons have to be stacked horizontally
        btn_area_layout = QHBoxLayout()
        btn_area_layout.addWidget(btn1)  # first button
        btn_area_layout.addWidget(btn2)  # second button
        btn_area_layout.addStretch()
        btn_area_layout.addWidget(btn3)  # third button

        self.setLayout(btn_area_layout)

    @pyqtSlot()
    def generate_script(self):
        write_ade_script_from_config(self)


