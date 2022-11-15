# defines the buttons in the main window
from PyQt5.QtWidgets import (QHBoxLayout, QPushButton, QWidget,
                             qApp, QStyle, QStackedWidget)



class MainButtons(QWidget):
    def __init__(self):
        super().__init__()

        # create the buttons
        btn1 = QPushButton("Run")
        btn1.setEnabled(False)  # implement run function later
        self.btn2 = QPushButton("Back")

        self.btn3 = QStackedWidget()  # button 3 has to change at the last tab
        self.next_btn = QPushButton("Next")
        self.generate_btn = QPushButton("Generate")
        self.btn3.addWidget(self.next_btn)
        self.btn3.addWidget(self.generate_btn)

        btn4 = QPushButton("Exit")

        # add icons to buttons
        btn1.setIcon(qApp.style().standardIcon(QStyle.SP_MediaPlay))
        self.btn2.setIcon(qApp.style().standardIcon(QStyle.SP_ArrowLeft))
        self.next_btn.setIcon(qApp.style().standardIcon(QStyle.SP_ArrowRight))
        self.generate_btn.setIcon(qApp.style().standardIcon(QStyle.SP_DialogSaveButton))
        btn4.setIcon(qApp.style().standardIcon(QStyle.SP_BrowserStop))

        # connect close to trigger function
        btn4.clicked.connect(qApp.quit)
        # other functions will be connected in the main window
        # so that the tabs can be changed

        # the buttons have to be stacked horizontally
        btn_area_layout = QHBoxLayout()
        btn_area_layout.addWidget(btn1)  # first button
        btn_area_layout.addWidget(self.btn2)  # back button
        btn_area_layout.addWidget(self.btn3)  # forward or generate button
        btn_area_layout.addStretch()
        btn_area_layout.addWidget(btn4)  # third button

        self.setLayout(btn_area_layout)




