from PyQt5.QtWidgets import QMainWindow, QVBoxLayout, QWidget
from adegui.work_area.work_area import WorkAreaTabs
from adegui.main_buttons import MainButtons


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("autodE GUI")

        # The layout for the main window => buttons at bottom,
        # tabs and work area at top
        main_layout = QVBoxLayout()  # create a vertical stacked layout
        main_layout.addWidget(WorkAreaTabs())  # first part of vertical layout <= main work area
        main_layout.addWidget(MainButtons())  # second part of vertical layout <= buttons

        # QMainWindow already has a fixed layout, so create a widget then
        # apply layout on it, and set it as the central widget
        main_widget = QWidget()
        main_widget.setLayout(main_layout)
        self.setCentralWidget(main_widget)
