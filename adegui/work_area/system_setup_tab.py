from PyQt5.QtCore import pyqtSlot, Qt
from PyQt5.QtWidgets import (QWidget, QSpinBox, QHBoxLayout,
                             QVBoxLayout, QGroupBox, QLineEdit)


class SystemSetupTab(QGroupBox):
    def __init__(self):
        super().__init__("System Settings")

        max_core_box = QLineEdit()
        n_cores_dial = QSpinBox()

        upper_bar_layout = QHBoxLayout()
        upper_bar_layout.addWidget(max_core_box)
        upper_bar_layout.addStretch()
        upper_bar = QWidget()
        upper_bar.setLayout(upper_bar_layout)

        lower_bar_layout = QHBoxLayout()
        lower_bar_layout.addWidget(n_cores_dial)
        lower_bar_layout.addStretch()
        lower_bar = QWidget()
        lower_bar.setLayout(lower_bar_layout)

        large_layout = QVBoxLayout()
        large_layout.addWidget(upper_bar)
        large_layout.addWidget(lower_bar)
        self.setLayout(large_layout)
