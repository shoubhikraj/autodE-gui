from PyQt5.QtCore import pyqtSlot, Qt
from PyQt5.QtWidgets import (QWidget, QSpinBox, QHBoxLayout,
                             QVBoxLayout, QGroupBox, QDoubleSpinBox,
                             QLabel)
from adegui import Config


class SystemSetupTab(QGroupBox):
    def __init__(self):
        super().__init__("System Settings")

        self.n_cores_dial = QSpinBox()
        self.n_cores_dial.setRange(1, 200)
        self.n_cores_dial.setValue(Config.ade_n_cores)
        self.n_cores_dial.valueChanged.connect(self.n_cores_changed)
        self.n_cores_dial.valueChanged.connect(self.update_explanation_label)

        self.max_core_dial = QDoubleSpinBox()
        self.max_core_dial.setRange(0.5, 20.0)
        self.max_core_dial.setSingleStep(0.5)
        self.max_core_dial.setValue(Config.ade_max_core_mem / 1024.0)
        self.max_core_dial.valueChanged.connect(self.max_core_mem_changed)
        self.max_core_dial.valueChanged.connect(self.update_explanation_label)

        self.mem_explanation = QLabel()
        self.update_explanation_label()  # initialize

        upper_bar_layout = QHBoxLayout()
        upper_bar_layout.addWidget(QLabel("Number of cores to use: "))
        upper_bar_layout.addWidget(self.n_cores_dial)
        upper_bar_layout.addStretch()
        upper_bar = QWidget()
        upper_bar.setLayout(upper_bar_layout)

        lower_bar_layout = QHBoxLayout()
        lower_bar_layout.addWidget(QLabel("Memory per core (GB): "))
        lower_bar_layout.addWidget(self.max_core_dial)
        lower_bar_layout.addWidget(self.mem_explanation)
        lower_bar_layout.addStretch()
        lower_bar = QWidget()
        lower_bar.setLayout(lower_bar_layout)

        large_layout = QVBoxLayout()
        large_layout.addWidget(upper_bar)
        large_layout.addWidget(lower_bar)
        large_layout.addStretch()
        self.setLayout(large_layout)

    @pyqtSlot()
    def max_core_mem_changed(self):
        # convert to GB
        Config.ade_max_core_mem = self.max_core_dial.value() * 1024.0

    @pyqtSlot()
    def n_cores_changed(self):
        Config.ade_n_cores = self.n_cores_dial.value()

    @pyqtSlot()
    def update_explanation_label(self):
        """
        The label shows the total amount of memory that
        would be used by autodE given the currently chosen
        n_cores and max_core
        """
        label = (f"  Total required memory = "
                 f"{self.n_cores_dial.value()} X "
                 f"{self.max_core_dial.value():.2f} GB = "
                 f"{self.max_core_dial.value()*self.n_cores_dial.value():.2f} GB")
        self.mem_explanation.setText(label)

