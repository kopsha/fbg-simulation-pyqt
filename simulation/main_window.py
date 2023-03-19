import random
from PySide6.QtWidgets import (
    QWidget,
    QPushButton,
    QLabel,
    QVBoxLayout,
    QStyle,
    QGridLayout,
)
from PySide6 import QtCore
from PySide6.QtCore import Qt


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.hellos = [
            "Hello world",
            "Hallo Welt",
            "Hei maailma",
            "Hola Mundo",
            "Привет мир",
        ]
        self.setup_ui()

    def setup_ui(self):
        self.layout = QVBoxLayout(self)

        self.icon_grid = self.make_icons_grid_layout()
        self.layout.addLayout(self.icon_grid)

        self.text = QLabel(self.hellos[0], alignment=Qt.AlignCenter)
        self.layout.addWidget(self.text)

        self.button = QPushButton("push")
        self.layout.addWidget(self.button)
        self.button.clicked.connect(self.on_click)

    def make_icons_grid_layout(self):
        icon_grid = QGridLayout()

        ICON_COLS = 5
        standard_icons = (
            si for si in QStyle.StandardPixmap if si.name.startswith("SP_")
        )
        for i, sp in enumerate(standard_icons):
            icon = self.style().standardIcon(sp)
            label = QLabel(self)
            label.setPixmap(icon.pixmap(32))
            icon_grid.addWidget(label, i // ICON_COLS, i % ICON_COLS)

        return icon_grid

    @QtCore.Slot()
    def on_click(self):
        self.text.setText(random.choice(self.hellos))
