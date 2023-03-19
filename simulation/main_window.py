import random
from PySide6.QtWidgets import QWidget, QPushButton, QLabel, QVBoxLayout, QStyle, QGridLayout
from PySide6 import QtCore
from PySide6.QtCore import Qt



class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.hellos = ["Hello world", "Hallo Welt", "Hei maailma", "Hola Mundo", "Привет мир"]

        icons = [
            res_name
            for res_name in dir(QStyle.StandardPixmap)
            if res_name.startswith("SP_")
        ]
        self.icon_grid = QGridLayout()

        ICON_COLS = 5
        for i, name in enumerate(icons):
            btn = QPushButton(name)
            pixmap = getattr(QStyle.StandardPixmap, name)
            icon = self.style().standardIcon(pixmap)
            btn.setIcon(icon)
            self.icon_grid.addWidget(btn, i//ICON_COLS, i%ICON_COLS)

        self.button = QPushButton("push")
        self.text = QLabel(self.hellos[0], alignment=Qt.AlignCenter)
        self.layout = QVBoxLayout(self)

        self.layout.addLayout(self.icon_grid)
        self.layout.addWidget(self.text)
        self.layout.addWidget(self.button)

        self.button.clicked.connect(self.on_click)

    @QtCore.Slot()
    def on_click(self):
        self.text.setText(random.choice(self.hellos))
