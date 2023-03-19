import random
from PySide6 import QtWidgets, QtCore



class MainWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.hellos = ["Hello world", "Hallo Welt", "Hei maailma", "Hola Mundo", "Привет мир"]

        self.button = QtWidgets.QPushButton("push")
        self.text = QtWidgets.QLabel(self.hellos[0], alignment=QtCore.Qt.AlignCenter)
        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.addWidget(self.text)
        self.layout.addWidget(self.button)

        self.button.clicked.connect(self.on_click)

    @QtCore.Slot()
    def on_click(self):
        self.text.setText(random.choice(self.hellos))
