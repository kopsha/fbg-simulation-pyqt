#!python3
import sys
from PySide6.QtWidgets import QApplication

from simulation.main_window import MainWindow


def main(argv):
    app = QApplication(argv)

    window = MainWindow()
    window.resize(1024, 768)
    window.show()

    # maybe call some post init stuff
    ret_code = app.exec()
    # maybe call some closing handlers

    return ret_code

if __name__ == "__main__":
    exit(main(sys.argv))
