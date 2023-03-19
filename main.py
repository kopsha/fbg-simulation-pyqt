#!/usr/bin/env python3
from PyQt5.QtWidgets import QApplication
import sys

from simulation.check_me import hello


def main(argv):
    print(f"{argv=}", hello())
    app = QApplication(argv)
    ret_code = app.exec_()

    if ret_code:
        print("Something went wrong")
    
    return ret_code


if __name__ == "__main__":
    exit(main(sys.argv))
