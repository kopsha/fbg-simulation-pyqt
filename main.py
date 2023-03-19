#!python3

import PySide6.QtCore
import sys

from simulation.check_me import hello


def main(argv):
    print(hello(), f"{argv=}")
    print(PySide6.__version__)
    print(PySide6.QtCore.__version__)


if __name__ == "__main__":
    exit(main(sys.argv))
