import numpy as np

from PySide6.QtWidgets import QDialog, QWidget, QHBoxLayout, QLabel, QVBoxLayout
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg as FigureCanvas,
)


class SpectrumView(QDialog):
    def __init__(self, parent: QWidget):
        super().__init__(parent)
        self.setWindowTitle(_("Plot FBG Spectrum"))
        self.setGeometry(512, 256, 1280, 768)

        self.setup_ui()

    def setup_ui(self):
        layout = QHBoxLayout(self)

        plot_layout = self.make_plot_layout()
        layout.addLayout(plot_layout, 78)

        side_layout = self.make_side_layout()
        layout.addLayout(side_layout, 22)

    def make_plot_layout(self):
        layout = QVBoxLayout()

        fig, ax = plt.subplots()
        x = np.linspace(0, 10, 100)
        ax.plot(x, np.sin(x), label="sin(x)")
        ax.plot(x, np.cos(x), label="cos(x)")
        ax.set_title("Sin and Cos Plots")
        ax.legend()

        self.canvas = FigureCanvas(fig)
        self.canvas.draw()

        layout.addWidget(self.canvas)

        return layout

    def make_side_layout(self):
        layout = QVBoxLayout()
        layout.addWidget(QLabel("Options"))
        return layout
