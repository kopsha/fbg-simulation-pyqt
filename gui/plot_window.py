import numpy as np

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog,
    QWidget,
    QHBoxLayout,
    QLabel,
    QVBoxLayout,
    QGroupBox,
    QCheckBox,
    QGridLayout,
    QSpinBox,
    QDoubleSpinBox,
    QComboBox,
    QPushButton,
    QTableView,
    QHeaderView,
)
from PySide6.QtGui import QStandardItem, QStandardItemModel
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
        layout.addLayout(plot_layout, 68)

        side_layout = self.make_side_layout()
        layout.addLayout(side_layout, 32)

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
        options_group = self.make_options_group()
        axis_group = self.make_axis_group()
        line_group = self.make_line_group()

        summary_group = self.make_summary_group()

        save_button = QPushButton(_("Save picture"))

        layout = QVBoxLayout()
        layout.addWidget(options_group)
        layout.addWidget(axis_group)
        layout.addWidget(line_group)
        layout.addWidget(summary_group)
        layout.addStretch()
        layout.addWidget(save_button)

        return layout

    def make_options_group(self):
        options_group = QGroupBox(_("Plot options"))
        options_group_layout = QVBoxLayout()

        legend = QCheckBox(_("Legend"), options_group, checked=True)
        grid = QCheckBox(_("Grid"), options_group, checked=False)
        split = QCheckBox(_("Split waveform (Birefringence)"), options_group, checked=False)

        options_group_layout.addWidget(legend)
        options_group_layout.addWidget(grid)
        options_group_layout.addWidget(split)
        options_group.setLayout(options_group_layout)

        return options_group

    def make_axis_group(self):
        axis_group = QGroupBox(_("Axis control"))
        axis_group_layout = QGridLayout()

        x_between = QLabel("↔", axis_group)
        x_label = QLabel(_("Wavelength"), axis_group)
        x_min = QSpinBox(
            axis_group,
            minimum=1000,
            maximum=2000,
            singleStep=10,
            value=1500,
            suffix=" nm",
            alignment=Qt.AlignmentFlag.AlignRight,
        )
        x_max = QSpinBox(
            axis_group,
            minimum=1000,
            maximum=2000,
            singleStep=10,
            value=1600,
            suffix=" nm",
            alignment=Qt.AlignmentFlag.AlignRight,
        )

        y_between = QLabel("↔", axis_group)
        y_label = QLabel(_("Reflectivity"), axis_group)
        y_min = QDoubleSpinBox(
            axis_group,
            minimum=-0.5,
            maximum=2.0,
            singleStep=0.1,
            value=0.0,
            suffix=" R",
            alignment=Qt.AlignmentFlag.AlignRight,
        )
        y_max = QDoubleSpinBox(
            axis_group,
            minimum=-0.5,
            maximum=2.0,
            singleStep=0.1,
            value=1.0,
            suffix=" R",
            alignment=Qt.AlignmentFlag.AlignRight,
        )

        axis_group_layout.addWidget(x_label, 0, 0)
        axis_group_layout.addWidget(x_min, 1, 0)
        axis_group_layout.addWidget(x_between, 1, 1)
        axis_group_layout.addWidget(x_max, 1, 2)

        axis_group_layout.addWidget(y_label, 2, 0)
        axis_group_layout.addWidget(y_min, 3, 0)
        axis_group_layout.addWidget(y_between, 3, 1)
        axis_group_layout.addWidget(y_max, 3, 2)

        axis_group_layout.setColumnStretch(0, 3)
        axis_group_layout.setColumnStretch(2, 3)
        axis_group.setLayout(axis_group_layout)

        return axis_group

    def make_line_group(self):
        line_group = QGroupBox(_("Line properties"))
        line_group_layout = QGridLayout()

        line_width_label = QLabel(_("Width"), line_group)
        line_width = QDoubleSpinBox(
            line_group,
            minimum=0.1,
            maximum=5.0,
            singleStep=0.1,
            value=0.5,
            suffix=" pt",
            alignment=Qt.AlignmentFlag.AlignRight,
        )
        underformed_color_label = QLabel(_("Undeformed color"), line_group)
        underformed_color = QComboBox(line_group)
        underformed_color.addItems("red,orange,blue,green,black".split(","))

        simulated_color_label = QLabel(_("Simulated color"), line_group)
        simulated_color = QComboBox(line_group)
        simulated_color.addItems("red,orange,blue,green,black".split(","))

        line_group_layout.addWidget(line_width_label, 0, 0)
        line_group_layout.addWidget(line_width, 0, 1)
        line_group_layout.addWidget(underformed_color_label, 1, 0)
        line_group_layout.addWidget(underformed_color, 1, 1)
        line_group_layout.addWidget(simulated_color_label, 2, 0)
        line_group_layout.addWidget(simulated_color, 2, 1)
        line_group.setLayout(line_group_layout)

        return line_group

    def make_summary_group(self):
        summary_group = QGroupBox(_("FBG Output"))
        summary_group_layout = QVBoxLayout()

        summary_table = QTableView()

        model = QStandardItemModel()
        summary_table.setModel(model)

        columns = [_("Wavelength"), _("Shift"), _("Width Variation")]
        model.setHorizontalHeaderLabels(columns)
        horizontal_header = summary_table.horizontalHeader()
        horizontal_header.setSectionResizeMode(QHeaderView.ResizeToContents)

        self.add_summary_row(model, ["15cm", "2mm", "1mm"])
        self.add_summary_row(model, ["20cm", "3mm", "1.5mm"])

        summary_group_layout.addWidget(summary_table)
        summary_group.setLayout(summary_group_layout)

        return summary_group

    def add_summary_row(self, model, data):
        items = [QStandardItem(str(datum)) for datum in data]
        model.appendRow(items)
