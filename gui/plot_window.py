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
    QFileDialog,
)
from PySide6.QtGui import QStandardItem, QStandardItemModel
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg as FigureCanvas,
)


class SpectrumView(QDialog):
    USE_COLORS = "orange,red,blue,green,purple,black".split(",")

    def __init__(self, parent: QWidget, data: dict):
        super().__init__(parent)
        self.setWindowTitle(_("Plot FBG Spectrum"))
        self.setGeometry(512, 256, 1280, 768)

        self.data = data

        self.setup_ui()

    def setup_ui(self):
        layout = QHBoxLayout(self)

        side_layout = self.make_side_layout()
        plot_layout = self.make_plot_layout()

        layout.addLayout(plot_layout, 68)
        layout.addLayout(side_layout, 32)

        self.refresh_plot()

    def make_plot_layout(self):
        layout = QVBoxLayout()

        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvas(self.fig)

        layout.addWidget(self.canvas)
        return layout

    def make_side_layout(self):
        options_group = self.make_options_group()
        axis_group = self.make_axis_group()
        line_group = self.make_line_group()

        summary_group = self.make_summary_group()

        save_button = QPushButton(_("Save plot figure"))
        save_button.clicked.connect(self.savePlotFigure)

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

        self.legend = QCheckBox(_("Legend"), options_group, checked=True)
        self.grid = QCheckBox(_("Grid"), options_group, checked=False)
        self.split = QCheckBox(_("Split waveform (Birefringence)"), options_group, checked=False)

        self.legend.clicked.connect(self.refresh_plot)
        self.grid.clicked.connect(self.refresh_plot)
        self.split.clicked.connect(self.refresh_plot)

        options_group_layout.addWidget(self.legend)
        options_group_layout.addWidget(self.grid)
        options_group_layout.addWidget(self.split)
        options_group.setLayout(options_group_layout)

        return options_group

    def make_axis_group(self):
        axis_group = QGroupBox(_("Axis control"))
        axis_group_layout = QGridLayout()

        x_between = QLabel("↔", axis_group)
        x_label = QLabel(_("Wavelength"), axis_group)
        self.xmin = QSpinBox(
            axis_group,
            minimum=1000,
            maximum=2000,
            singleStep=10,
            value=1490,
            suffix=" nm",
            alignment=Qt.AlignmentFlag.AlignRight,
        )
        self.xmax = QSpinBox(
            axis_group,
            minimum=1000,
            maximum=2000,
            singleStep=10,
            value=1610,
            suffix=" nm",
            alignment=Qt.AlignmentFlag.AlignRight,
        )

        y_between = QLabel("↔", axis_group)
        y_label = QLabel(_("Reflectivity"), axis_group)
        self.ymin = QDoubleSpinBox(
            axis_group,
            minimum=-0.5,
            maximum=2.0,
            singleStep=0.1,
            value=-0.1,
            suffix=" R",
            alignment=Qt.AlignmentFlag.AlignRight,
        )
        self.ymax = QDoubleSpinBox(
            axis_group,
            minimum=-0.5,
            maximum=2.0,
            singleStep=0.1,
            value=1.1,
            suffix=" R",
            alignment=Qt.AlignmentFlag.AlignRight,
        )

        self.xmin.valueChanged.connect(self.refresh_plot)
        self.xmax.valueChanged.connect(self.refresh_plot)
        self.ymin.valueChanged.connect(self.refresh_plot)
        self.ymax.valueChanged.connect(self.refresh_plot)

        axis_group_layout.addWidget(x_label, 0, 0)
        axis_group_layout.addWidget(self.xmin, 1, 0)
        axis_group_layout.addWidget(x_between, 1, 1)
        axis_group_layout.addWidget(self.xmax, 1, 2)

        axis_group_layout.addWidget(y_label, 2, 0)
        axis_group_layout.addWidget(self.ymin, 3, 0)
        axis_group_layout.addWidget(y_between, 3, 1)
        axis_group_layout.addWidget(self.ymax, 3, 2)

        axis_group_layout.setColumnStretch(0, 3)
        axis_group_layout.setColumnStretch(2, 3)
        axis_group.setLayout(axis_group_layout)

        return axis_group

    def make_line_group(self):
        line_group = QGroupBox(_("Line properties"))
        line_group_layout = QGridLayout()

        line_width_label = QLabel(_("Width"), line_group)
        self.line_width = QDoubleSpinBox(
            line_group,
            minimum=0.1,
            maximum=5.0,
            singleStep=0.1,
            value=0.5,
            suffix=" pt",
            alignment=Qt.AlignmentFlag.AlignRight,
        )
        undeformed_color_label = QLabel(_("Undeformed color"), line_group)
        self.undeformed_color = QComboBox(line_group)
        self.undeformed_color.addItems(self.USE_COLORS)
        self.undeformed_color.setCurrentText("blue")

        simulated_color_label = QLabel(_("Simulated color"), line_group)
        self.simulated_color = QComboBox(line_group)
        self.simulated_color.addItems(self.USE_COLORS)
        self.simulated_color.setCurrentText("orange")

        self.line_width.valueChanged.connect(self.refresh_plot)
        self.undeformed_color.currentTextChanged.connect(self.refresh_plot)
        self.simulated_color.currentTextChanged.connect(self.refresh_plot)

        line_group_layout.addWidget(line_width_label, 0, 0)
        line_group_layout.addWidget(self.line_width, 0, 1)
        line_group_layout.addWidget(undeformed_color_label, 1, 0)
        line_group_layout.addWidget(self.undeformed_color, 1, 1)
        line_group_layout.addWidget(simulated_color_label, 2, 0)
        line_group_layout.addWidget(self.simulated_color, 2, 1)
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

        self.read_data_rows(model)

        summary_group_layout.addWidget(summary_table)
        summary_group.setLayout(summary_group_layout)

        return summary_group

    def read_data_rows(self, model):
        for i in range(self.data["params"]["fbg_count"]):
            key = f"FBG{i+1}"
            wave_length = "{} nm".format(self.data["params"]["original_wavelengths"][i])
            wave_shift = "{:.15f} nm".format(self.data["summary"][key]["wave_shift"])
            wave_width = "{:.15f} nm".format(self.data["summary"][key]["wave_width"])
            self.add_summary_row(model, [wave_length, wave_shift, wave_width])

    def add_summary_row(self, model, data):
        items = [QStandardItem(str(datum)) for datum in data]
        for item in items:
            item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        model.appendRow(items)

    def refresh_plot(self):
        self.ax.clear()

        if self.data["undeformed"]:
            self.ax.plot(
                self.data["undeformed"]["wavelength"],
                self.data["undeformed"]["reflec"],
                color=self.undeformed_color.currentText(),
                linewidth=self.line_width.value(),
                label=_("Undeformed FBG Spectrum"),
            )

        self.ax.plot(
            self.data["deformed"]["wavelength"],
            self.data["deformed"]["reflec"],
            color=self.simulated_color.currentText(),
            linewidth=self.line_width.value(),
            label=_("Deformed FBG Spectrum"),
        )

        if self.split.isChecked():
            self.ax.plot(
                self.data["deformed"]["Y_split"]["wavelength"],
                self.data["deformed"]["Y_split"]["reflec"],
                color="red",
                linewidth=self.line_width.value(),
                label=_("Y-Wave Contribution"),
            )
            self.ax.plot(
                self.data["deformed"]["Z_split"]["wavelength"],
                self.data["deformed"]["Z_split"]["reflec"],
                color="green",
                linewidth=self.line_width.value(),
                label=_("Z-Wave Contribution"),
            )

        if self.legend.isChecked():
            self.ax.legend()

        if self.grid.isChecked():
            self.ax.grid()

        self.ax.set_xlabel("{} [nm]".format(_("Wavelength")))
        self.ax.set_ylabel("{} [R]".format(_("Reflectivity")))
        self.ax.set_xlim(xmin=self.xmin.value(), xmax=self.xmax.value())
        self.ax.set_ylim(ymin=self.ymin.value(), ymax=self.ymax.value())
        self.canvas.draw()

    def savePlotFigure(self):
        to_file, filter = QFileDialog.getSaveFileName(
            self,
            caption=_("Save FBG Spectrum Plot"),
            dir="./sample",
            filter="{} (*.png)".format(_("Images"))
        )
        if to_file:
            self.fig.savefig(to_file)
