import random
from PySide6.QtWidgets import (
    QWidget,
    QPushButton,
    QLabel,
    QHBoxLayout,
    QVBoxLayout,
    QLineEdit,
    QGridLayout,
    QTextEdit,
    QFileDialog,
    QCheckBox,
    QGroupBox,
    QRadioButton,
)
from PySide6 import QtCore
from PySide6.QtCore import Qt
from PySide6.QtGui import QTextCursor


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
        self.layout = QHBoxLayout(self)

        main_layout = self.make_main_layout()
        self.layout.addLayout(main_layout, 78)

        side_layout = self.make_side_layout()
        self.layout.addLayout(side_layout, 22)

        self.println("Aplicatia a fost initializata cu succes.")

    def make_side_layout(self):
        title = QLabel(self)
        title.setText("Jurnal mesaje")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.console = QTextEdit(self)
        self.console.setReadOnly(True)

        clear_button = QPushButton("Sterge jurnal", self)
        clear_button.clicked.connect(self.console.clear)

        layout = QVBoxLayout()
        layout.addWidget(title)
        layout.addWidget(self.console)
        layout.addWidget(clear_button)

        return layout

    def make_main_layout(self):
        grid = QGridLayout()

        section_1 = self.make_loader_section(1)
        grid.addLayout(section_1, 1, 1)

        section_2 = self.make_deform_types_section(2)
        grid.addLayout(section_2, 1, 2)

        section_3 = self.make_parameters_section(3)
        grid.addLayout(section_3, 2, 1)

        return grid

    def make_loader_section(self, section_id: int):
        title = QLabel(f"({section_id}) Incarca datele de tensiune si intindere")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)

        load_button = QPushButton("Alege fisier")
        load_button.clicked.connect(self.load_file)

        self.filepath = QLineEdit(self)
        self.filepath.setReadOnly(True)

        row = QHBoxLayout()
        row.addWidget(self.filepath)
        row.addWidget(load_button)

        si_units_group = QGroupBox("Daca distantele nu sunt exprimate in milimetri (mm)")
        use_si_units = QCheckBox("Aplica conversia din (m) in (mm)", si_units_group)
        group_layout = QVBoxLayout()
        group_layout.addWidget(use_si_units)
        si_units_group.setLayout(group_layout)

        layout = QVBoxLayout()
        layout.addWidget(title)
        layout.addLayout(row)
        layout.addWidget(si_units_group)
        layout.addStretch()

        return layout

    def make_deform_types_section(self, section_id: int):
        title = QLabel(f"({section_id}) Alege tipul simularii")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)

        strain_type_group = QGroupBox("Alege tipul de intindere")
        strain_group_layout = QVBoxLayout()
        no_strain = QRadioButton("Fara intindere", strain_type_group)
        no_strain.setChecked(True)
        uniform_strain = QRadioButton("Cu intindere longitudinala uniforma", strain_type_group)
        non_uniform_strain = QRadioButton("Cu intindere longitudinala neuniforma", strain_type_group)
        strain_group_layout.addWidget(no_strain)
        strain_group_layout.addWidget(uniform_strain)
        strain_group_layout.addWidget(non_uniform_strain)
        strain_type_group.setLayout(strain_group_layout)

        stress_type_group = QGroupBox("Alege tipul de tensiune")
        stress_group_layout = QVBoxLayout()
        no_stress = QRadioButton("Fara tensiune", stress_type_group)
        no_stress.setChecked(True)
        included_stress = QRadioButton("Cu tensiune transversala", stress_type_group)
        stress_group_layout.addWidget(no_stress)
        stress_group_layout.addWidget(included_stress)
        stress_type_group.setLayout(stress_group_layout)

        emulation_group = QGroupBox("Optiuni emulare")
        row1 = QHBoxLayout()
        has_emulate_temperature = QCheckBox("Emuleaza temperatura model", emulation_group)
        unit_label = QLabel("[K]", emulation_group)
        emulate_temperature = QLineEdit("293.15", emulation_group, alignment=Qt.AlignmentFlag.AlignRight)
        emulate_temperature.setEnabled(False)
        has_emulate_temperature.toggled.connect(emulate_temperature.setEnabled)
        row1.addWidget(has_emulate_temperature, stretch=3)
        row1.addWidget(unit_label)
        row1.addWidget(emulate_temperature)

        row2 = QHBoxLayout()
        has_host_expansion = QCheckBox(emulation_group)
        has_host_expansion.setText("Coeficientul de dilatatie termica (host)")
        unit_label = QLabel("[K<sup>-1</sup>]", emulation_group)
        host_expansion_coefficient = QLineEdit("5e-5", emulation_group, alignment=Qt.AlignmentFlag.AlignRight)
        host_expansion_coefficient.setEnabled(False)
        has_host_expansion.toggled.connect(host_expansion_coefficient.setEnabled)
        row2.addWidget(has_host_expansion, stretch=3)
        row2.addWidget(unit_label)
        row2.addWidget(host_expansion_coefficient)

        emulation_group_layout = QVBoxLayout()
        emulation_group_layout.addLayout(row1)
        emulation_group_layout.addLayout(row2)
        emulation_group.setLayout(emulation_group_layout)

        layout = QVBoxLayout()
        layout.addWidget(title)
        layout.addWidget(strain_type_group)
        layout.addWidget(stress_type_group)
        layout.addWidget(emulation_group)
        layout.addStretch()

        return layout

    def make_parameters_section(self, section_id: int):
        title = QLabel(f"({section_id}) Parametri simulare", alignment=Qt.AlignmentFlag.AlignCenter)

        row1, simulation_resolution = self.make_float_parameter("Rezoluția simulării", "[nm]", "0.05")
        row2, max_bandwidth = self.make_float_parameter("Lățime de bandă maximă", "[nm]", "1500.00")
        row3, min_bandwidth = self.make_float_parameter("Lățime de bandă minimă", "[nm]", "1600.00")
        row4, ambient_temperature = self.make_float_parameter("Temperatura ambientală", "[K]", "293.15")

        advanded_group = QGroupBox("Atributele fibrei (mod avansat)", checkable=True, checked=False)

        row5, initial_refractive_index = self.make_float_parameter("Indicele de refracție inițială", "[n<sub>eff</sub>]", "1.46")
        row6, mean_change_refractive_index = self.make_float_parameter("Variația medie în indicele de refracție", "[δn<sub>edd</sub>]", "4.5e-4")
        row7, fringe_visibility = self.make_float_parameter("Vizibilitatea franjelor", "%", "1.0")
        row8, directional_refractive_p11 = self.make_float_parameter("Constanta fotoelastică normală (Pockel)", "p<sub>11</sub>", "0.121")
        row9, directional_refractive_p12 = self.make_float_parameter("Constanta fotoelastică de tăiere (Pockel)", "p<sub>12</sub>", "0.270")
        row10, youngs_mod = self.make_float_parameter("Modulul de elasticitate (Young)", "[Pa]", "75e9")
        row11, poissons_coefficient = self.make_float_parameter("Coeficientul Poisson", "", "0.17")
        row12, fiber_expansion_coefficient = self.make_float_parameter("Coeficientul de dilatație termică", "[K<sup>-1</sup>]", "0.55e-6")
        row13, host_expansion_coefficient = self.make_float_parameter("Coeficientul de dilatație termică (host)", "[K<sup>-1</sup>]", "8.3e-6")

        advanced_group_layout = QVBoxLayout()
        advanced_group_layout.addLayout(row5)
        advanced_group_layout.addLayout(row6)
        advanced_group_layout.addLayout(row7)
        advanced_group_layout.addLayout(row8)
        advanced_group_layout.addLayout(row9)
        advanced_group_layout.addLayout(row10)
        advanced_group_layout.addLayout(row11)
        advanced_group_layout.addLayout(row12)
        advanced_group_layout.addLayout(row13)
        advanded_group.setLayout(advanced_group_layout)

        layout = QVBoxLayout()
        layout.addWidget(title)
        layout.addLayout(row1)
        layout.addLayout(row2)
        layout.addLayout(row3)
        layout.addLayout(row4)
        layout.addWidget(advanded_group)
        layout.addStretch()

        return layout

    def make_float_parameter(self, display_text: str, unit_text, value_text):
        row = QHBoxLayout()
        label = QLabel(display_text)
        unit_label = QLabel(unit_text)
        value = QLineEdit(value_text, alignment=Qt.AlignmentFlag.AlignRight)
        row.addWidget(label, stretch=3)
        row.addWidget(unit_label)
        row.addWidget(value)
        return row, value

    def println(self, text: str):
        if isinstance(text, str):
            show_text = text
        else:
            show_text = repr(text)
        self.console.append(show_text)
        self.console.moveCursor(QTextCursor.MoveOperation.End)

    def load_file(self):
        fullpath, _ = QFileDialog.getOpenFileName(self, "Incarca datele din", "./sample", "text (*.txt)")
        self.filepath.setText(fullpath)

    @QtCore.Slot()
    def on_click(self):
        self.text.setText(random.choice(self.hellos))
