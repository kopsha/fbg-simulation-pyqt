from PySide6.QtCore import QThread, Signal
from osa.simulator import OsaSimulator


class WorkerThread(QThread):
    progress = Signal(int)

    def __init__(self, params: dict) -> None:
        super().__init__()

        self.units = params.pop("units")
        self.include_undeformed_signal = params.pop("has_reflected_signal")
        self.strain_type = params.pop("strain_type")
        self.stress_type = params.pop("stress_type")
        self.datafile = params.pop("filepath")
        self.params = params
        self.error_message = ""
        self.data = None

    def run(self):
        self.progress.emit(13)

        undeformed_data = None
        deformed_data = None
        summary_data = None

        try:
            simu = OsaSimulator(**self.params)
            simu.from_file(filepath=self.datafile, units=self.units)
            self.progress.emit(21)

            if self.include_undeformed_signal:
                undeformed_data = simu.undeformed_fbg()
                self.progress.emit(34)

            self.progress.emit(55)
            deformed_data = simu.deformed_fbg(
                strain_type=self.strain_type,
                stress_type=self.stress_type,
            )
            self.progress.emit(89)

            summary_data = simu.compute_fbg_shifts_and_widths(
                strain_type=self.strain_type,
                stress_type=self.stress_type,
            )
            self.progress.emit(100)

        except Exception as err:
            self.error_message = str(err)

        self.data = dict(
            undeformed=undeformed_data,
            deformed=deformed_data,
            summary=summary_data,
        )
