"""
Simulation of reflected FBG spectrum using coupled-mode theory.
TODO: review with [Ben Frey](https://github.com/benfrey)
"""
import numpy as np
from enum import Enum, auto
from collections import defaultdict


class SiUnits(Enum):
    METERS = auto()
    MILLIMETERS = auto()


class OsaSimulator:
    def __init__(
        self, fbg_count: int, fbg_length: float, tolerance: float, fbg_positions: list
    ):
        """
        Prepares an OSA simulation, with following parameters:

        :param int fbg_count: Number of FBG per optical fibre
        :param float fbg_length: Length of the Gratting
        :param float tolerance: Tolerance in the FBG length
        :param list positions: Position of each FBG from the begginng of the Path

        NB: Everything in this class will be expressed in (mm)
        """
        self.data = None
        self.fbg_count = fbg_count
        self.fbg_length = fbg_length
        self.fbg_positions = fbg_positions
        self.tolerance = tolerance

    def from_file(self, filepath: str, units=SiUnits.MILLIMETERS):
        """
        Loads FBG data from text file expecting the following fiels:
            - 'x', 'LE11', 'LE22', 'LE33', 'S11', 'S22', 'S33', 'T'

        :param str filepath: full path to the data file
        :param bool from_meters: enable conversion from meters
        """
        from_meters = units == SiUnits.METERS
        fields = dict(
            x=defaultdict(format="f8", multiply=1000 if from_meters else 1),
            LE11=defaultdict(format="f8"),
            LE22=defaultdict(format="f8"),
            LE33=defaultdict(format="f8"),
            S11=defaultdict(format="f8", divide=(10**6) if from_meters else 1),
            S22=defaultdict(format="f8", divide=(10**6) if from_meters else 1),
            S33=defaultdict(format="f8", divide=(10**6) if from_meters else 1),
            T=defaultdict(format="f8"),
        )
        dtypes = dict(
            names=list(fields.keys()), formats=[v["format"] for v in fields.values()]
        )
        raw_data = np.genfromtxt(filepath, dtypes, comments="%")

        data = dict()
        for b in range(self.fbg_count):
            b_key = f"FBG{b+1}"
            data[b_key] = dict()
            for key in fields:
                if key.startswith("S"):
                    factor = fields[key]["divide"]
                    data[b_key][key] = raw_data[key] / factor
                elif key == "x":
                    factor = fields[key]["multiply"]
                    data[b_key][key] = raw_data[key] * factor
                else:
                    data[b_key][key] = raw_data[key]

        self.data = data
        return data
