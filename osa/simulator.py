"""
Simulation of reflected FBG spectrum using coupled-mode theory.
TODO: review with [Ben Frey](https://github.com/benfrey)
"""
import numpy as np


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

        """

        self.data = None
        from_meters = True
        self.fbg_count = fbg_count
        self.fbg_length = fbg_length
        self.fbg_positions = fbg_positions
        self.tolerance = tolerance

    def from_file(self, filepath: str, from_meters=True):
        """
        Loads FBG data from text file expecting the following fiels:
            - 'x', 'LE11', 'LE22', 'LE33', 'S11', 'S22', 'S33', 'T'

        :param str filepath: full path to the data file
        :param bool from_meters: enable conversion from meters

        NB: Everything in this class will be expressed in (mm)
        """

        fields = dict(
            x=dict(format="f8", si_conv=1000 if from_meters else 1),
            LE11=dict(format="f8", si_conv=1),
            LE22=dict(format="f8", si_conv=1),
            LE33=dict(format="f8", si_conv=1),
            S11=dict(format="f8", si_conv=1 / 10**6 if from_meters else 1),
            S22=dict(format="f8", si_conv=1 / 10**6 if from_meters else 1),
            S33=dict(format="f8", si_conv=1 / 10**6 if from_meters else 1),
            T=dict(format="f8", si_conv=1),
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
                conv = fields[key]["si_conv"]
                data[b_key][key] = conv * raw_data[key]

        self.data = data
        return data
