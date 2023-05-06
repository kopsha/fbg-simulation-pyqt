"""
Simulation of reflected FBG spectrum using coupled-mode theory.
TODO: review with [Ben Frey](https://github.com/benfrey)
"""
import numpy as np
from cmath import pi, sqrt, cosh, sinh
from enum import Enum, auto


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
        self.fbg_count = int(fbg_count)
        self.fbg_length = np.float64(fbg_length)
        self.tolerance = np.float64(tolerance)
        self.fbg_positions = np.array(fbg_positions, dtype=np.float64)

    def from_file(self, filepath: str, units=SiUnits.MILLIMETERS):
        """
        Loads FBG data from text file expecting the following fiels:
            - 'x', 'LE11', 'LE22', 'LE33', 'S11', 'S22', 'S33', 'T'

        :param str filepath: full path to the data file
        :param bool from_meters: enable conversion from meters
        """
        from_meters = units == SiUnits.METERS
        fields = dict(
            x=dict(format="f8", multiply=1000 if from_meters else 1),
            LE11=dict(format="f8"),
            LE22=dict(format="f8"),
            LE33=dict(format="f8"),
            S11=dict(format="f8", divide=(10**6) if from_meters else 1),
            S22=dict(format="f8", divide=(10**6) if from_meters else 1),
            S33=dict(format="f8", divide=(10**6) if from_meters else 1),
            T=dict(format="f8"),
        )
        dtypes = dict(
            names=list(fields.keys()), formats=[v["format"] for v in fields.values()]
        )
        raw_data = np.genfromtxt(filepath, dtypes, comments="%")

        fbg = dict()  # Former FBG array (it holds fbg data after all)
        for b in range(self.fbg_count):
            b_key = f"FBG{b+1}"
            fbg[b_key] = dict()
            for key in fields:
                if key.startswith("S"):
                    factor = fields[key]["divide"]
                    fbg[b_key][key] = raw_data[key] / factor
                elif key == "x":
                    factor = fields[key]["multiply"]
                    fbg[b_key][key] = raw_data[key] * factor
                else:
                    fbg[b_key][key] = raw_data[key]

        self.fbg = fbg
        return fbg

    def sigma(self, period, wavelen):
        first = (1.0 / wavelen) - (1.0 / (2.0 * self.initial_refractive_index * period))
        second = 2.0 * pi * self.mean_change_refractive_index / wavelen
        result = 2.0 * pi * self.initial_refractive_index * first + second
        return result

    def kaa(self, wavelen):
        result = (
            pi * self.fringe_visibility * self.mean_change_refractive_index / wavelen
        )
        return result

    def transfer_matrix(self, period, wavelen):
        f1 = np.identity(2)
        M = 20
        deltz = self.fbg_length * (10**6) / M
        for _ in range(M):
            sig = self.sigma(period, wavelen=wavelen)
            kaa = self.kaa(wavelen=wavelen)
            gamma = sqrt(kaa**2 - sig**2)

            f11 = complex(cosh(gamma * deltz), -(sig / gamma) * sinh(gamma * deltz))
            f22 = complex(cosh(gamma * deltz), (sig / gamma) * sinh(gamma * deltz))
            f12 = complex(0, -(kaa / gamma) * sinh(gamma * deltz))
            f21 = complex(0, +(kaa / gamma) * sinh(gamma * deltz))

            ftx = np.array([[f11, f12], [f21, f22]])
            f1 = np.dot(f1, ftx)

        return f1

    def undeformed_fbg(
        self,
        resolution: float,
        min_bandwidth: float,
        max_bandwidth: float,
        initial_refractive_index: float,
        mean_change_refractive_index: float,
        fringe_visibility: float,
        original_wavelengths: list,
    ):
        assert len(original_wavelengths) == self.fbg_count

        # TODO: move this in constructor, or even to class variables
        self.initial_refractive_index = initial_refractive_index
        self.mean_change_refractive_index = mean_change_refractive_index
        self.fringe_visibility = fringe_visibility

        original_fbg_periods = [
            wavelen / (2.0 * initial_refractive_index)
            for wavelen in original_wavelengths
        ]
        # Empty Original Reflec spectrum
        o_reflect = dict(
            wavelength=list(),
            reflec=list(),
        )

        # Cycle all the FBG sensors, but using only the periods of this cycle
        for period in original_fbg_periods:
            for wl in np.arange(min_bandwidth, max_bandwidth, resolution):
                f1 = self.transfer_matrix(period, wavelen=wl)
                PO = f1[0, 0]
                NO = f1[1, 0]
                reflectivity = abs(NO / PO) ** 2

                o_reflect["wavelength"].append(wl)
                o_reflect["reflec"].append(reflectivity)

        return o_reflect
