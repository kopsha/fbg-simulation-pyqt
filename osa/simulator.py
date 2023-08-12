"""
Simulation of reflected FBG spectrum using coupled-mode theory.
TODO: review with [Ben Frey](https://github.com/benfrey)
"""
import numpy as np
from cmath import pi, sqrt, cosh, sinh
from enum import IntEnum, auto


class SiUnits(IntEnum):
    METERS = auto()
    MILLIMETERS = auto()


class StrainTypes(IntEnum):
    NONE = 0
    UNIFORM = 1
    NON_UNIFORM = 2


class OsaSimulator:
    def __init__(
        self,
        fbg_count: int,
        fbg_length: float,
        tolerance: float,
        fbg_positions: list,
        original_wavelengths: list,
        initial_refractive_index: float,
        directional_refractive_p11: float,
        directional_refractive_p12: float,
        poissons_coefficient: float,
        emulate_temperature: float,
    ):
        """
        Prepares an OSA simulation, with following parameters:

        :param int fbg_count: Number of FBG per optical fibre
        :param float fbg_length: Length of the Gratting
        :param float tolerance: Tolerance in the FBG length
        :param list positions: Position of each FBG from the begginng of the Path
        :param float initial_refractive_index: Initial effective refractive index (neff)
        :param float directional_refractive_p11: Pockel’s normal photoelastic constant
        :param float directional_refractive_p12: Pockel’s shear photoelastic constant
        :param float poissons_coefficient: Poisson's Coefficient of fiber
        :param float emulate_temperature: Theoretical emulated temperature of fiber

        NB: Everything in this class will be expressed in (mm)
        """
        self.data = None
        self.fbg_count = int(fbg_count)
        self.fbg_length = np.float64(fbg_length)
        self.tolerance = np.float64(tolerance)
        self.fbg_positions = np.array(fbg_positions, dtype=np.float64)

        assert len(original_wavelengths) == self.fbg_count
        self.original_wavelengths = original_wavelengths
        self.APFBG = original_wavelengths[: self.fbg_count] / (
            2.0 * initial_refractive_index
        )

        self.initial_refractive_index = initial_refractive_index
        self.directional_refractive_p11 = directional_refractive_p11
        self.directional_refractive_p12 = directional_refractive_p12
        self.poissons_coefficient = poissons_coefficient
        self.emulate_temperature = emulate_temperature

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
        mean_change_refractive_index: float,
        fringe_visibility: float,
    ):
        # TODO: move this in constructor, or even to class variables
        self.mean_change_refractive_index = mean_change_refractive_index
        self.fringe_visibility = fringe_visibility

        original_fbg_periods = [
            wavelen / (2.0 * self.initial_refractive_index)
            for wavelen in self.original_wavelengths
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

    def deformed_fbg(
        self,
        strain_type: StrainTypes,
        ambient_temperature: float,
        thermo_optic: float,
        fiber_expansion_coefficient: float,
        host_expansion_coefficient: float,
    ):
        """
        :param IntEnum strain_type: 0 for none, 1 for uniform, 2 for non-uniform
        :param float ambient_temperature: Base in which our reference temperature is set.
        :param float thermo_optic: Thermo optic coefficient of fiber
        :param float fiber_expansion_coefficient: Thermal expansion coefficient of fiber material
        :param float host_expansion_coefficient : Thermal expansion coefficient of host material
        """
        # Calculate photoelastic coef from directional coefs.
        self.photo_elastic_param = (self.initial_refractive_index**2 / 2) * (
            self.directional_refractive_p12
            - self.poissons_coefficient
            * (self.directional_refractive_p11 + self.directional_refractive_p12)
        )

        # Determine if we need to emulate a theoretical temperature value
        if self.emulate_temperature != -1.0:
            for i in range(self.fbg_count):
                key = f"FBG{i+1}"
                self.fbg[key]["T"][:] = self.emulate_temperature

        # Empty Reflec spectrum- Two waves (individual contributions to account for any transverse stress)
        self.YReflect = {}  # Y wave
        self.YReflect["wavelength"] = []
        self.YReflect["reflec"] = []

        self.ZReflect = {}  # Z wave
        self.ZReflect["wavelength"] = []
        self.ZReflect["reflec"] = []

        # Empty Reflec spectrum- Composite wave of Y and Z contributions
        self.DReflect = {}  # Composite wave
        self.DReflect["wavelength"] = []
        self.DReflect["reflec"] = []

        # Cycle all the FBG sensors
        for i in range(self.fbg_count):
            key = f"FBG{i+1}"
            # Sections the gratting is divided-- Transfer Matrix
            M = len(self.fbg[key]["x"])
            # FBG increment size (nm)
            deltz = (self.fbg_length * (10.0**6)) / M

            ## STRAIN
            fbg_period = []
            if strain_type == StrainTypes.NONE:  # no longitudinal strain
                fbg_period = [self.APFBG[i] for _ in range(M)]
            elif strain_type == StrainTypes.UNIFORM:  # uniform longitudinal strain
                # Strain average over FBG region
                strain_mean = np.mean(self.fbg[key]["LE11"])
                # Temperature average over FBG region
                temp_mean = np.mean(self.fbg[key]["T"])

                # Calculate wavelength at uniform strain and temperature
                thermo_dynamic_part = (
                    fiber_expansion_coefficient
                    + (1 - self.photo_elastic_param)
                    * (host_expansion_coefficient - fiber_expansion_coefficient)
                    + thermo_optic
                )

                new_wave_lengths = self.original_wavelengths[i] * (
                    1
                    + (1 - self.photo_elastic_param) * strain_mean
                    + thermo_dynamic_part * (temp_mean - ambient_temperature)
                )

                fbg_period = [
                    new_wave_lengths / (2.0 * self.initial_refractive_index)
                    for _ in range(M)
                ]

            # elif strain_type == StrainTypes.NON_UNIFORM:
            #     ## Case of "non-uniform longitudinal strain" build the grating period changed ---
            #     for j in np.arange(0, M):
            #         TempnewWavelength = self.FBGOriginalWavel[i] * (
            #             1
            #             + (1 - self.PhotoElasticParam)
            #             * self.FBGArray["FBG" + str(i + 1)]["LE11"][j]
            #             + (
            #                 self.FiberThermalExpansionCoefficient
            #                 + (1 - self.PhotoElasticParam)
            #                 * (
            #                     self.HostThermalExpansionCoefficient
            #                     - self.FiberThermalExpansionCoefficient
            #                 )
            #                 + self.ThermoOptic
            #             )
            #             * (
            #                 self.FBGArray["FBG" + str(i + 1)]["T"][j]
            #                 - self.AmbientTemperature
            #             )
            #         )  # weavelength at nonuniform strain and temperature
            #         # print(self.FBGArray['FBG'+str(i+1)]['T'][0]-self.AmbientTemperature)
            #         FBGperiod.append(
            #             TempnewWavelength / (2.0 * self.InitialRefractiveIndex)
            #         )  # Grating period
            else:
                raise ValueError(f"{strain_type} is not a valid strain_type")

        return []
