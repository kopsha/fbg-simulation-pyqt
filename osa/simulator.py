"""
Simulation of reflected FBG spectrum using coupled-mode theory.
TODO: review with [Ben Frey](https://github.com/benfrey)
"""
import numpy as np
from cmath import pi, sqrt, cosh, sinh
from enum import IntEnum


class SiUnits(IntEnum):
    METERS = 0
    MILLIMETERS = 1


class StrainTypes(IntEnum):
    NONE = 0
    UNIFORM = 1
    NON_UNIFORM = 2


class StressTypes(IntEnum):
    NONE = 0
    INCLUDED = 1


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
        resolution: float,
        min_bandwidth: float,
        max_bandwidth: float,
        mean_change_refractive_index: float,
        fringe_visibility: float,
        ambient_temperature: float,
        thermo_optic: float,
        fiber_expansion_coefficient: float,
        host_expansion_coefficient: float,
        youngs_mod: float,
    ):
        """
        Prepares an OSA simulation with the specified parameters.

        Parameters
        ----------
        fbg_count : int
            Number of FBG per optical fibre.
        fbg_length : float
            Length of the Grating.
        tolerance : float
            Tolerance in the FBG length.
        positions : list
            Position of each FBG from the beginning of the Path.
        initial_refractive_index : float
            Initial effective refractive index (neff).
        mean_change_refractive_index : float
            Mean induced change in the refractive index (dneff).
        directional_refractive_p11 : float
            Pockel’s normal photoelastic constant.
        directional_refractive_p12 : float
            Pockel’s shear photoelastic constant.
        poissons_coefficient : float
            Poisson's Coefficient of fiber.
        emulate_temperature : float
            Theoretical emulated temperature of fiber.
        resolution : float
            Simulation resolution - Wavelength increment.
        min_bandwidth : float
            Light Minimum Bandwidth.
        max_bandwidth : float
            Light Maximum Bandwidth.
        fringe_visibility : float
            Fringe Visibility (FV).
        ambient_temperature : float
            Base in which our reference temperature is set.
        thermo_optic : float
            Thermo optic coefficient of fiber.
        fiber_expansion_coefficient : float
            Thermal expansion coefficient of fiber material.
        host_expansion_coefficient : float
            Thermal expansion coefficient of host material.
        youngs_mod : float
            Young's modulus of fiber.

        Note
        ----
        Everything in this class will be expressed in (mm).
        """

        self.fbg_count = int(fbg_count)
        self.fbg_length = np.float64(fbg_length)
        self.tolerance = np.float64(tolerance)
        self.fbg_positions = np.array(fbg_positions, dtype=np.float64)
        self.initial_refractive_index = np.float64(initial_refractive_index)
        self.mean_change_refractive_index = np.float64(mean_change_refractive_index)
        self.directional_refractive_p11 = np.float64(directional_refractive_p11)
        self.directional_refractive_p12 = np.float64(directional_refractive_p12)
        self.poissons_coefficient = np.float64(poissons_coefficient)
        self.emulate_temperature = np.float64(emulate_temperature)
        self.resolution = np.float64(resolution)
        self.min_bandwidth = np.float64(min_bandwidth)
        self.max_bandwidth = np.float64(max_bandwidth)
        self.fringe_visibility = np.float64(fringe_visibility)
        self.ambient_temperature = np.float64(ambient_temperature)
        self.thermo_optic = np.float64(thermo_optic)
        self.fiber_expansion_coefficient = np.float64(fiber_expansion_coefficient)
        self.host_expansion_coefficient = np.float64(host_expansion_coefficient)
        self.youngs_mod = np.float64(youngs_mod)

        assert len(original_wavelengths) == self.fbg_count
        self.original_wavelengths = np.array(original_wavelengths, dtype=np.float64)
        self.grating_periods = self.original_wavelengths[: self.fbg_count] / (
            2.0 * initial_refractive_index
        )
        self.original_fbg_periods = [
            wavelen / (2.0 * self.initial_refractive_index)
            for wavelen in self.original_wavelengths
        ]

        self.directional_refractive_p11 = directional_refractive_p11
        self.directional_refractive_p12 = directional_refractive_p12
        self.poissons_coefficient = poissons_coefficient
        self.emulate_temperature = emulate_temperature
        self.resolution = resolution
        self.min_bandwidth = min_bandwidth
        self.max_bandwidth = max_bandwidth
        self.fringe_visibility = fringe_visibility

        self.ambient_temperature = ambient_temperature
        self.thermo_optic = thermo_optic
        self.fiber_expansion_coefficient = fiber_expansion_coefficient
        self.host_expansion_coefficient = host_expansion_coefficient
        self.youngs_mod = youngs_mod

        self._compute_photo_elastic_coefficient()

    def _compute_photo_elastic_coefficient(self):
        """
        Calculate photoelastic coefficient from directional coefficients.
        """
        self.photo_elastic_param = (self.initial_refractive_index**2 / 2) * (
            self.directional_refractive_p12
            - self.poissons_coefficient
            * (self.directional_refractive_p11 + self.directional_refractive_p12)
        )

    def from_file(self, filepath: str, units=SiUnits.MILLIMETERS) -> dict:
        """
        Loads FBG data from a text file.

        The file is expected to have the following fields:
            - 'x', 'LE11', 'LE22', 'LE33', 'S11', 'S22', 'S33', 'T'

        Parameters
        ----------
        filepath : str
            Full path to the data file.
        units : SiUnits, optional
            Units of the data (default is SiUnits.MILLIMETERS).
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
        dtypes = dict(names=list(fields.keys()), formats=[v["format"] for v in fields.values()])
        raw_data = np.genfromtxt(filepath, dtypes, comments="%")

        fbg = dict()
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

    def sigma(self, period: float, wavelen: float, dneff: float) -> float:
        """
        Compute the sigma value based on period, wavelength, and change in effective refractive index.

        Parameters
        ----------
        period : float
            FBG grating period.
        wavelen : float
            Wavelength under consideration.
        dneff : float
            Change in effective refractive index.

        Returns
        -------
        float
            Computed sigma value.
        """
        refractive_term = self.initial_refractive_index + dneff
        inverse_wavelength_difference = (1.0 / wavelen) - (1.0 / (2.0 * refractive_term * period))
        mean_refraction_effect = 2.0 * pi * self.mean_change_refractive_index / wavelen

        return 2.0 * pi * refractive_term * inverse_wavelength_difference + mean_refraction_effect

    def kaa(self, wavelen: float) -> float:
        """
        Compute the kaa value based on the provided wavelength.

        Parameters
        ----------
        wavelen : float
            Wavelength under consideration.

        Returns
        -------
        float
            Computed kaa value.
        """
        return pi * self.fringe_visibility * self.mean_change_refractive_index / wavelen

    def transfer_matrix(
        self, count: int, wavelen: float, use_period, use_dneff: list = []
    ) -> np.ndarray:
        """
        Calculate the transfer matrix for the FBG based on various parameters.

        Parameters
        ----------
        count : int
            Number of divisions to consider for the FBG length.
        wavelen : float
            Wavelength under consideration.
        use_period : float or list
            FBG grating period. Can be a float (for uniform period) or a list (for varying period).
        use_dneff : list, optional
            List of changes in effective refractive index.

        Returns
        -------
        np.ndarray
            Calculated transfer matrix.
        """
        trx_mat = np.identity(2)
        delta_z = self.fbg_length * (10**6) / count

        for z in range(count):
            period = use_period if isinstance(use_period, float) else use_period[z]
            dneff = use_dneff[z] if len(use_dneff) else 0.0
            sig = self.sigma(period, wavelen=wavelen, dneff=dneff)
            kaa_value = self.kaa(wavelen=wavelen)
            gamma = sqrt(kaa_value**2 - sig**2)

            f11 = complex(cosh(gamma * delta_z), -(sig / gamma) * sinh(gamma * delta_z))
            f22 = complex(cosh(gamma * delta_z), (sig / gamma) * sinh(gamma * delta_z))
            f12 = complex(0, -(kaa_value / gamma) * sinh(gamma * delta_z))
            f21 = complex(0, +(kaa_value / gamma) * sinh(gamma * delta_z))

            section_transfer_matrix = np.array([[f11, f12], [f21, f22]])
            trx_mat = np.dot(trx_mat, section_transfer_matrix)

        return trx_mat

    def undeformed_fbg(self) -> dict:
        """
        Calculate the undeformed (original) reflection spectrum of the FBG.

        Returns
        -------
        dict
            Dictionary containing the undeformed reflection spectrum with keys 'wavelength' and 'reflec'.
        """
        reflection_spectrum = {"wavelength": [], "reflec": []}
        # ASSUMPTION: all columns have the same length
        M = len(self.fbg["FBG1"]["x"])
        # Cycle through all the FBG sensors using only the periods for this cycle
        for period in self.original_fbg_periods:
            # Number of sections the grating is divided into for the Transfer Matrix
            wavelengths = np.arange(self.min_bandwidth, self.max_bandwidth, self.resolution)

            for wl in wavelengths:
                trx_mat = self.transfer_matrix(count=M, use_period=period, wavelen=wl)
                PO, NO = trx_mat[0, 0], trx_mat[1, 0]
                reflectivity = abs(NO / PO) ** 2

                reflection_spectrum["wavelength"].append(wl)
                reflection_spectrum["reflec"].append(reflectivity)

        return reflection_spectrum

    def deformed_fbg(self, strain_type: StrainTypes, stress_type: StressTypes) -> dict:
        """
        Calculate the deformed reflection spectrum of the FBG considering strain and stress effects.

        Parameters
        ----------
        strain_type : StrainTypes
            Type of strain to consider: NONE, UNIFORM, or NON_UNIFORM.
        stress_type : StressTypes
            Type of stress to consider: NONE or INCLUDED.

        Returns
        -------
        dict
            Dictionary containing the deformed reflection spectrum with keys 'wavelength' and 'reflec'.
        """
        # Apply a theoretical temperature value if specified
        if self.emulate_temperature != -1.0:
            for i in range(self.fbg_count):
                self.fbg[f"FBG{i+1}"]["T"][:] = self.emulate_temperature

        y_reflection = {"wavelength": [], "reflec": []}
        z_reflection = {"wavelength": [], "reflec": []}
        combined_reflection = {"wavelength": [], "reflec": []}

        # Compute the thermo-dynamic part
        thermo_dynamic_effect = (
            self.fiber_expansion_coefficient
            + (1 - self.photo_elastic_param)
            * (self.host_expansion_coefficient - self.fiber_expansion_coefficient)
            + self.thermo_optic
        )

        # Iterate over all the FBG sensors
        for i in range(self.fbg_count):
            sensor_data = self.fbg[f"FBG{i+1}"]
            M = len(sensor_data["x"])

            # Determine the FBG grating period based on the strain type
            if strain_type == StrainTypes.NONE:
                fbg_period = [self.grating_periods[i] for _ in range(M)]
            elif strain_type == StrainTypes.UNIFORM:
                strain_avg = np.mean(sensor_data["LE11"])
                temp_avg = np.mean(sensor_data["T"])
                new_wavelength = self.original_wavelengths[i] * (
                    1
                    + (1 - self.photo_elastic_param) * strain_avg
                    + thermo_dynamic_effect * (temp_avg - self.ambient_temperature)
                )
                fbg_period = [
                    new_wavelength / (2.0 * self.initial_refractive_index) for _ in range(M)
                ]
            elif strain_type == StrainTypes.NON_UNIFORM:
                fbg_period = [
                    self.original_wavelengths[i]
                    * (
                        1
                        + (1 - self.photo_elastic_param) * sensor_data["LE11"][j]
                        + thermo_dynamic_effect * (sensor_data["T"][j] - self.ambient_temperature)
                    )
                    / (2.0 * self.initial_refractive_index)
                    for j in range(M)
                ]
            else:
                raise ValueError(f"Invalid strain_type: {strain_type}.")

            # Determine the change in effective refractive index based on the stress type
            self.dneff_y = np.zeros(M)
            self.dneff_z = np.zeros(M)
            if stress_type == StressTypes.INCLUDED:
                coef = -(self.initial_refractive_index**3.0) / (2 * self.youngs_mod)
                factor1 = (
                    self.directional_refractive_p11
                    - 2 * self.poissons_coefficient * self.directional_refractive_p12
                )
                factor2 = (
                    (1 - self.poissons_coefficient) * self.directional_refractive_p12
                    - self.poissons_coefficient * self.directional_refractive_p11
                )
                self.dneff_y = coef * (
                    factor1 * sensor_data["S22"]
                    + factor2 * (sensor_data["S33"] + sensor_data["S11"])
                )
                self.dneff_z = coef * (
                    factor1 * sensor_data["S33"]
                    + factor2 * (sensor_data["S22"] + sensor_data["S11"])
                )
            elif stress_type != StressTypes.NONE:
                raise ValueError(f"Invalid stress_type: {stress_type}.")

            # Simulate for Y and Z waves
            wavelengths = np.arange(self.min_bandwidth, self.max_bandwidth, self.resolution)
            for wl in wavelengths:
                trx_y = self.transfer_matrix(
                    count=M, use_period=fbg_period, wavelen=wl, use_dneff=self.dneff_y
                )
                PO_y, NO_y = trx_y[0, 0], trx_y[1, 0]
                reflectivity_y = abs(NO_y / PO_y) ** 2
                y_reflection["wavelength"].append(wl)
                y_reflection["reflec"].append(reflectivity_y)

                trx_z = self.transfer_matrix(
                    count=M, use_period=fbg_period, wavelen=wl, use_dneff=self.dneff_z
                )
                PO_z, NO_z = trx_z[0, 0], trx_z[1, 0]
                reflectivity_z = abs(NO_z / PO_z) ** 2
                z_reflection["wavelength"].append(wl)
                z_reflection["reflec"].append(reflectivity_z)

        # Combine the Y and Z wave reflections
        combined_reflection["wavelength"] = y_reflection["wavelength"]
        combined_reflection["reflec"] = np.add(
            np.divide(y_reflection["reflec"], 2.0), np.divide(z_reflection["reflec"], 2.0)
        )

        return combined_reflection

    def _calculate_strain_stats(self):
        """
        Calculates statistics related to strain and stress for each FBG sensor.

        For each FBG sensor, the function computes the average, maximum, and minimum
        values for different strain and stress types (LE11, LE22, LE33, S11, S22, S33, T).

        Returns
        -------
        dict
            A dictionary with the computed statistics for each FBG sensor.
        """

        stat_types = ["AV", "Max", "Min"]
        measures = ["LE11", "LE22", "LE33", "S11", "S22", "S33", "T"]
        operations = {"AV": np.mean, "Max": np.max, "Min": np.min}

        stats_dict = {
            f"FBG{i+1}": {
                f"{stat}-{measure}": operations[stat](self.fbg[f"FBG{i+1}"][measure])
                if not (stat == "AV" and measure == "T" and self.emulate_temperature != -1.0)
                else self.emulate_temperature
                for measure in measures
                for stat in stat_types
            }
            for i in range(self.fbg_count)
        }

        return stats_dict

    def _calculate_wave_shift(
        self, strain_type: StrainTypes, fbg_stat: dict, original_wavelength: float
    ) -> float:
        """
        Calculate the wavelength shift for a given strain type.

        The method computes the wavelength shift based on the provided strain type. For `StrainTypes.NONE`,
        the shift is determined solely by the thermo-optic effect. For other strain types, the shift
        is influenced by various factors including photo-elastic parameters, fiber and host expansion
        coefficients, and thermo-optic effects.

        Parameters
        ----------
        strain_type : StrainTypes
            Type of strain to consider. Options are: NONE, UNIFORM, or NON_UNIFORM.
        fbg_stat : dict
            Dictionary containing statistics related to Fiber Bragg Gratings (FBG).
            Expected to contain keys like "AV-T" for average temperature and "AV-LE11" for average longitudinal strain.
        original_wavelength : float
            The original wavelength of the FBG before any strain or temperature effects.

        Returns
        -------
        float
            The calculated wavelength shift based on the provided parameters.

        Notes
        -----
        - The calculation is based on the assumption that strain and temperature changes are small
        enough to be considered linear.
        - The method uses class attributes such as `self.thermo_optic`, `self.photo_elastic_param`,
        `self.fiber_expansion_coefficient`, `self.host_expansion_coefficient`, and `self.ambient_temperature`
        for its calculations.
        """
        temperature_diff = fbg_stat["AV-T"] - self.ambient_temperature

        if strain_type == StrainTypes.NONE:
            return original_wavelength * self.thermo_optic * temperature_diff
        else:
            common_term = (1 - self.photo_elastic_param) * fbg_stat["AV-LE11"] + (
                self.fiber_expansion_coefficient
                + (1 - self.photo_elastic_param)
                * (self.host_expansion_coefficient - self.fiber_expansion_coefficient)
                + self.thermo_optic
            ) * temperature_diff
            return original_wavelength * common_term

    def _calculate_grating_periods(
        self, strain_type: StrainTypes, fbg_stat: dict, original_wavelength: float
    ) -> tuple:
        """
        Calculate the grating periods for the provided strain type.

        For `StrainTypes.NONE` and `StrainTypes.UNIFORM`, both the maximum and
        minimum grating periods are the same. For `StrainTypes.NON_UNIFORM`, the
        grating periods are determined separately for the maximum and minimum strains.

        Parameters
        ----------
        strain_type : StrainTypes
            Type of strain to consider. Options are: NONE, UNIFORM, or NON_UNIFORM.
        fbg_stat : dict
            Dictionary containing statistics related to Fiber Bragg Gratings (FBG).
            Expected to contain keys like "AV-T" for average temperature, "AV-LE11" for average
            longitudinal strain, "Max-LE11" for maximum longitudinal strain, and "Min-LE11" for minimum
            longitudinal strain.
        original_wavelength : float
            The original wavelength of the FBG before any strain or temperature effects.

        Returns
        -------
        tuple
            A tuple containing two values: the maximum grating period and the minimum grating period.

        Notes
        -----
        The calculation is based on the assumption that strain and temperature
        changes are small enough to be considered linear.
        """
        temperature_diff = fbg_stat["AV-T"] - self.ambient_temperature

        common_term = (
            original_wavelength
            / (2 * self.initial_refractive_index)
            * (
                1
                + (1 - self.photo_elastic_param) * fbg_stat["AV-LE11"]
                + self.thermo_optic * temperature_diff
            )
        )

        if strain_type in [StrainTypes.NONE, StrainTypes.UNIFORM]:
            return common_term, common_term
        elif strain_type == StrainTypes.NON_UNIFORM:
            grating_period_max = (
                original_wavelength
                / (2 * self.initial_refractive_index)
                * (
                    1
                    + (1 - self.photo_elastic_param) * fbg_stat["Max-LE11"]
                    + self.thermo_optic * temperature_diff
                )
            )
            grating_period_min = (
                original_wavelength
                / (2 * self.initial_refractive_index)
                * (
                    1
                    + (1 - self.photo_elastic_param) * fbg_stat["Min-LE11"]
                    + self.thermo_optic * temperature_diff
                )
            )
            return grating_period_max, grating_period_min

    def _calculate_peak_width(
        self,
        strain_type: StrainTypes,
        stress_type: StressTypes,
        fbg_stat: dict,
        original_wavelength: float,
        grating_period_max: float,
        grating_period_min: float,
    ) -> float:
        """
        Computes the peak width of a Fiber Bragg Grating (FBG) based on strain and
        stress considerations. The width is influenced by non-uniform strain (if
        present) and by the transverse stress (if included).

        Parameters
        ----------
        strain_type : StrainTypes
            Type of strain to consider. Options are: NONE, UNIFORM, or NON_UNIFORM.
        stress_type : StressTypes
            Type of stress to consider. Options are: NONE or INCLUDED.
        fbg_stat : dict
            Dictionary containing statistics related to Fiber Bragg Gratings (FBG).
            Expected to contain keys like "AV-S33" and "AV-S22" for stress values
            in different directions.
        original_wavelength : float
            The original wavelength of the FBG before any strain or temperature effects.
        grating_period_max : float
            Maximum grating period calculated based on strain effects.
        grating_period_min : float
            Minimum grating period calculated based on strain effects.

        Returns
        -------
        float
            The calculated peak width for the FBG based on the provided strain and stress types.

        Notes
        -----
        The method considers the peak width contributions only if the corresponding
        strain or stress type is specified.
        """
        peak_width_1 = (
            2 * self.initial_refractive_index * (grating_period_max - grating_period_min)
            if strain_type == StrainTypes.NON_UNIFORM
            else 0
        )
        peak_width_2 = (
            (
                abs(fbg_stat["AV-S33"] - fbg_stat["AV-S22"])
                * (original_wavelength / (2 * self.initial_refractive_index))
                * (
                    (1 + self.poissons_coefficient) * self.directional_refractive_p12
                    - (1 + self.poissons_coefficient) * self.directional_refractive_p11
                )
                * self.initial_refractive_index**3
                / self.youngs_mod
            )
            if stress_type == StressTypes.INCLUDED
            else 0
        )

        return peak_width_1 + peak_width_2

    def compute_fbg_shifts_and_widths(
        self, strain_type: StrainTypes, stress_type: StressTypes
    ) -> dict:
        """
        Calculate the wavelength shift and peak width for a set of FBGs.

        This method computes the wavelength shift and peak width for each FBG in
        the set based on the given strain and stress types. The calculations are
        determined using helper methods for wave shift, grating periods, and peak
        width.

        Parameters
        ----------
        strain_type : StrainTypes
            Type of strain to consider. Options are: NONE, UNIFORM, or NON_UNIFORM.
        stress_type : StressTypes
            Type of stress to consider. Options are: NONE or INCLUDED.

        Returns
        -------
        dict
            A dictionary containing the calculated wavelength shift and peak width
            for each FBG. The keys are FBG identifiers (e.g., "FBG1", "FBG2", ...),
            and the values are dictionaries with keys "wave_shift" and "wave_width",
            each containing a list of corresponding values.

        Raises
        ------
        ValueError
            If the provided strain_type or stress_type is not valid.
        """
        if strain_type not in StrainTypes._value2member_map_:
            raise ValueError(f"{strain_type} is not a valid strain type.")

        if stress_type not in StressTypes._value2member_map_:
            raise ValueError(f"{strain_type} is not a valid stress type.")

        fbg_stats = self._calculate_strain_stats()
        output = {f"FBG{i+1}": {"wave_shift": [], "wave_width": []} for i in range(self.fbg_count)}

        for i in range(self.fbg_count):
            key = f"FBG{i+1}"
            fbg_stat = fbg_stats[key]
            original_wavelength = self.original_wavelengths[i]

            wavelength_shift = self._calculate_wave_shift(
                strain_type, fbg_stat, original_wavelength
            )
            grating_period_max, grating_period_min = self._calculate_grating_periods(
                strain_type, fbg_stat, original_wavelength
            )
            peak_width_total = self._calculate_peak_width(
                strain_type,
                stress_type,
                fbg_stat,
                self.original_wavelengths[i],
                grating_period_max,
                grating_period_min,
            )

            output[key]["wave_shift"].append(wavelength_shift)
            output[key]["wave_width"].append(peak_width_total)

        return output
