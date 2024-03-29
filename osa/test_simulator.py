"""
Testing the every function of the new simulator against the old one
The new implmentation is just a refactoring to a more pythonic state of the
old code, the computations should stay exactly the same.
"""
import pytest
import numpy as np

from old_simulator import OSASimulation
from simulator import OsaSimulator, SiUnits, StrainTypes, StressTypes


@pytest.fixture
def init_params():
    params = dict(
        fbg_count=5,
        fbg_length=9,
        tolerance=0.01,
        fbg_positions=list(),
        initial_refractive_index=1.46,
        directional_refractive_p11=0,
        directional_refractive_p12=0,
        poissons_coefficient=0,
        resolution=0.05,
        min_bandwidth=1500.0,
        max_bandwidth=1600.0,
        mean_change_refractive_index=4.5e-4,
        fringe_visibility=np.float64(1),
        emulate_temperature=293.15,  # 373.15,
        original_wavelengths=[1500.0, 1525.0, 1550.0, 1575.0, 1600.0],
        ambient_temperature=293.15,
        thermo_optic=8.3e-6,
        fiber_expansion_coefficient=10e-6,  # Internet says 0.5E-6 to 1E-6
        host_expansion_coefficient=10e-6,
        youngs_mod=75e9,
    )
    return params


def assert_fbg_equals(left: dict, right: dict):
    assert left.keys() == right.keys()
    for fb_left, fb_right in zip(left.values(), right.values()):
        assert isinstance(fb_left, dict)
        assert isinstance(fb_right, dict)

        assert fb_left.keys() == fb_right.keys()
        for values_left, values_right in zip(fb_left.values(), fb_right.values()):
            assert (values_left == values_right).all()


def test_valid_data_from_file(init_params):
    """Compare FBG data with old simulator"""
    units = SiUnits.MILLIMETERS
    ref_sim = OSASimulation(
        filename="sample/tut-export-limited.txt",
        NumberFBG=init_params["fbg_count"],
        FBGLength=init_params["fbg_length"],
        Tolerance=init_params["tolerance"],
        SkipRow=8,
        FBGPosition=init_params["fbg_positions"],
        InputUnits=units,
    )
    ref_data = ref_sim.FBGArray

    simu = OsaSimulator(**init_params)
    data = simu.from_file("sample/tut-export-limited.txt", units=units)
    assert_fbg_equals(ref_data, data)


def test_undeformed_fbg(init_params):
    units = SiUnits.MILLIMETERS
    ## Prepare reference simulation
    ref_sim = OSASimulation(
        filename="sample/tut-export-limited.txt",
        NumberFBG=init_params["fbg_count"],
        FBGLength=init_params["fbg_length"],
        Tolerance=init_params["tolerance"],
        SkipRow=8,
        FBGPosition=init_params["fbg_positions"],
        InputUnits=units,  # 0 -> meters, 1 -> mm
    )
    ref_sim.UndeformedFBG(
        SimulationResolution=init_params["resolution"],
        MinBandWidth=init_params["min_bandwidth"],
        MaxBandWidth=init_params["max_bandwidth"],
        InitialRefractiveIndex=init_params["initial_refractive_index"],
        MeanChangeRefractiveIndex=init_params["mean_change_refractive_index"],
        FringeVisibility=init_params["fringe_visibility"],
        DirectionalRefractiveP11=init_params["directional_refractive_p11"],
        DirectionalRefractiveP12=init_params["directional_refractive_p12"],
        PoissonsCoefficient=init_params["poissons_coefficient"],
        FBGOriginalWavel=init_params["original_wavelengths"],
    )
    ref_data = ref_sim.OReflect

    ## Prepare simulation to be tested
    simu = OsaSimulator(**init_params)
    simu.from_file("sample/tut-export-limited.txt", units=units)
    data = simu.undeformed_fbg()

    ## Compare results
    assert data["wavelength"] == ref_data["wavelength"]
    assert data["reflec"] == ref_data["reflec"]


def test_deformed_fbg(init_params):
    # initial params
    units = SiUnits.MILLIMETERS

    ## Prepare reference simulation
    ref_sim = OSASimulation(
        filename="sample/tut-export-limited.txt",
        NumberFBG=init_params["fbg_count"],
        FBGLength=init_params["fbg_length"],
        Tolerance=init_params["tolerance"],
        SkipRow=8,
        FBGPosition=init_params["fbg_positions"],
        InputUnits=units,  # 0 -> meters, 1 -> mm
    )
    ref_sim.DeformedFBG(
        SimulationResolution=init_params["resolution"],
        MinBandWidth=init_params["min_bandwidth"],
        MaxBandWidth=init_params["max_bandwidth"],
        AmbientTemperature=init_params["ambient_temperature"],
        InitialRefractiveIndex=init_params["initial_refractive_index"],
        MeanChangeRefractiveIndex=init_params["mean_change_refractive_index"],
        FringeVisibility=init_params["fringe_visibility"],
        DirectionalRefractiveP11=init_params["directional_refractive_p11"],
        DirectionalRefractiveP12=init_params["directional_refractive_p12"],
        YoungsModule=init_params["youngs_mod"],
        PoissonsCoefficient=init_params["poissons_coefficient"],
        ThermoOptic=init_params["thermo_optic"],
        EmulateTemperature=init_params["emulate_temperature"],
        FiberThermalExpansionCoefficient=init_params["fiber_expansion_coefficient"],
        HostThermalExpansionCoefficient=init_params["host_expansion_coefficient"],
        FBGOriginalWavel=init_params["original_wavelengths"],
        StrainType=StrainTypes.NON_UNIFORM,
        StressType=StressTypes.INCLUDED,
    )
    ref_data = ref_sim.DReflect
    ref_Y_data = ref_sim.YReflect
    ref_Z_data = ref_sim.ZReflect

    ## Prepare simulation to be tested
    simu = OsaSimulator(**init_params)
    simu.from_file("sample/tut-export-limited.txt", units=units)
    data = simu.deformed_fbg(
        strain_type=StrainTypes.NON_UNIFORM,
        stress_type=StressTypes.INCLUDED,
    )

    assert data["wavelength"] == ref_data["wavelength"]
    assert np.array_equal(data["reflec"], ref_data["reflec"])

    assert np.array_equal(data["Y_split"], ref_Y_data["reflec"])
    assert np.array_equal(data["Z_split"], ref_Z_data["reflec"])


def test_output_sum(init_params):
    # initial params
    units = SiUnits.MILLIMETERS

    ## Prepare reference simulation
    ref_sim = OSASimulation(
        filename="sample/tut-export-limited.txt",
        NumberFBG=init_params["fbg_count"],
        FBGLength=init_params["fbg_length"],
        Tolerance=init_params["tolerance"],
        SkipRow=8,
        FBGPosition=init_params["fbg_positions"],
        InputUnits=units,  # 0 -> meters, 1 -> mm
    )
    ref_sim.FBGOutputSum(
        AmbientTemperature=init_params["ambient_temperature"],
        InitialRefractiveIndex=init_params["initial_refractive_index"],
        FringeVisibility=init_params["fringe_visibility"],
        DirectionalRefractiveP11=init_params["directional_refractive_p11"],
        DirectionalRefractiveP12=init_params["directional_refractive_p12"],
        YoungsModule=init_params["youngs_mod"],
        PoissonsCoefficient=init_params["poissons_coefficient"],
        ThermoOptic=init_params["thermo_optic"],
        StrainType=StrainTypes.NON_UNIFORM,
        StressType=StressTypes.INCLUDED,
        EmulateTemperature=init_params["emulate_temperature"],
        FiberThermalExpansionCoefficient=init_params["fiber_expansion_coefficient"],
        HostThermalExpansionCoefficient=init_params["host_expansion_coefficient"],
        FBGOriginalWavel=init_params["original_wavelengths"],
    )
    ref_data = ref_sim.FBGOutSum
    ref_fbg_stats = ref_sim._FBGmaxmin

    ## Prepare simulation to be tested
    simu = OsaSimulator(**init_params)
    simu.from_file("sample/tut-export-limited.txt", units=units)
    data = simu.compute_fbg_shifts_and_widths(
        strain_type=StrainTypes.NON_UNIFORM,
        stress_type=StressTypes.INCLUDED,
    )
    fbg_stats = simu._fbg_stats

    ## Assert starting stats are equal too
    for fkey in fbg_stats.keys():
        stats = fbg_stats[fkey]
        ref_stats = ref_fbg_stats[fkey]
        for measure in stats.keys():
            assert np.isclose(stats[measure], ref_stats[measure])

    assert data.keys() == ref_data.keys()
    for fkey in data.keys():
        wave_shift = data[fkey]["wave_shift"]
        assert len(ref_data[fkey]["WaveShift"]) == 1
        ref_wave_shift = ref_data[fkey]["WaveShift"][0]
        assert np.isclose(wave_shift, ref_wave_shift)

        wave_shift = data[fkey]["wave_width"]
        assert len(ref_data[fkey]["WaveWidth"]) == 1
        ref_wave_shift = ref_data[fkey]["WaveWidth"][0]
        assert np.isclose(wave_shift, ref_wave_shift)
