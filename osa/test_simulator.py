"""
Testing the every function of the new simulator against the old one
The new implmentation is just a refactoring to a more pythonic state of the
old code, the computations should stay exactly the same.
"""

from simulator import OsaSimulator, SiUnits
from old_simulator import OSASimulation
import numpy as np


def assert_fbg_equals(left: dict, right: dict):
    assert left.keys() == right.keys()
    for fb_left, fb_right in zip(left.values(), right.values()):
        assert isinstance(fb_left, dict)
        assert isinstance(fb_right, dict)

        assert fb_left.keys() == fb_right.keys()
        for values_left, values_right in zip(fb_left.values(), fb_right.values()):
            assert (values_left == values_right).all()


def test_valid_data_from_file():
    """Compare FBG data with old simulator"""
    # params
    fbg_count = 5
    fbg_length = 9
    tolerance = 0.01
    fbg_positions = list()
    units = SiUnits.MILLIMETERS

    ref_sim = OSASimulation(
        filename="sample/tut-export.txt",
        NumberFBG=fbg_count,
        FBGLength=fbg_length,
        Tolerance=tolerance,
        SkipRow=8,
        FBGPosition=fbg_positions,
        InputUnits=int(units == SiUnits.MILLIMETERS),  # 0 -> meters, 1 -> mm
    )
    ref_data = ref_sim.FBGArray

    simu = OsaSimulator(
        fbg_count=fbg_count,
        fbg_length=fbg_length,
        tolerance=tolerance,
        fbg_positions=fbg_positions,
    )
    data = simu.from_file("sample/tut-export.txt", units=units)
    assert_fbg_equals(ref_data, data)


def test_undeformed_fbg():
    # initial params
    fbg_count = 5
    fbg_length = 9
    tolerance = 0.01
    fbg_positions = list()
    units = SiUnits.MILLIMETERS

    ref_sim = OSASimulation(
        filename="sample/tut-export.txt",
        NumberFBG=fbg_count,
        FBGLength=fbg_length,
        Tolerance=tolerance,
        SkipRow=8,
        FBGPosition=fbg_positions,
        InputUnits=int(units == SiUnits.MILLIMETERS),  # 0 -> meters, 1 -> mm
    )

    simu = OsaSimulator(
        fbg_count=fbg_count,
        fbg_length=fbg_length,
        tolerance=tolerance,
        fbg_positions=fbg_positions,
    )
    simu.from_file("sample/tut-export.txt", units=units)

    # undeformed fbg params
    resolution = np.float64(0.05)
    min_bandwidth = np.float64(1500.0)
    max_bandwidth = np.float64(1600.0)
    initial_refractive_index = np.float64(1.46)
    mean_change_refractive_index = np.float64("4.5E-4")
    fringe_visibility = np.float64(1)
    original_wavelengths = [
        1500.0,
        1525.0,
        1550.0,
        1575.0,
        1600.0,
    ]

    ref_sim.UndeformedFBG(
        SimulationResolution=resolution,
        MinBandWidth=min_bandwidth,
        MaxBandWidth=max_bandwidth,
        InitialRefractiveIndex=initial_refractive_index,
        MeanChangeRefractiveIndex=mean_change_refractive_index,
        FringeVisibility=fringe_visibility,
        DirectionalRefractiveP11=0,
        DirectionalRefractiveP12=0,
        PoissonsCoefficient=0,
        FBGOriginalWavel=original_wavelengths,
    )
    ref_data = ref_sim.OReflect

    data = simu.undeformed_fbg(
        resolution=resolution,
        min_bandwidth=min_bandwidth,
        max_bandwidth=max_bandwidth,
        initial_refractive_index=initial_refractive_index,
        original_wavelengths=original_wavelengths,
        mean_change_refractive_index=mean_change_refractive_index,
        fringe_visibility=fringe_visibility,
    )

    assert data["wavelength"] == ref_data["wavelength"]
    assert data["reflec"] == ref_data["reflec"]
