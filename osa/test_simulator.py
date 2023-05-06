from simulator import OsaSimulator, SiUnits
from old_simulator import OSASimulation


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
