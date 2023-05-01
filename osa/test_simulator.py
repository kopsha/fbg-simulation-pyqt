from simulator import OsaSimulator


def test_load_file():
    simu = OsaSimulator(fbg_count=5, fbg_length=9, tolerance=0.01, fbg_positions=[])
    simu.from_file("sample/tut-export.txt")
    assert False
