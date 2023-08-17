from .simulator import OsaSimulator, SiUnits, StrainTypes, StressTypes


def simulator_worker(params: dict):
    print("started")

    units = params.pop("units")
    include_underformed_signal = params.pop("has_reflected_signal")
    strain_type = params.pop("strain_type")
    stress_type = params.pop("stress_type")
    datafile = params.pop("filepath")

    print(params)

    simu = OsaSimulator(**params)
    simu.from_file(filepath=datafile, units=units)

    if include_underformed_signal:
        underformed_data = simu.undeformed_fbg()
        print("reflected signal", underformed_data["reflec"][:5])

    print("finished")
