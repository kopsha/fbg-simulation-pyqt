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
        print("undeformed reflected signal", underformed_data["reflec"][:5])

    deformed_data = simu.deformed_fbg(
        strain_type=StrainTypes.NON_UNIFORM,
        stress_type=StressTypes.INCLUDED,
    )
    print("deformed reflected signal", deformed_data["reflec"][:5])

    summary_data = simu.compute_fbg_shifts_and_widths(
        strain_type=StrainTypes.NON_UNIFORM,
        stress_type=StressTypes.INCLUDED,
    )
    print(f"{summary_data=}")

    print("finished")
