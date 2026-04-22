import numpy as np

from orbit_viewer import OrbitElements, make_orbit_viewer, propagate_two_body, sample_times


MU_SUN_AU_DAY2 = 2.9591220828559093e-4


def main() -> None:
    times = sample_times(-365.0, 365.0, 900)

    apophis = OrbitElements(
        a=0.9223803173917017,
        e=0.1911663355386932,
        i=3.340958441017069,
        raan=203.8996515621043,
        argp=126.6728325163065,
        nu0=25.0,
        epoch=0.0,
        mu=MU_SUN_AU_DAY2,
        angle_unit="deg",
        distance_unit="au",
        time_unit="day",
    )

    earth = OrbitElements(
        a=1.00000011,
        e=0.01671022,
        i=0.00005,
        raan=-11.26064,
        argp=114.20783,
        M0=100.46435,
        epoch=0.0,
        mu=MU_SUN_AU_DAY2,
        angle_unit="deg",
        distance_unit="au",
        time_unit="day",
    )

    apophis_positions = propagate_two_body(apophis, times)
    earth_positions = propagate_two_body(earth, times)

    fig = make_orbit_viewer(
        apophis_positions,
        apophis,
        reference_orbits={"Earth": earth_positions},
        title="Orbit Viewer - 99942 Apophis",
    )
    fig.write_html("orbit_viewer_demo.html", include_plotlyjs="cdn")
    print("Archivo generado: orbit_viewer_demo.html")


if __name__ == "__main__":
    main()
