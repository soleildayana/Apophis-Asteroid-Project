from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import plotly.graph_objects as go
import pymcel as pc


@dataclass(frozen=True)
class OrbitElements:
    a: float
    e: float
    i: float
    raan: float
    argp: float
    epoch: float
    mu: float
    M0: Optional[float] = None
    nu0: Optional[float] = None
    angle_unit: str = "deg"
    distance_unit: str = "au"
    time_unit: str = "day"
    frame: str = "ecliptic"

    def __post_init__(self) -> None:
        if not (0 <= self.e < 1):
            raise ValueError("Only elliptical orbits are supported (0 <= e < 1).")
        if (self.M0 is None) == (self.nu0 is None):
            raise ValueError("Provide exactly one of M0 or nu0.")
        if self.angle_unit not in {"deg", "rad"}:
            raise ValueError("angle_unit must be 'deg' or 'rad'.")

    def to_radians(self, value: float) -> float:
        return np.deg2rad(value) if self.angle_unit == "deg" else value


def sample_times(t_start: float, t_end: float, n: int) -> np.ndarray:
    if n < 2:
        raise ValueError("n must be >= 2")
    return np.linspace(t_start, t_end, n)


def _rotation_perifocal_to_ecliptic(raan: float, i: float, argp: float) -> np.ndarray:
    cO, sO = np.cos(raan), np.sin(raan)
    ci, si = np.cos(i), np.sin(i)
    cw, sw = np.cos(argp), np.sin(argp)
    return np.array(
        [
            [cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
            [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
            [sw * si, cw * si, ci],
        ]
    )


def _true_to_mean_anomaly(nu: float, e: float) -> float:
    E = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(nu / 2), np.sqrt(1 + e) * np.cos(nu / 2))
    return E - e * np.sin(E)


def _eccentric_to_true_anomaly(E: float, e: float) -> float:
    return 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))


def orbital_plane_normal(i: float, raan: float, angle_unit: str = "deg") -> np.ndarray:
    if angle_unit == "deg":
        i = np.deg2rad(i)
        raan = np.deg2rad(raan)
    if angle_unit != "deg" and angle_unit != "rad":
        raise ValueError("angle_unit must be 'deg' or 'rad'.")
    n_hat = np.array([np.sin(i) * np.sin(raan), -np.sin(i) * np.cos(raan), np.cos(i)])
    return n_hat / np.linalg.norm(n_hat)


def propagate_two_body(elements: OrbitElements, times: np.ndarray) -> np.ndarray:
    times = np.asarray(times, dtype=float)

    i = elements.to_radians(elements.i)
    raan = elements.to_radians(elements.raan)
    argp = elements.to_radians(elements.argp)

    if elements.M0 is not None:
        M0 = elements.to_radians(elements.M0)
    else:
        nu0 = elements.to_radians(elements.nu0 if elements.nu0 is not None else 0.0)
        M0 = _true_to_mean_anomaly(nu0, elements.e)

    n = np.sqrt(elements.mu / elements.a**3)
    M = M0 + n * (times - elements.epoch)

    positions_perifocal = np.zeros((times.size, 3), dtype=float)
    for idx, m in enumerate(M):
        E, _, _ = pc.kepler_newton(float(m), elements.e, G0=float(m), delta=1e-12)
        x = elements.a * (np.cos(E) - elements.e)
        y = elements.a * np.sqrt(1 - elements.e**2) * np.sin(E)
        positions_perifocal[idx] = [x, y, 0.0]

    rotation = _rotation_perifocal_to_ecliptic(raan, i, argp)
    return (rotation @ positions_perifocal.T).T


def make_orbit_viewer(
    positions: np.ndarray,
    elements: OrbitElements,
    reference_orbits: Optional[dict[str, np.ndarray]] = None,
    title: str = "Orbit Viewer",
) -> go.Figure:
    reference_orbits = reference_orbits or {}

    positions = np.asarray(positions)
    max_radius = float(np.max(np.linalg.norm(positions, axis=1)))
    for ref_positions in reference_orbits.values():
        max_radius = max(max_radius, float(np.max(np.linalg.norm(ref_positions, axis=1))))
    extent = 1.25 * max_radius

    i = elements.to_radians(elements.i)
    raan = elements.to_radians(elements.raan)
    n_hat = orbital_plane_normal(elements.i, elements.raan, angle_unit=elements.angle_unit)
    node_hat = np.array([np.cos(raan), np.sin(raan), 0.0])
    node_hat /= np.linalg.norm(node_hat)
    k_hat = np.array([0.0, 0.0, 1.0])
    p_hat = np.cross(k_hat, node_hat)
    p_hat /= np.linalg.norm(p_hat)
    p_orb = np.cos(i) * p_hat + np.sin(i) * k_hat

    fig = go.Figure()

    fig.add_trace(
        go.Surface(
            x=np.array([[-extent, extent], [-extent, extent]]),
            y=np.array([[-extent, -extent], [extent, extent]]),
            z=np.array([[0.0, 0.0], [0.0, 0.0]]),
            opacity=0.13,
            showscale=False,
            colorscale=[[0, "#66a3ff"], [1, "#66a3ff"]],
            name="Eclíptica",
        )
    )

    coeff = np.array([[-1, 1], [-1, 1]], dtype=float) * extent
    orbital_plane = np.zeros((2, 2, 3))
    for r in range(2):
        for c in range(2):
            orbital_plane[r, c] = coeff[r, c] * node_hat + coeff[c, r] * p_orb
    fig.add_trace(
        go.Surface(
            x=orbital_plane[:, :, 0],
            y=orbital_plane[:, :, 1],
            z=orbital_plane[:, :, 2],
            opacity=0.12,
            showscale=False,
            colorscale=[[0, "#d291ff"], [1, "#d291ff"]],
            name="Plano orbital",
        )
    )

    ring_angles = np.linspace(0, 2 * np.pi, 80)
    blend = np.linspace(0, 1, 20)
    X = np.zeros((blend.size, ring_angles.size))
    Y = np.zeros_like(X)
    Z = np.zeros_like(X)
    for j, phi in enumerate(ring_angles):
        cphi, sphi = np.cos(phi), np.sin(phi)
        p_ecl = extent * (cphi * node_hat + sphi * p_hat)
        p_obj = extent * (cphi * node_hat + sphi * p_orb)
        for i_b, beta in enumerate(blend):
            p_mid = (1 - beta) * p_ecl + beta * p_obj
            X[i_b, j], Y[i_b, j], Z[i_b, j] = p_mid
    fig.add_trace(
        go.Surface(
            x=X,
            y=Y,
            z=Z,
            opacity=0.3,
            showscale=False,
            colorscale=[[0, "#ffcf66"], [1, "#ffcf66"]],
            name="Banda entre planos",
        )
    )

    fig.add_trace(
        go.Scatter3d(
            x=positions[:, 0],
            y=positions[:, 1],
            z=positions[:, 2],
            mode="lines",
            name="Órbita objeto",
            line=dict(color="#e6e6e6", width=6),
        )
    )

    for orbit_name, orbit_positions in reference_orbits.items():
        fig.add_trace(
            go.Scatter3d(
                x=orbit_positions[:, 0],
                y=orbit_positions[:, 1],
                z=orbit_positions[:, 2],
                mode="lines",
                name=f"Órbita {orbit_name}",
                line=dict(width=5),
            )
        )

    node_line = np.vstack([-extent * node_hat, extent * node_hat])
    fig.add_trace(
        go.Scatter3d(
            x=node_line[:, 0],
            y=node_line[:, 1],
            z=node_line[:, 2],
            mode="lines",
            line=dict(color="#ffd60a", width=7),
            name="Línea de nodos",
        )
    )

    normal_line = np.vstack([np.zeros(3), extent * n_hat])
    fig.add_trace(
        go.Scatter3d(
            x=normal_line[:, 0],
            y=normal_line[:, 1],
            z=normal_line[:, 2],
            mode="lines+markers",
            line=dict(color="#7CFC00", width=6),
            marker=dict(size=3),
            name="Normal plano orbital",
        )
    )

    fig.add_trace(
        go.Scatter3d(
            x=[0], y=[0], z=[0], mode="markers+text", text=["Sun"], textposition="top center", name="Sol",
            marker=dict(size=7, color="#fff14f")
        )
    )

    axis_title = f"Distancia [{elements.distance_unit}]"
    fig.update_layout(
        title=title,
        paper_bgcolor="black",
        plot_bgcolor="black",
        scene=dict(
            xaxis=dict(title=axis_title, backgroundcolor="black", gridcolor="#333", color="white"),
            yaxis=dict(title=axis_title, backgroundcolor="black", gridcolor="#333", color="white"),
            zaxis=dict(title=axis_title, backgroundcolor="black", gridcolor="#333", color="white"),
            bgcolor="black",
            aspectmode="data",
        ),
        legend=dict(bgcolor="rgba(0,0,0,0.2)", font=dict(color="white")),
    )

    return fig


__all__ = [
    "OrbitElements",
    "sample_times",
    "propagate_two_body",
    "orbital_plane_normal",
    "make_orbit_viewer",
]
