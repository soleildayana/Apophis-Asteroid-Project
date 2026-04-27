from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import numpy as np
import plotly.graph_objects as go
import pymcel as pc


@dataclass(frozen=True)
class OrbitElements:
    """Keplerian orbital elements for a two-body problem.

    Parameters
    ----------
    a : float
        Semi-major axis (in *distance_unit*).
    e : float
        Eccentricity (0 <= e < 1, elliptical only).
    i : float
        Inclination (in *angle_unit*).
    raan : float
        Right ascension of the ascending node Ω (in *angle_unit*).
    argp : float
        Argument of periapsis ω (in *angle_unit*).
    epoch : float
        Reference epoch (in *time_unit*). M0 / nu0 is evaluated at this epoch.
    mu : float
        Gravitational parameter of the central body (distance_unit³/time_unit²).
    M0 : float, optional
        Mean anomaly at *epoch* (in *angle_unit*). Exactly one of M0/nu0 must be given.
    nu0 : float, optional
        True anomaly at *epoch* (in *angle_unit*). Exactly one of M0/nu0 must be given.
    angle_unit : str
        ``'deg'`` (default) or ``'rad'``.
    distance_unit : str
        Label only – e.g. ``'au'``.
    time_unit : str
        Label only – e.g. ``'day'``.
    frame : str
        Reference frame label (default ``'ecliptic'``).
    """

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
        """Convert *value* from the element's angle_unit to radians."""
        return np.deg2rad(value) if self.angle_unit == "deg" else value

    def period(self) -> float:
        """Orbital period in the element's time_unit."""
        return 2 * np.pi * np.sqrt(self.a**3 / self.mu)


def sample_times(t_start: float, t_end: float, n: int) -> np.ndarray:
    """Return *n* evenly-spaced times between *t_start* and *t_end*."""
    if n < 2:
        raise ValueError("n must be >= 2")
    return np.linspace(t_start, t_end, n)


# ──────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ──────────────────────────────────────────────────────────────────────────────

def _rotation_perifocal_to_ecliptic(raan: float, i: float, argp: float) -> np.ndarray:
    """3×3 rotation matrix from the perifocal frame to the ecliptic frame."""
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


def _orbit_ticks(
    positions: np.ndarray, n_ticks: int = 90
) -> tuple:
    """Picket-fence tick segments from orbit positions down to the ecliptic (z=0).

    Returns ``(xs, ys, zs)`` lists with ``None`` separators, ready for
    a ``go.Scatter3d`` with ``mode='lines'``.
    """
    n = len(positions)
    step = max(1, n // n_ticks)
    pts = positions[::step]
    xs, ys, zs = [], [], []
    for p in pts:
        xs += [float(p[0]), float(p[0]), None]
        ys += [float(p[1]), float(p[1]), None]
        zs += [float(p[2]), 0.0, None]
    return xs, ys, zs


# ──────────────────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────────────────

def orbital_plane_normal(i: float, raan: float, angle_unit: str = "deg") -> np.ndarray:
    """Unit normal vector of the orbital plane in the ecliptic frame.

    Parameters
    ----------
    i : float
        Inclination (in *angle_unit*).
    raan : float
        Right ascension of the ascending node Ω (in *angle_unit*).
    angle_unit : str
        ``'deg'`` or ``'rad'``.
    """
    if angle_unit not in {"deg", "rad"}:
        raise ValueError("angle_unit must be 'deg' or 'rad'.")
    if angle_unit == "deg":
        i = np.deg2rad(i)
        raan = np.deg2rad(raan)
    n_hat = np.array([np.sin(i) * np.sin(raan), -np.sin(i) * np.cos(raan), np.cos(i)])
    return n_hat / np.linalg.norm(n_hat)


def propagate_two_body(elements: OrbitElements, times: np.ndarray) -> np.ndarray:
    """Propagate a Keplerian two-body orbit.

    Parameters
    ----------
    elements : OrbitElements
    times : array_like
        Evaluation times (in the element's *time_unit*).

    Returns
    -------
    positions : ndarray, shape (N, 3)
        Ecliptic-frame Cartesian positions.
    """
    times = np.asarray(times, dtype=float)

    i_rad = elements.to_radians(elements.i)
    raan_rad = elements.to_radians(elements.raan)
    argp_rad = elements.to_radians(elements.argp)

    if elements.M0 is not None:
        M0 = elements.to_radians(elements.M0)
    else:
        nu0 = elements.to_radians(elements.nu0)
        M0 = _true_to_mean_anomaly(nu0, elements.e)

    n_motion = np.sqrt(elements.mu / elements.a**3)
    M_arr = M0 + n_motion * (times - elements.epoch)

    positions_pf = np.zeros((times.size, 3), dtype=float)
    for idx, m in enumerate(M_arr):
        E, _, _ = pc.kepler_newton(float(m), elements.e, G0=float(m), delta=1e-12)
        positions_pf[idx, 0] = elements.a * (np.cos(E) - elements.e)
        positions_pf[idx, 1] = elements.a * np.sqrt(1.0 - elements.e**2) * np.sin(E)

    rotation = _rotation_perifocal_to_ecliptic(raan_rad, i_rad, argp_rad)
    return (rotation @ positions_pf.T).T


def make_orbit_viewer(
    body_name: str,
    body_positions: np.ndarray,
    body_elements: OrbitElements,
    body_color: str = "#e0e0e0",
    body_epoch_pos: Optional[np.ndarray] = None,
    reference_bodies: Optional[List[dict]] = None,
    epoch_label: str = "",
    title: str = "Orbit Viewer",
    n_ticks: int = 90,
    show_wedge: bool = True,
) -> go.Figure:
    """Build a JPL-VOP-style interactive 3D orbit viewer.

    Parameters
    ----------
    body_name : str
        Name of the main body (e.g. ``'99942 Apophis'``).
    body_positions : ndarray, shape (N, 3)
        Ecliptic Cartesian positions for one complete orbit of the main body.
    body_elements : OrbitElements
        Elements used to derive the node line, orbital plane, and wedge.
    body_color : str
        Hex color for the main body's orbit and tick marks.
    body_epoch_pos : ndarray, shape (3,), optional
        Current (epoch) position of the main body; displayed as a white dot.
    reference_bodies : list of dict, optional
        Each dict must have ``'name'`` (str), ``'positions'`` (N×3 array),
        ``'color'`` (str), and optionally ``'epoch_pos'`` (shape-(3,) array).
    epoch_label : str
        Text appended to the info annotation, e.g. a date string.
    title : str
        Figure title.
    n_ticks : int
        Number of picket-fence tick marks drawn per orbit.
    show_wedge : bool
        Whether to render the shaded wedge between the ecliptic and orbital plane.

    Returns
    -------
    fig : plotly.graph_objects.Figure
    """
    reference_bodies = reference_bodies or []
    body_positions = np.asarray(body_positions)

    all_pos = [body_positions] + [np.asarray(b["positions"]) for b in reference_bodies]
    max_r = max(float(np.max(np.linalg.norm(p, axis=1))) for p in all_pos)
    extent = 1.3 * max_r

    i_rad = body_elements.to_radians(body_elements.i)
    raan_rad = body_elements.to_radians(body_elements.raan)
    n_hat = orbital_plane_normal(body_elements.i, body_elements.raan, angle_unit=body_elements.angle_unit)
    node_hat = np.array([np.cos(raan_rad), np.sin(raan_rad), 0.0])
    node_hat /= np.linalg.norm(node_hat)
    k_hat = np.array([0.0, 0.0, 1.0])
    p_hat = np.cross(k_hat, node_hat)
    p_hat /= np.linalg.norm(p_hat)
    p_orb = np.cos(i_rad) * p_hat + np.sin(i_rad) * k_hat

    fig = go.Figure()

    # ── Ecliptic grid ──────────────────────────────────────────────────────────
    n_grid = 14
    for gv in np.linspace(-extent, extent, n_grid):
        for xs, ys in [
            ([gv, gv], [-extent, extent]),
            ([-extent, extent], [gv, gv]),
        ]:
            fig.add_trace(
                go.Scatter3d(
                    x=xs, y=ys, z=[0.0, 0.0],
                    mode="lines",
                    line=dict(color="#2a2a2a", width=1),
                    showlegend=False,
                    hoverinfo="skip",
                )
            )

    # ── Node line (yellow) ─────────────────────────────────────────────────────
    nl = np.vstack([-extent * node_hat, extent * node_hat])
    fig.add_trace(
        go.Scatter3d(
            x=nl[:, 0], y=nl[:, 1], z=nl[:, 2],
            mode="lines",
            line=dict(color="#ffd700", width=4),
            name="Línea de nodos",
        )
    )

    # ── Orbital-plane normal (green) ───────────────────────────────────────────
    nm_end = 0.6 * extent * n_hat
    fig.add_trace(
        go.Scatter3d(
            x=[0.0, nm_end[0]], y=[0.0, nm_end[1]], z=[0.0, nm_end[2]],
            mode="lines",
            line=dict(color="#7CFC00", width=4),
            name="Normal plano orbital",
        )
    )

    # ── Wedge between ecliptic and orbital plane ───────────────────────────────
    if show_wedge and abs(body_elements.i) > 0.01:
        ring_a = np.linspace(0, 2 * np.pi, 120)
        blend = np.linspace(0, 1, 15)
        Xw = np.zeros((blend.size, ring_a.size))
        Yw = np.zeros_like(Xw)
        Zw = np.zeros_like(Xw)
        for j, phi in enumerate(ring_a):
            cphi, sphi = np.cos(phi), np.sin(phi)
            p_ecl = extent * (cphi * node_hat + sphi * p_hat)
            p_obj = extent * (cphi * node_hat + sphi * p_orb)
            for ib, beta in enumerate(blend):
                Xw[ib, j], Yw[ib, j], Zw[ib, j] = (1 - beta) * p_ecl + beta * p_obj
        fig.add_trace(
            go.Surface(
                x=Xw, y=Yw, z=Zw,
                opacity=0.18,
                showscale=False,
                colorscale=[[0, "#ffcf66"], [1, "#ffcf66"]],
                name="Banda entre planos",
            )
        )

    # ── Reference body orbits ──────────────────────────────────────────────────
    for body in reference_bodies:
        bpos = np.asarray(body["positions"])
        bcolor = body.get("color", "#808080")
        bname = body["name"]

        fig.add_trace(
            go.Scatter3d(
                x=bpos[:, 0], y=bpos[:, 1], z=bpos[:, 2],
                mode="lines",
                name=bname,
                line=dict(color=bcolor, width=4),
            )
        )
        tx, ty, tz = _orbit_ticks(bpos, n_ticks=n_ticks)
        fig.add_trace(
            go.Scatter3d(
                x=tx, y=ty, z=tz,
                mode="lines",
                line=dict(color=bcolor, width=1),
                showlegend=False,
                hoverinfo="skip",
            )
        )
        ep = body.get("epoch_pos")
        if ep is not None:
            ep = np.asarray(ep).ravel()
            fig.add_trace(
                go.Scatter3d(
                    x=[ep[0]], y=[ep[1]], z=[ep[2]],
                    mode="markers+text",
                    text=[bname],
                    textposition="top center",
                    textfont=dict(color=bcolor, size=11),
                    marker=dict(size=5, color="white"),
                    showlegend=False,
                )
            )

    # ── Main body orbit ────────────────────────────────────────────────────────
    fig.add_trace(
        go.Scatter3d(
            x=body_positions[:, 0],
            y=body_positions[:, 1],
            z=body_positions[:, 2],
            mode="lines",
            name=body_name,
            line=dict(color=body_color, width=5),
        )
    )
    tx, ty, tz = _orbit_ticks(body_positions, n_ticks=n_ticks)
    fig.add_trace(
        go.Scatter3d(
            x=tx, y=ty, z=tz,
            mode="lines",
            line=dict(color=body_color, width=1),
            showlegend=False,
            hoverinfo="skip",
        )
    )
    if body_epoch_pos is not None:
        ep = np.asarray(body_epoch_pos).ravel()
        fig.add_trace(
            go.Scatter3d(
                x=[ep[0]], y=[ep[1]], z=[ep[2]],
                mode="markers+text",
                text=[body_name],
                textposition="top center",
                textfont=dict(color=body_color, size=11),
                marker=dict(size=5, color="white"),
                showlegend=False,
            )
        )

    # ── Sun marker ────────────────────────────────────────────────────────────
    fig.add_trace(
        go.Scatter3d(
            x=[0], y=[0], z=[0],
            mode="markers+text",
            text=["Sol"],
            textposition="top center",
            textfont=dict(color="#ffd700", size=12),
            marker=dict(size=8, color="#fff14f"),
            name="Sol",
        )
    )

    # ── Info annotation ───────────────────────────────────────────────────────
    info_lines = [body_name]
    if body_epoch_pos is not None:
        ep_arr = np.asarray(body_epoch_pos).ravel()
        sun_dist = float(np.linalg.norm(ep_arr))
        info_lines.append(f"Distancia al Sol: {sun_dist:.3f} au")
        for b in reference_bodies:
            if b["name"].lower() in ("tierra", "earth"):
                ep_earth = b.get("epoch_pos")
                if ep_earth is not None:
                    d_earth = float(np.linalg.norm(ep_arr - np.asarray(ep_earth).ravel()))
                    info_lines.append(f"Distancia a la Tierra: {d_earth:.3f} au")
                break
    if epoch_label:
        info_lines.append(epoch_label)

    fig.add_annotation(
        text="<br>".join(info_lines),
        xref="paper", yref="paper",
        x=0.01, y=0.04,
        showarrow=False,
        font=dict(color="white", size=11, family="monospace"),
        align="left",
        bgcolor="rgba(0,0,0,0)",
    )

    # ── Layout ────────────────────────────────────────────────────────────────
    dist_unit = body_elements.distance_unit
    ax = dict(backgroundcolor="black", gridcolor="#1a1a1a", color="white", showticklabels=True)
    fig.update_layout(
        title=dict(text=title, font=dict(color="white")),
        paper_bgcolor="black",
        plot_bgcolor="black",
        scene=dict(
            xaxis=dict(title=f"X [{dist_unit}]", **ax),
            yaxis=dict(title=f"Y [{dist_unit}]", **ax),
            zaxis=dict(title=f"Z [{dist_unit}]", **ax),
            bgcolor="black",
            aspectmode="data",
        ),
        legend=dict(bgcolor="rgba(0,0,0,0.3)", font=dict(color="white")),
        font=dict(color="white"),
    )
    return fig


__all__ = [
    "OrbitElements",
    "sample_times",
    "propagate_two_body",
    "orbital_plane_normal",
    "make_orbit_viewer",
]
