#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
postprocess_cylinder_flow_field.py
==================================

Render per-VTP flow-field images for the 3D Eulerian fully compressible
cylinder case (Ma = 0.3, Re = 100).

The script reads ``Cylinder3DFluid_*.vtp`` frames produced by
``test_3d_eulerian_compressible_flow_around_cylinder_LG``, draws the cylinder as
an opaque solid obstacle, and writes one PNG per (VTP frame, field) pair under
``results/flow_field_frames/``.

Fields (``--fields``, default ``speed,mach,pressure,density``):
  speed     |u|, absolute
  mach      |u| / sqrt(gamma*p/rho), the true local Mach number
  pressure  p / p_inf
  density   rho / rho_inf
  vorticity |curl(u)|, only if the solver wrote a vorticity array

The default view is a 3D scatter on the fixed particle cloud plus the opaque
cylinder.  ``grid-slice`` remains available when a conventional XY cell field
is preferred.

Simulation time is read from each VTP's ``TimeValue`` FieldData scalar
(written by SPHinXsys) rather than reverse-engineered from the frame
index, so the script is robust against restart runs that reuse the same
output directory.

Example
-------
::

    # Render every Cylinder3DFluid_*.vtp in output/ to results/flow_field_frames/
    python postprocess_cylinder_flow_field.py \
        --data_dir output \
        --results_dir results \
        --config config.ini

    # Render a single mid-span XY grid field instead of the default 3D view
    python postprocess_cylinder_flow_field.py --mode grid-slice
"""

from __future__ import annotations

import argparse
import configparser
import math
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pyvista as pv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)


# ---------------------------------------------------------------------------
# Defaults / configuration parsing
# ---------------------------------------------------------------------------

DEFAULT_FRAMES_DIRNAME = "flow_field_frames"
DEFAULT_FLUID_PREFIX = "Cylinder3DFluid_"
DEFAULT_VELOCITY_KEY = "Velocity"
DEFAULT_PARTICLE_PROJECTED_AREA_FRACTION = 0.35


@dataclass
class CylinderCaseConfig:
    """Lightweight view of the case geometry / physical parameters."""

    DL: float
    DH: float
    DW: float
    dp: float
    cylinder_center_x: float
    cylinder_center_y: float
    cylinder_radius: float
    gamma: float
    rho_inf: float
    c_inf: float
    mach_inf: float
    re: float

    @property
    def D(self) -> float:
        return 2.0 * self.cylinder_radius

    @property
    def u_inf(self) -> float:
        """Freestream speed, derived exactly as the C++ config does."""
        return self.mach_inf * self.c_inf

    @property
    def p_inf(self) -> float:
        """Freestream pressure p = rho*c^2/gamma."""
        return self.rho_inf * self.c_inf * self.c_inf / self.gamma

    def to_dict(self) -> Dict[str, float]:
        """Return init-only fields; the derived properties are excluded so the
        dict can be fed straight back into ``__init__`` after a process hop
        (e.g. ProcessPoolExecutor worker).
        """
        return {
            "DL": self.DL,
            "DH": self.DH,
            "DW": self.DW,
            "dp": self.dp,
            "cylinder_center_x": self.cylinder_center_x,
            "cylinder_center_y": self.cylinder_center_y,
            "cylinder_radius": self.cylinder_radius,
            "gamma": self.gamma,
            "rho_inf": self.rho_inf,
            "c_inf": self.c_inf,
            "mach_inf": self.mach_inf,
            "re": self.re,
        }


def parse_float(value: str, name: str) -> float:
    """Parse a float with whitespace tolerance; raise with context on failure."""
    try:
        return float(str(value).strip())
    except (TypeError, ValueError) as exc:
        raise ValueError(f"failed to parse {name}={value!r}: {exc}") from exc


def load_case_config(config_path: Path) -> CylinderCaseConfig:
    """Load geometry + physical fields from the case's ``config.ini``."""
    if not config_path.exists():
        raise FileNotFoundError(f"config.ini not found: {config_path}")

    parser = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    parser.read(config_path, encoding="utf-8")

    def need(section: str, key: str) -> str:
        if section not in parser:
            raise KeyError(f"missing [{section}] in {config_path}")
        if key not in parser[section]:
            raise KeyError(f"missing key '{key}' in [{section}] of {config_path}")
        return parser[section][key]

    return CylinderCaseConfig(
        DL=parse_float(need("geometry", "dl"), "dl"),
        DH=parse_float(need("geometry", "dh"), "dh"),
        DW=parse_float(need("geometry", "dw"), "dw"),
        dp=parse_float(need("geometry", "global_resolution"), "global_resolution"),
        cylinder_center_x=parse_float(need("geometry", "cylinder_center_x"), "cylinder_center_x"),
        cylinder_center_y=parse_float(need("geometry", "cylinder_center_y"), "cylinder_center_y"),
        cylinder_radius=parse_float(need("geometry", "cylinder_radius"), "cylinder_radius"),
        gamma=parse_float(need("physical", "gamma"), "gamma"),
        rho_inf=parse_float(need("physical", "rho_inf"), "rho_inf"),
        c_inf=parse_float(need("physical", "c_inf"), "c_inf"),
        mach_inf=parse_float(need("physical", "mach_inf"), "mach_inf"),
        re=parse_float(need("physical", "re"), "re"),
    )


# ---------------------------------------------------------------------------
# Frame discovery / time decoding
# ---------------------------------------------------------------------------

def extract_frame_index(vtp_path: Path) -> int:
    """Extract the iteration / step counter encoded in the VTP filename.

    SPHinXsys uses the body-name + fixed-width integer pattern, e.g.
    ``Cylinder3DFluid_0010000040.vtp`` or
    ``Cylinder3DFluid_ite_0000000000.vtp``. The body name itself may
    contain digits (``Cylinder3``), so we anchor on the trailing numeric
    token that follows the optional ``ite_`` marker instead of using a
    bare ``re.search(r"(\\d+)", ...)`` (which would always match the
    first ``3`` in ``Cylinder3`` and yield frame_index = 3 for every
    file).
    """
    stem = vtp_path.stem
    match = re.search(r"(?:^|_)ite_(\d+)$", stem) or re.search(r"_(\d+)$", stem)
    if not match:
        raise ValueError(f"cannot parse frame index from {vtp_path.name}")
    return int(match.group(1))


def read_time_value(vtp_file: Path) -> Optional[float]:
    """Return the SPHinXsys ``TimeValue`` field stored in a VTP, if any.

    SPHinXsys encodes the simulation time as a top-level FieldData scalar
    on every VTP it writes, which is far more reliable than reverse-
    engineering the iteration count from the filename (which depends on
    the simulation's output cadence).
    """
    mesh = pv.read(vtp_file)
    try:
        if "TimeValue" in mesh.field_data:
            arr = np.asarray(mesh.field_data["TimeValue"]).reshape(-1)
            if arr.size:
                return float(arr[0])
    finally:
        del mesh
    return None


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

def load_fluid_scalar(
    vtp_file: Path,
    cfg: CylinderCaseConfig,
    field: str,
    velocity_key: str,
) -> Tuple[np.ndarray, np.ndarray, str, Optional[Tuple[float, float]]]:
    """Return ``(position, values, colorbar_label, default_range)`` for one frame.

    Compressible fields are reported non-dimensionally against the freestream so
    a plot is readable without knowing the absolute scale:
      speed     -> |u| (absolute, matches the historical behaviour)
      mach      -> |u| / sqrt(gamma*p/rho), the true local Mach number
      pressure  -> p / p_inf
      density   -> rho / rho_inf
      vorticity -> |curl(u)|, only if the solver wrote a vorticity array
    """
    mesh = pv.read(vtp_file)
    try:
        available = list(mesh.point_data.keys())
        position = np.asarray(mesh.points, dtype=np.float64)

        def need_array(key: str) -> np.ndarray:
            if key not in mesh.point_data:
                raise KeyError(
                    f"field '{key}' missing in {vtp_file.name}; available: {available}"
                )
            return np.asarray(mesh.point_data[key], dtype=np.float64)

        if field == "speed":
            velocity = need_array(velocity_key)
            values = np.linalg.norm(velocity[:, :3], axis=1)
            return position, values, "Velocity magnitude |u|", None

        if field == "mach":
            velocity = need_array(velocity_key)
            pressure = need_array("Pressure").reshape(-1)
            density = need_array("Density").reshape(-1)
            if np.any(density <= 0.0) or np.any(pressure <= 0.0):
                raise ValueError(
                    f"non-positive density or pressure in {vtp_file.name}; "
                    "the local sound speed is undefined"
                )
            sound_speed = np.sqrt(cfg.gamma * pressure / density)
            values = np.linalg.norm(velocity[:, :3], axis=1) / sound_speed
            return position, values, "Local Mach number", (0.0, max(0.6, 2.0 * cfg.mach_inf))

        if field == "pressure":
            pressure = need_array("Pressure").reshape(-1)
            values = pressure / cfg.p_inf
            # Frame-independent range so an animation does not flicker. For an
            # isentropic subsonic field p/p_inf stays within roughly
            # (1 +- gamma*Ma^2), widened for the stagnation/suction peaks.
            half_width = max(0.05, 2.0 * cfg.gamma * cfg.mach_inf * cfg.mach_inf)
            return position, values, r"$p / p_\infty$", (1.0 - half_width, 1.0 + half_width)

        if field == "density":
            density = need_array("Density").reshape(-1)
            values = density / cfg.rho_inf
            # Same reasoning; density varies by ~Ma^2 in a subsonic field.
            half_width = max(0.05, 2.0 * cfg.mach_inf * cfg.mach_inf)
            return position, values, r"$\rho / \rho_\infty$", (1.0 - half_width, 1.0 + half_width)

        if field == "vorticity":
            vorticity_key = next(
                (k for k in ("Vorticity", "VorticityVector", "VelocityCurl") if k in mesh.point_data),
                None,
            )
            if vorticity_key is None:
                raise KeyError(
                    f"no vorticity array in {vtp_file.name} (looked for Vorticity / "
                    f"VorticityVector / VelocityCurl); available: {available}. "
                    "Add the field to the solver's output list to plot it."
                )
            vorticity = need_array(vorticity_key)
            values = (
                np.linalg.norm(vorticity[:, :3], axis=1)
                if vorticity.ndim == 2 and vorticity.shape[1] >= 3
                else np.abs(vorticity.reshape(-1))
            )
            return position, values, r"$|\nabla \times u|$", None

        raise ValueError(f"unknown field {field!r}")
    finally:
        del mesh


def maybe_downsample(
    position: np.ndarray,
    values: np.ndarray,
    max_points: int,
    seed: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Uniform random downsample to keep the scatter under ``max_points``."""
    n = position.shape[0]
    if max_points <= 0 or n <= max_points:
        return position, values
    rng = np.random.default_rng(seed)
    sample_idx = np.sort(rng.choice(n, size=max_points, replace=False))
    return position[sample_idx], values[sample_idx]


# ---------------------------------------------------------------------------
# Plotting primitives
# ---------------------------------------------------------------------------

def _format_time_tag(time_value: float) -> str:
    """Filesystem-safe time tag, e.g. ``t00012p345678s``."""
    return f"{time_value:010.6f}".replace(".", "p") + "s"


def _scalar_color_range(
    values: np.ndarray,
    explicit: Optional[Tuple[float, float]],
    field_default: Optional[Tuple[float, float]],
) -> Tuple[float, float]:
    """Resolve the colour range: CLI override > field default > data extent."""
    if explicit is not None:
        vmin, vmax = explicit
        if vmin >= vmax:
            raise ValueError(f"invalid color range [{vmin}, {vmax}]")
        return float(vmin), float(vmax)
    if field_default is not None:
        return field_default
    # Falling through to the data extent means this field autoscales PER FRAME,
    # so an animation of it will flicker. Only 'vorticity' lands here, because
    # its magnitude has no a-priori freestream scale.
    if values.size == 0:
        return 0.0, 1.0
    vmin = float(np.nanmin(values))
    vmax = float(np.nanmax(values))
    if not math.isfinite(vmin) or not math.isfinite(vmax) or vmin >= vmax:
        # A uniform field would give a degenerate colour bar.
        return vmin - 0.5, vmin + 0.5
    return vmin, vmax


def _velocity_color_range(
    cfg: CylinderCaseConfig,
    explicit: Optional[Tuple[float, float]],
) -> Tuple[float, float]:
    if explicit is not None:
        vmin, vmax = explicit
        if vmin >= vmax:
            raise ValueError(f"invalid color range [{vmin}, {vmax}]")
        return float(vmin), float(vmax)
    if cfg.u_inf <= 0.0:
        raise ValueError(f"u_inf must be positive to derive default color range, got {cfg.u_inf}")
    # Default: 0 .. 1.6 * u_inf (Re=100 cylinder flow rarely exceeds ~1.5 * u_inf)
    return 0.0, 1.6 * float(cfg.u_inf)


def _draw_cylinder_solid_xy(ax, cfg: CylinderCaseConfig) -> None:
    """Draw the mid-span cylinder section as an opaque solid obstacle."""
    ax.add_patch(
        Circle(
            (cfg.cylinder_center_x, cfg.cylinder_center_y),
            cfg.cylinder_radius,
            facecolor="#1f2933",
            edgecolor="#080c10",
            linewidth=0.9,
            zorder=5,
        )
    )


def _draw_cylinder_solid_3d(ax, cfg: CylinderCaseConfig) -> None:
    """Draw the finite spanwise cylinder as an opaque shaded surface."""
    theta = np.linspace(0.0, 2.0 * np.pi, 96)
    z = np.linspace(0.0, cfg.DW, 24)
    theta_grid, z_grid = np.meshgrid(theta, z)
    x_grid = cfg.cylinder_center_x + cfg.cylinder_radius * np.cos(theta_grid)
    y_grid = cfg.cylinder_center_y + cfg.cylinder_radius * np.sin(theta_grid)
    ax.plot_surface(
        x_grid,
        y_grid,
        z_grid,
        color="#1f2933",
        linewidth=0.0,
        antialiased=True,
        shade=True,
        alpha=1.0,
        zorder=1,
    )
    radial = np.linspace(0.0, cfg.cylinder_radius, 32)
    radial_grid, cap_theta = np.meshgrid(radial, theta)
    cap_x = cfg.cylinder_center_x + radial_grid * np.cos(cap_theta)
    cap_y = cfg.cylinder_center_y + radial_grid * np.sin(cap_theta)
    for cap_z in (0.0, cfg.DW):
        ax.plot_surface(
            cap_x,
            cap_y,
            np.full_like(cap_x, cap_z),
            color="#1f2933",
            linewidth=0.0,
            antialiased=True,
            shade=True,
            alpha=1.0,
            zorder=1,
        )


def _occluded_by_cylinder(
    position: np.ndarray,
    cfg: CylinderCaseConfig,
    elev: float,
    azim: float,
) -> np.ndarray:
    """Return points hidden by the finite solid cylinder in the orthographic view.

    ``mplot3d`` sorts a scatter collection and a surface only as whole objects,
    not point-by-point. Cast an orthographic ray from each fluid point toward
    the camera so that flow behind the analytic cylinder cannot be painted on
    top of its opaque surface.
    """
    elev_rad = math.radians(elev)
    azim_rad = math.radians(azim)
    to_camera = np.array(
        [
            math.cos(elev_rad) * math.cos(azim_rad),
            math.cos(elev_rad) * math.sin(azim_rad),
            math.sin(elev_rad),
        ],
        dtype=np.float64,
    )
    rel_x = position[:, 0] - cfg.cylinder_center_x
    rel_y = position[:, 1] - cfg.cylinder_center_y
    rel_z = position[:, 2]
    dx, dy, dz = to_camera
    radius_squared = cfg.cylinder_radius * cfg.cylinder_radius

    # Side-surface intersections of p + t * to_camera, t > 0.
    a = dx * dx + dy * dy
    b = 2.0 * (rel_x * dx + rel_y * dy)
    c = rel_x * rel_x + rel_y * rel_y - radius_squared
    discriminant = b * b - 4.0 * a * c
    side_hit = np.zeros(position.shape[0], dtype=bool)
    if a > np.finfo(float).eps:
        valid = discriminant >= 0.0
        root = np.sqrt(np.maximum(discriminant, 0.0))
        for t in ((-b - root) / (2.0 * a), (-b + root) / (2.0 * a)):
            z_hit = rel_z + t * dz
            side_hit |= valid & (t > 1.0e-12) & (z_hit >= 0.0) & (z_hit <= cfg.DW)

    # End-cap intersections complete the finite-cylinder visibility test.
    cap_hit = np.zeros(position.shape[0], dtype=bool)
    if abs(dz) > np.finfo(float).eps:
        for cap_z in (0.0, cfg.DW):
            t = (cap_z - rel_z) / dz
            x_hit = rel_x + t * dx
            y_hit = rel_y + t * dy
            cap_hit |= (t > 1.0e-12) & (x_hit * x_hit + y_hit * y_hit <= radius_squared)

    return side_hit | cap_hit


def estimate_volume_marker_size(
    cfg: CylinderCaseConfig,
    *,
    figure_width_in: float,
    axes_width_fraction: float,
) -> float:
    """Estimate scatter area from one SPH particle's physical volume.

    A particle occupies ``V=dp^3``.  For plotting, replace that cube by the
    equal-volume sphere and use 35% of its projected disk area. The reduction
    retains a particle-cloud appearance without turning the full 3D depth
    stack into an opaque colour sheet. The orthographic view has a uniform
    physical-to-screen scale, estimated from the widest in-plane domain span
    and the allocated axes width. Matplotlib's scatter ``s`` is measured in
    points squared, so this returns the corresponding screen-space area.
    """
    particle_volume = cfg.dp ** 3
    particle_radius = (3.0 * particle_volume / (4.0 * math.pi)) ** (1.0 / 3.0)
    widest_in_plane_span = max(cfg.DL, cfg.DH)
    points_per_length = 72.0 * figure_width_in * axes_width_fraction / widest_in_plane_span
    return DEFAULT_PARTICLE_PROJECTED_AREA_FRACTION * math.pi * (
        particle_radius * points_per_length
    ) ** 2


def plot_3d_cloud(
    *,
    position: np.ndarray,
    speed: np.ndarray,
    cfg: CylinderCaseConfig,
    time_value: float,
    output_path: Path,
    marker_size: Optional[float],
    elev: float,
    azim: float,
    color_range: Tuple[float, float],
    color_label: str = "Velocity magnitude |u|",
) -> Dict[str, float]:
    """Render an opaque high-contrast 3D scalar cloud around a solid cylinder."""
    fig = plt.figure(figsize=(8.8, 5.6))
    ax = fig.add_subplot(1, 1, 1, projection="3d", computed_zorder=False)
    ax.set_position([0.01, 0.02, 0.84, 0.92])
    effective_marker_size = (
        marker_size
        if marker_size is not None
        else estimate_volume_marker_size(cfg, figure_width_in=8.8, axes_width_fraction=0.84)
    )

    _draw_cylinder_solid_3d(ax, cfg)
    hidden = _occluded_by_cylinder(position, cfg, elev, azim)
    visible_position = position[~hidden]
    visible_speed = speed[~hidden]
    sc = ax.scatter(
        visible_position[:, 0],
        visible_position[:, 1],
        visible_position[:, 2],
        c=visible_speed,
        cmap="turbo",
        marker="o",
        s=effective_marker_size,
        vmin=color_range[0],
        vmax=color_range[1],
        edgecolors="none",
        linewidths=0.0,
        depthshade=False,
        alpha=1.0,
        zorder=2,
    )

    # The analytic surface is deliberate: particle-wall points made the obstacle
    # translucent and visually noisy at this resolution. Automatic mplot3d
    # z-sorting is disabled, so only ray-visible points are painted over it.

    ax.set_title(
        f"3D Flow Field (t = {time_value:.4f}s, Ma = {cfg.mach_inf:.2f}, Re = {cfg.re:.0f})",
        pad=4,
    )
    ax.view_init(elev=elev, azim=azim)
    ax.set_proj_type("ortho")
    ax.set_box_aspect([cfg.DL, cfg.DH, cfg.DW])
    ax.set_axis_off()

    cbar = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.01, shrink=0.84)
    cbar.set_label(color_label)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    return {
        "n_total": int(position.shape[0]),
        "speed_max": float(np.nanmax(speed)) if speed.size else 0.0,
        "speed_mean": float(np.nanmean(speed)) if speed.size else 0.0,
    }


def plot_xy_slice(
    *,
    position: np.ndarray,
    speed: np.ndarray,
    cfg: CylinderCaseConfig,
    time_value: float,
    output_path: Path,
    slice_half_width: float,
    marker_size: Optional[float],
    color_range: Tuple[float, float],
    color_label: str = "Velocity magnitude |u|",
) -> Dict[str, float]:
    """Render an XY slice of the selected scalar at z = DW/2 within ±slice_half_width."""
    z_mid = 0.5 * cfg.DW
    mask = np.abs(position[:, 2] - z_mid) <= slice_half_width
    pos = position[mask]
    spd = speed[mask]

    fig, ax = plt.subplots(figsize=(8.8, 4.4))
    ax.set_position([0.06, 0.10, 0.84, 0.82])
    effective_marker_size = (
        marker_size
        if marker_size is not None
        else estimate_volume_marker_size(cfg, figure_width_in=8.8, axes_width_fraction=0.84)
    )

    sc = ax.scatter(
        pos[:, 0],
        pos[:, 1],
        c=spd,
        cmap="turbo",
        marker="o",
        s=effective_marker_size,
        vmin=color_range[0],
        vmax=color_range[1],
        edgecolors="none",
        linewidths=0.0,
    )

    _draw_cylinder_solid_xy(ax, cfg)

    ax.set_xlim(0.0, cfg.DL)
    ax.set_ylim(0.0, cfg.DH)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (m, streamwise)")
    ax.set_ylabel("y (m, vertical)")
    ax.set_title(
        f"XY Slice @ z = {z_mid:.3f} ± {slice_half_width:.3f}  |  "
        f"t = {time_value:.4f}s  |  Ma = {cfg.mach_inf:.2f}  |  Re = {cfg.re:.0f}"
    )
    ax.grid(True, alpha=0.25, linewidth=0.4)

    cbar = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.015)
    cbar.set_label(color_label)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    return {
        "n_total": int(position.shape[0]),
        "n_slice": int(pos.shape[0]),
        "speed_max": float(np.nanmax(spd)) if spd.size else 0.0,
        "speed_mean": float(np.nanmean(spd)) if spd.size else 0.0,
    }


def _fixed_xy_layer_grid(
    position: np.ndarray,
    values: np.ndarray,
    cfg: CylinderCaseConfig,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Map the fixed Eulerian mid-span particle layer to native XY cells.

    No interpolation is used: every particle in the closest z layer owns its
    Cartesian cell.  Missing cells remain masked, which exposes the cylinder
    interior until the analytic solid is drawn over it.
    """
    z_values = np.unique(position[:, 2])
    if z_values.size == 0:
        raise ValueError("cannot build a grid slice from an empty position array")
    z_layer = float(z_values[np.argmin(np.abs(z_values - 0.5 * cfg.DW))])
    z_tolerance = max(1.0e-12, 0.25 * cfg.dp)
    layer_mask = np.abs(position[:, 2] - z_layer) <= z_tolerance
    layer_position = position[layer_mask]
    layer_values = values[layer_mask]
    if layer_position.size == 0:
        raise ValueError("closest mid-span layer contains no fluid particles")

    x_values = np.unique(layer_position[:, 0])
    y_values = np.unique(layer_position[:, 1])
    grid = np.full((y_values.size, x_values.size), np.nan, dtype=float)
    x_index = np.searchsorted(x_values, layer_position[:, 0])
    y_index = np.searchsorted(y_values, layer_position[:, 1])
    grid[y_index, x_index] = layer_values
    return x_values, y_values, grid, z_layer


def plot_xy_grid_slice(
    *,
    position: np.ndarray,
    speed: np.ndarray,
    cfg: CylinderCaseConfig,
    time_value: float,
    output_path: Path,
    color_range: Tuple[float, float],
    color_label: str = "Velocity magnitude |u|",
) -> Dict[str, float]:
    """Render a conventional XY field from the fixed Eulerian particle grid."""
    x_values, y_values, grid, z_layer = _fixed_xy_layer_grid(position, speed, cfg)
    cmap = plt.colormaps["turbo"].copy()
    cmap.set_bad((1.0, 1.0, 1.0, 0.0))

    fig, ax = plt.subplots(figsize=(8.8, 4.4))
    ax.set_position([0.06, 0.10, 0.84, 0.82])
    sc = ax.pcolormesh(
        x_values,
        y_values,
        np.ma.masked_invalid(grid),
        shading="nearest",
        cmap=cmap,
        vmin=color_range[0],
        vmax=color_range[1],
        antialiased=False,
    )
    _draw_cylinder_solid_xy(ax, cfg)

    ax.set_xlim(0.0, cfg.DL)
    ax.set_ylim(0.0, cfg.DH)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (m, streamwise)")
    ax.set_ylabel("y (m, vertical)")
    ax.set_title(
        f"XY Grid Field @ z = {z_layer:.3f}  |  t = {time_value:.4f}s  |  "
        f"Ma = {cfg.mach_inf:.2f}  |  Re = {cfg.re:.0f}"
    )
    ax.grid(False)

    cbar = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.015)
    cbar.set_label(color_label)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    finite_values = grid[np.isfinite(grid)]
    return {
        "n_total": int(position.shape[0]),
        "n_slice": int(finite_values.size),
        "speed_max": float(np.nanmax(finite_values)) if finite_values.size else 0.0,
        "speed_mean": float(np.nanmean(finite_values)) if finite_values.size else 0.0,
    }


# ---------------------------------------------------------------------------
# Worker
# ---------------------------------------------------------------------------

def render_frame_task(task: Dict[str, object]) -> Dict[str, object]:
    """Render one frame in a worker process."""
    fluid_file = Path(str(task["fluid_file"]))
    output_path = Path(str(task["output_path"]))
    cfg_dict = dict(task["cfg_dict"])
    cfg = CylinderCaseConfig(**cfg_dict)

    field = str(task["field"])
    position, values, color_label, field_default = load_fluid_scalar(
        fluid_file, cfg, field, velocity_key=str(task["velocity_key"])
    )
    n_total = int(position.shape[0])
    mode = str(task["mode"])
    if mode == "grid-slice":
        position_d, values_d = position, values
    else:
        position_d, values_d = maybe_downsample(
            position, values, max_points=int(task["max_points"]), seed=int(task["seed"])
        )
    explicit_range = task.get("color_range")
    color_range = (
        _velocity_color_range(cfg, explicit_range)
        if field == "speed"
        else _scalar_color_range(values, explicit_range, field_default)
    )
    marker_size_raw = task["marker_size"]
    marker_size = None if marker_size_raw is None else float(marker_size_raw)

    if mode == "grid-slice":
        stats = plot_xy_grid_slice(
            position=position_d,
            speed=values_d,
            cfg=cfg,
            time_value=float(task["time_value"]),
            output_path=output_path,
            color_range=color_range,
            color_label=color_label,
        )
    elif mode == "xy-slice":
        stats = plot_xy_slice(
            position=position_d,
            speed=values_d,
            cfg=cfg,
            time_value=float(task["time_value"]),
            output_path=output_path,
            slice_half_width=float(task["slice_half_width"]),
            marker_size=marker_size,
            color_range=color_range,
            color_label=color_label,
        )
    else:
        stats = plot_3d_cloud(
            position=position_d,
            speed=values_d,
            cfg=cfg,
            time_value=float(task["time_value"]),
            output_path=output_path,
            marker_size=marker_size,
            elev=float(task["elev"]),
            azim=float(task["azim"]),
            color_range=color_range,
            color_label=color_label,
        )

    return {
        "frame_no": int(task["frame_no"]),
        "total_frames": int(task["total_frames"]),
        "field": field,
        "fluid_file": fluid_file.name,
        "output_name": output_path.name,
        "n_total": n_total,
        "n_rendered": int(position_d.shape[0]),
        "speed_max": float(stats.get("speed_max", 0.0)),
        "speed_mean": float(stats.get("speed_mean", 0.0)),
        "time_value": float(task["time_value"]),
    }


# ---------------------------------------------------------------------------
# CLI / main
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Render per-VTP flow-field images for the 3D Eulerian "
                    "fully compressible cylinder case.",
    )
    parser.add_argument("--data_dir", type=Path, default=Path("output"), help="Directory containing VTP files")
    parser.add_argument("--results_dir", type=Path, default=Path("results"), help="Output directory for rendered PNGs")
    parser.add_argument("--config", type=Path, default=Path("config.ini"), help="Path to config.ini")
    parser.add_argument(
        "--mode",
        choices=("grid-slice", "3d-cloud", "xy-slice"),
        default=None,
        help="Rendering mode: native fixed-grid XY field, full 3D scatter, or mid-span XY scatter. "
             "Defaults to [postprocess] flow_field_mode in config.ini.",
    )
    parser.add_argument("--fluid_prefix", type=str, default=DEFAULT_FLUID_PREFIX, help="Filename prefix for fluid VTPs")
    parser.add_argument("--velocity_key", type=str, default=DEFAULT_VELOCITY_KEY, help="Point-data field name for velocity")
    parser.add_argument(
        "--fields",
        type=str,
        default="speed,mach,pressure,density",
        help="Comma-separated scalars to render: speed, mach, pressure, density, vorticity. "
             "Compressible fields are normalised by the freestream (p/p_inf, rho/rho_inf); "
             "'vorticity' requires the solver to write a vorticity array.",
    )
    parser.add_argument(
        "--frames_dir_name",
        type=str,
        default=DEFAULT_FRAMES_DIRNAME,
        help="Sub-directory under --results_dir where PNGs are written.",
    )
    parser.add_argument(
        "--max_points",
        type=int,
        default=None,
        help="Cap on rendered fluid particles (downsampled uniformly if exceeded). 0 disables. "
             "Defaults to [postprocess] flow_field_max_points in config.ini.",
    )
    parser.add_argument(
        "--marker_size",
        type=float,
        default=None,
        help="Override scatter marker area in points squared. By default it is estimated "
             "from one particle's volume (dp^3) and the plot scale.",
    )
    parser.add_argument("--elev", type=float, default=22.0, help="3D camera elevation (degrees, only for 3d-cloud mode).")
    parser.add_argument("--azim", type=float, default=-58.0, help="3D camera azimuth (degrees, only for 3d-cloud mode).")
    parser.add_argument(
        "--slice_half_width",
        type=float,
        default=None,
        help="Half-thickness of XY slice (units of dp). Defaults to [postprocess] xy_slice_half_width_dp.",
    )
    parser.add_argument(
        "--vmin",
        type=float,
        default=None,
        help="Override lower color-bar bound (default 0).",
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=None,
        help="Override upper color-bar bound (default 1.6 * u_inf).",
    )
    parser.add_argument(
        "--include_initial",
        action="store_true",
        help="Render the initial condition VTP frame (ite_0000000000.vtp) as well. "
             "Defaults to NOT include initial unless this flag is set or "
             "[postprocess] flow_field_skip_initial=false.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, (os.cpu_count() or 4) - 1),
        help="Number of frames to render in parallel.",
    )
    parser.add_argument(
        "--overwrite_existing",
        action="store_true",
        help="Re-render frames whose PNG already exists.",
    )
    parser.add_argument("--verbose", action="store_true", help="Print detailed per-frame logs.")
    return parser


def _postprocess_defaults(config_path: Path) -> Dict[str, object]:
    """Read the ``[postprocess]`` section of ``config.ini`` as a dict.

    Missing keys are filled with the script's hard-coded defaults so the
    caller always receives a complete set of overrides.
    """
    defaults: Dict[str, object] = {
        "flow_field_mode": "3d-cloud",
        "flow_field_max_points": 200000,
        "flow_field_skip_initial": True,
        "xy_slice_half_width_dp": 0.5,
    }
    if not config_path.exists():
        return defaults
    parser = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    parser.read(config_path, encoding="utf-8")
    if "postprocess" not in parser:
        return defaults
    section = parser["postprocess"]
    if "flow_field_mode" in section:
        defaults["flow_field_mode"] = section["flow_field_mode"].strip()
    if "flow_field_max_points" in section:
        try:
            defaults["flow_field_max_points"] = int(section["flow_field_max_points"])
        except ValueError:
            pass
    if "flow_field_skip_initial" in section:
        defaults["flow_field_skip_initial"] = parse_bool_like(section["flow_field_skip_initial"])
    if "xy_slice_half_width_dp" in section:
        try:
            defaults["xy_slice_half_width_dp"] = float(section["xy_slice_half_width_dp"])
        except ValueError:
            pass
    return defaults


def parse_bool_like(value: str) -> bool:
    """Parse a forgiving boolean string (true/false/yes/no/1/0, case insensitive)."""
    norm = str(value).strip().lower()
    if norm in {"true", "1", "yes", "on"}:
        return True
    if norm in {"false", "0", "no", "off"}:
        return False
    raise ValueError(f"cannot parse boolean value: {value!r}")


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    cfg = load_case_config(args.config)
    pp_defaults = _postprocess_defaults(args.config)

    # Resolve CLI / config defaults before doing any work.
    if args.mode is None:
        mode = str(pp_defaults["flow_field_mode"]).strip().lower()
        if mode not in {"grid-slice", "3d-cloud", "xy-slice"}:
            raise ValueError(
                "[postprocess].flow_field_mode must be 'grid-slice', '3d-cloud', or "
                f"'xy-slice', got {mode!r}"
            )
        args.mode = mode
    if args.max_points is None:
        args.max_points = int(pp_defaults["flow_field_max_points"])
    if args.max_points < 0:
        raise ValueError(f"--max_points must be >= 0, got {args.max_points}")
    if args.marker_size is not None and args.marker_size <= 0.0:
        raise ValueError(f"--marker_size must be > 0, got {args.marker_size}")
    if args.slice_half_width is None:
        args.slice_half_width = float(pp_defaults["xy_slice_half_width_dp"])
    if args.slice_half_width < 0:
        raise ValueError(f"--slice_half_width must be >= 0, got {args.slice_half_width}")
    skip_initial_cfg = bool(pp_defaults["flow_field_skip_initial"])
    # CLI flag takes precedence when explicitly set; argparse defaults are False.
    args.include_initial = bool(args.include_initial) or not skip_initial_cfg
    color_range: Optional[Tuple[float, float]]
    if args.vmin is not None or args.vmax is not None:
        if args.vmin is None or args.vmax is None:
            raise ValueError("--vmin and --vmax must be supplied together")
        color_range = (float(args.vmin), float(args.vmax))
    else:
        color_range = None

    if not args.data_dir.exists():
        raise FileNotFoundError(f"data_dir not found: {args.data_dir}")

    fluid_files = sorted(args.data_dir.glob(f"{args.fluid_prefix}*.vtp"))
    if not fluid_files:
        raise FileNotFoundError(f"no {args.fluid_prefix}*.vtp files found in {args.data_dir}")
    if not args.include_initial:
        fluid_files = [p for p in fluid_files if "ite_0000000000" not in p.name]
    if not fluid_files:
        raise FileNotFoundError("all fluid VTPs were filtered out; pass --include_initial to keep t=0 frame")

    frames_meta: List[Dict[str, object]] = []
    for f in fluid_files:
        idx = extract_frame_index(f)
        time_value = read_time_value(f)
        if time_value is None:
            # Fall back to "unknown time" so the file still renders.
            time_value = float("nan")
        frames_meta.append(
            {
                "fluid_file": f,
                "frame_index": idx,
                "time_value": time_value,
            }
        )
    frames_meta.sort(key=lambda m: int(m["frame_index"]))

    output_dir = args.results_dir / args.frames_dir_name
    output_dir.mkdir(parents=True, exist_ok=True)

    known_fields = ("speed", "mach", "pressure", "density", "vorticity")
    fields = [f.strip().lower() for f in str(args.fields).split(",") if f.strip()]
    if not fields:
        raise ValueError("--fields must name at least one scalar")
    unknown = [f for f in fields if f not in known_fields]
    if unknown:
        raise ValueError(f"unknown --fields entries {unknown}; choose from {list(known_fields)}")
    if color_range is not None and len(fields) > 1:
        raise ValueError("--vmin/--vmax apply to a single field; pass one entry in --fields")

    cfg_dict = cfg.to_dict()
    tasks: List[Dict[str, object]] = []
    skipped = 0
    total = len(frames_meta) * len(fields)
    frame_no = 0
    for meta in frames_meta:
        t_value = float(meta["time_value"])
        time_tag = "t_unknown" if math.isnan(t_value) else _format_time_tag(t_value)
        for field in fields:
            frame_no += 1
            output_name = f"{args.mode}_{field}_{time_tag}.png"
            output_path = output_dir / output_name
            if output_path.exists() and not args.overwrite_existing:
                skipped += 1
                if args.verbose:
                    print(f"[FRAME {frame_no}/{total}] Skip existing: {output_path.name}", flush=True)
                continue
            tasks.append(
                {
                    "frame_no": frame_no,
                    "total_frames": total,
                    "fluid_file": str(meta["fluid_file"]),
                    "output_path": str(output_path),
                    "time_value": t_value,
                    "cfg_dict": cfg_dict,
                    "field": field,
                    "velocity_key": args.velocity_key,
                    "max_points": int(args.max_points),
                    "seed": 20260708 + frame_no,
                    "mode": args.mode,
                    "marker_size": args.marker_size,
                    "elev": float(args.elev),
                    "azim": float(args.azim),
                    "slice_half_width": float(args.slice_half_width) * float(cfg.dp),
                    "color_range": color_range,
                }
            )

    if args.verbose:
        print(f"[INFO] Mode: {args.mode}  Fields: {fields}", flush=True)
        print(
            f"[INFO] Freestream: Ma={cfg.mach_inf:.3f} u_inf={cfg.u_inf:.6f} "
            f"p_inf={cfg.p_inf:.6f} rho_inf={cfg.rho_inf:.3f} gamma={cfg.gamma:.3f}",
            flush=True,
        )
        print(
            f"[INFO] Speed color range: {color_range if color_range else _velocity_color_range(cfg, None)}",
            flush=True,
        )
        print(f"[INFO] Output dir: {output_dir}", flush=True)
        print(f"[INFO] Queued frames for rendering: {len(tasks)} (skipped existing: {skipped})", flush=True)
        print(f"[INFO] Parallel workers: {max(1, int(args.workers))}", flush=True)
        for meta in frames_meta:
            t = float(meta["time_value"])
            print(
                f"[INFO] frame_idx={int(meta['frame_index'])}  "
                f"fluid={Path(str(meta['fluid_file'])).name}  "
                f"t={'unknown' if math.isnan(t) else f'{t:.6f}s'}",
                flush=True,
            )

    processed = 0
    if tasks:
        workers = max(1, int(args.workers))
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futures = {ex.submit(render_frame_task, task): task for task in tasks}
            for fut in as_completed(futures):
                result = fut.result()
                processed += 1
                if args.verbose:
                    print(
                        f"[FRAME {result['frame_no']}/{result['total_frames']}] "
                        f"t={result['time_value']:.4f}s  field={result['field']}  "
                        f"particles={result['n_total']:,} rendered={result['n_rendered']:,}  "
                        f"max={result['speed_max']:.4f} mean={result['speed_mean']:.4f}  "
                        f"-> {result['output_name']}",
                        flush=True,
                    )

    print(
        f"[OK] {args.mode} frames written to {output_dir} "
        f"(processed={processed}, skipped={skipped}, total={total})",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
