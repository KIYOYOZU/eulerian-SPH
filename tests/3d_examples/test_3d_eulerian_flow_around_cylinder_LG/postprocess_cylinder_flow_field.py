#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
postprocess_cylinder_flow_field.py
==================================

Render per-VTP flow-field images for the 3D Eulerian LG cylinder case.

The script reads ``Cylinder3DFluid_*.vtp`` frames produced by
``test_3d_eulerian_flow_around_cylinder_LG``, optionally overlays the
``Cylinder_*.vtp`` wall particles as a faint outline, and writes a PNG
per VTP frame under ``results/flow_field_frames/``.

The default view is a 3D scatter of velocity magnitude on the particle
cloud plus the cylinder surface; use ``--mode xy-slice`` to switch to a
2D mid-span XY slice coloured by velocity magnitude, which is much
cheaper to render for smoke / quick visual sanity checks.

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

    # Only render a single mid-span XY slice
    python postprocess_cylinder_flow_field.py --mode xy-slice
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
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)


# ---------------------------------------------------------------------------
# Defaults / configuration parsing
# ---------------------------------------------------------------------------

DEFAULT_FRAMES_DIRNAME = "flow_field_frames"
DEFAULT_FLUID_PREFIX = "Cylinder3DFluid_"
DEFAULT_WALL_PREFIX = "Cylinder_"
DEFAULT_VELOCITY_KEY = "Velocity"


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
    rho0_f: float
    u_f: float
    re: float

    @property
    def D(self) -> float:
        return 2.0 * self.cylinder_radius

    def to_dict(self) -> Dict[str, float]:
        """Return init-only fields; the derived ``D`` property is excluded so
        the dict can be fed straight back into ``__init__`` after a process
        hop (e.g. ProcessPoolExecutor worker).
        """
        return {
            "DL": self.DL,
            "DH": self.DH,
            "DW": self.DW,
            "dp": self.dp,
            "cylinder_center_x": self.cylinder_center_x,
            "cylinder_center_y": self.cylinder_center_y,
            "cylinder_radius": self.cylinder_radius,
            "rho0_f": self.rho0_f,
            "u_f": self.u_f,
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
        rho0_f=parse_float(need("physical", "rho0_f"), "rho0_f"),
        u_f=parse_float(need("physical", "u_f"), "u_f"),
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

def load_fluid_frame(vtp_file: Path, velocity_key: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(position, velocity, speed)`` for one fluid VTP frame."""
    mesh = pv.read(vtp_file)
    if velocity_key not in mesh.point_data:
        raise KeyError(
            f"velocity key '{velocity_key}' missing in {vtp_file.name}; "
            f"available: {list(mesh.point_data.keys())}"
        )
    position = np.asarray(mesh.points, dtype=np.float64)
    velocity = np.asarray(mesh.point_data[velocity_key], dtype=np.float64)
    if velocity.ndim != 2 or velocity.shape[1] < 3:
        raise ValueError(f"unexpected velocity shape {velocity.shape} in {vtp_file.name}")
    speed = np.linalg.norm(velocity[:, :3], axis=1)
    del mesh
    return position, velocity, speed


def load_wall_frame(vtp_file: Path) -> Optional[np.ndarray]:
    """Return wall positions for a Cylinder_*.vtp, or ``None`` if missing."""
    if not vtp_file.exists():
        return None
    mesh = pv.read(vtp_file)
    pos = np.asarray(mesh.points, dtype=np.float64)
    del mesh
    return pos


def maybe_downsample(
    position: np.ndarray,
    speed: np.ndarray,
    max_points: int,
    seed: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Stratified downsample to keep the scatter under ``max_points``."""
    n = position.shape[0]
    if max_points <= 0 or n <= max_points:
        return position, speed
    rng = np.random.default_rng(seed)
    sample_idx = np.sort(rng.choice(n, size=max_points, replace=False))
    return position[sample_idx], speed[sample_idx]


# ---------------------------------------------------------------------------
# Plotting primitives
# ---------------------------------------------------------------------------

def _format_time_tag(time_value: float) -> str:
    """Filesystem-safe time tag, e.g. ``t00012p345678s``."""
    return f"{time_value:010.6f}".replace(".", "p") + "s"


def _velocity_color_range(
    cfg: CylinderCaseConfig,
    explicit: Optional[Tuple[float, float]],
) -> Tuple[float, float]:
    if explicit is not None:
        vmin, vmax = explicit
        if vmin >= vmax:
            raise ValueError(f"invalid color range [{vmin}, {vmax}]")
        return float(vmin), float(vmax)
    if cfg.u_f <= 0.0:
        raise ValueError(f"u_f must be positive to derive default color range, got {cfg.u_f}")
    # Default: 0 .. 1.6 * u_f (Re=100 cylinder flow rarely exceeds ~1.5 * u_inf)
    return 0.0, 1.6 * float(cfg.u_f)


def _draw_cylinder_wireframe(ax, cfg: CylinderCaseConfig, color: str = "black", lw: float = 0.6) -> None:
    """Draw a circle representing the cylinder at z=mid-span."""
    theta = np.linspace(0.0, 2.0 * np.pi, 96)
    x = cfg.cylinder_center_x + cfg.cylinder_radius * np.cos(theta)
    y = cfg.cylinder_center_y + cfg.cylinder_radius * np.sin(theta)
    ax.plot(x, y, color=color, linewidth=lw, zorder=5)


def plot_3d_cloud(
    *,
    position: np.ndarray,
    speed: np.ndarray,
    cfg: CylinderCaseConfig,
    wall_position: Optional[np.ndarray],
    time_value: float,
    output_path: Path,
    marker_size: float,
    elev: float,
    azim: float,
    color_range: Tuple[float, float],
) -> Dict[str, float]:
    """Render a 3D scatter of velocity magnitude + cylinder wall outline."""
    fig = plt.figure(figsize=(8.8, 5.6))
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    ax.set_position([0.01, 0.02, 0.84, 0.92])

    sc = ax.scatter(
        position[:, 0],
        position[:, 1],
        position[:, 2],
        c=speed,
        cmap="turbo",
        marker="o",
        s=marker_size,
        vmin=color_range[0],
        vmax=color_range[1],
        edgecolors="none",
        linewidths=0.0,
        depthshade=False,
        alpha=0.9,
    )

    # Draw a cylinder wireframe as a stack of circles along the spanwise axis
    # so the obstacle is visually obvious even when the wall particles are not
    # loaded. We always render this (it is cheap and informative).
    z_levels = np.linspace(0.0, cfg.DW, 7)
    theta = np.linspace(0.0, 2.0 * np.pi, 64)
    circle_x = cfg.cylinder_center_x + cfg.cylinder_radius * np.cos(theta)
    circle_y = cfg.cylinder_center_y + cfg.cylinder_radius * np.sin(theta)
    for z in z_levels:
        ax.plot(
            circle_x,
            circle_y,
            np.full_like(circle_x, z),
            color="black",
            linewidth=0.6,
            alpha=0.7,
        )
    # Two vertical "spines" to make the cylinder read as a 3D solid.
    for ang in (0.0, 0.5 * np.pi, np.pi, 1.5 * np.pi):
        ax.plot(
            [cfg.cylinder_center_x + cfg.cylinder_radius * np.cos(ang)] * 2,
            [cfg.cylinder_center_y + cfg.cylinder_radius * np.sin(ang)] * 2,
            [0.0, cfg.DW],
            color="black",
            linewidth=0.6,
            alpha=0.7,
        )

    if wall_position is not None and wall_position.size > 0:
        ax.scatter(
            wall_position[:, 0],
            wall_position[:, 1],
            wall_position[:, 2],
            c="black",
            s=0.05,
            alpha=0.4,
            depthshade=False,
            linewidths=0.0,
        )

    ax.set_title(
        f"3D Flow Field (t = {time_value:.4f}s, Re = {cfg.re:.0f})",
        pad=4,
    )
    ax.view_init(elev=elev, azim=azim)
    ax.set_box_aspect([cfg.DL, cfg.DH, cfg.DW])
    ax.set_axis_off()

    cbar = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.01, shrink=0.84)
    cbar.set_label("Velocity magnitude (m/s)")

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
    marker_size: float,
    color_range: Tuple[float, float],
) -> Dict[str, float]:
    """Render an XY slice of the flow field at z = DW/2 within ±slice_half_width."""
    z_mid = 0.5 * cfg.DW
    mask = np.abs(position[:, 2] - z_mid) <= slice_half_width
    pos = position[mask]
    spd = speed[mask]

    fig, ax = plt.subplots(figsize=(8.8, 4.4))
    ax.set_position([0.06, 0.10, 0.84, 0.82])

    sc = ax.scatter(
        pos[:, 0],
        pos[:, 1],
        c=spd,
        cmap="turbo",
        marker="o",
        s=marker_size,
        vmin=color_range[0],
        vmax=color_range[1],
        edgecolors="none",
        linewidths=0.0,
    )

    _draw_cylinder_wireframe(ax, cfg, color="black", lw=0.9)

    ax.set_xlim(0.0, cfg.DL)
    ax.set_ylim(0.0, cfg.DH)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (m, streamwise)")
    ax.set_ylabel("y (m, vertical)")
    ax.set_title(
        f"XY Slice @ z = {z_mid:.3f} m ± {slice_half_width:.3f} m  |  "
        f"t = {time_value:.4f}s  |  Re = {cfg.re:.0f}"
    )
    ax.grid(True, alpha=0.25, linewidth=0.4)

    cbar = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.015)
    cbar.set_label("Velocity magnitude (m/s)")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    return {
        "n_total": int(position.shape[0]),
        "n_slice": int(pos.shape[0]),
        "speed_max": float(np.nanmax(spd)) if spd.size else 0.0,
        "speed_mean": float(np.nanmean(spd)) if spd.size else 0.0,
    }


# ---------------------------------------------------------------------------
# Worker
# ---------------------------------------------------------------------------

def render_frame_task(task: Dict[str, object]) -> Dict[str, object]:
    """Render one frame in a worker process."""
    fluid_file = Path(str(task["fluid_file"]))
    wall_file = Path(str(task["wall_file"])) if task.get("wall_file") else None
    output_path = Path(str(task["output_path"]))
    cfg_dict = dict(task["cfg_dict"])
    cfg = CylinderCaseConfig(**cfg_dict)

    position, _velocity, speed = load_fluid_frame(fluid_file, velocity_key=str(task["velocity_key"]))
    n_total = int(position.shape[0])
    position_d, speed_d = maybe_downsample(
        position, speed, max_points=int(task["max_points"]), seed=int(task["seed"])
    )
    wall_position = load_wall_frame(wall_file) if wall_file is not None else None

    color_range = _velocity_color_range(cfg, task.get("color_range"))

    if str(task["mode"]) == "xy-slice":
        stats = plot_xy_slice(
            position=position_d,
            speed=speed_d,
            cfg=cfg,
            time_value=float(task["time_value"]),
            output_path=output_path,
            slice_half_width=float(task["slice_half_width"]),
            marker_size=float(task["marker_size"]),
            color_range=color_range,
        )
    else:
        stats = plot_3d_cloud(
            position=position_d,
            speed=speed_d,
            cfg=cfg,
            wall_position=wall_position,
            time_value=float(task["time_value"]),
            output_path=output_path,
            marker_size=float(task["marker_size"]),
            elev=float(task["elev"]),
            azim=float(task["azim"]),
            color_range=color_range,
        )

    return {
        "frame_no": int(task["frame_no"]),
        "total_frames": int(task["total_frames"]),
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
        description="Render per-VTP flow-field images for the 3D Eulerian LG cylinder case.",
    )
    parser.add_argument("--data_dir", type=Path, default=Path("output"), help="Directory containing VTP files")
    parser.add_argument("--results_dir", type=Path, default=Path("results"), help="Output directory for rendered PNGs")
    parser.add_argument("--config", type=Path, default=Path("config.ini"), help="Path to config.ini")
    parser.add_argument(
        "--mode",
        choices=("3d-cloud", "xy-slice"),
        default=None,
        help="Rendering mode: full 3D scatter or mid-span XY slice. "
             "Defaults to [postprocess] flow_field_mode in config.ini.",
    )
    parser.add_argument("--fluid_prefix", type=str, default=DEFAULT_FLUID_PREFIX, help="Filename prefix for fluid VTPs")
    parser.add_argument("--wall_prefix", type=str, default=DEFAULT_WALL_PREFIX, help="Filename prefix for cylinder wall VTPs")
    parser.add_argument("--velocity_key", type=str, default=DEFAULT_VELOCITY_KEY, help="Point-data field name for velocity")
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
        default=0.35,
        help="Scatter marker size for the chosen --mode.",
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
        help="Override upper color-bar bound (default 1.6 * u_f).",
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
        if mode not in {"3d-cloud", "xy-slice"}:
            raise ValueError(
                f"[postprocess].flow_field_mode must be '3d-cloud' or 'xy-slice', got {mode!r}"
            )
        args.mode = mode
    if args.max_points is None:
        args.max_points = int(pp_defaults["flow_field_max_points"])
    if args.max_points < 0:
        raise ValueError(f"--max_points must be >= 0, got {args.max_points}")
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

    wall_files_by_index: Dict[int, Path] = {}
    for wf in args.data_dir.glob(f"{args.wall_prefix}*.vtp"):
        try:
            wall_files_by_index[extract_frame_index(wf)] = wf
        except ValueError:
            continue

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
                "wall_file": wall_files_by_index.get(idx),
                "time_value": time_value,
            }
        )
    frames_meta.sort(key=lambda m: int(m["frame_index"]))

    output_dir = args.results_dir / args.frames_dir_name
    output_dir.mkdir(parents=True, exist_ok=True)

    cfg_dict = cfg.to_dict()
    tasks: List[Dict[str, object]] = []
    skipped = 0
    total = len(frames_meta)
    for frame_no, meta in enumerate(frames_meta, start=1):
        t_value = float(meta["time_value"])
        time_tag = "t_unknown" if math.isnan(t_value) else _format_time_tag(t_value)
        output_name = f"{args.mode}_{time_tag}.png"
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
                "wall_file": str(meta["wall_file"]) if meta["wall_file"] else None,
                "output_path": str(output_path),
                "time_value": t_value,
                "cfg_dict": cfg_dict,
                "velocity_key": args.velocity_key,
                "max_points": int(args.max_points),
                "seed": 20260708 + frame_no,
                "mode": args.mode,
                "marker_size": float(args.marker_size),
                "elev": float(args.elev),
                "azim": float(args.azim),
                "slice_half_width": float(args.slice_half_width) * float(cfg.dp),
                "color_range": color_range,
            }
        )

    if args.verbose:
        print(f"[INFO] Mode: {args.mode}", flush=True)
        print(
            f"[INFO] Color range: {color_range if color_range else _velocity_color_range(cfg, None)}",
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
                        f"t={result['time_value']:.4f}s  particles={result['n_total']:,} "
                        f"rendered={result['n_rendered']:,}  speed_max={result['speed_max']:.3f} "
                        f"speed_mean={result['speed_mean']:.3f}  -> {result['output_name']}",
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