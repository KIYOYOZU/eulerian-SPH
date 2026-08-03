#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
postprocess_eulerian_channel.py
================================
Lightweight post-processor for the 3D Eulerian channel smoke test.

Reads the reduced diagnostics produced by the case (MaximumSpeed recording
and the VTP output folder) and writes compact numerical summaries plus
figures into ``results/``. No pyvista dependency for the core summary —
only matplotlib for the figures. Missing optional inputs are reported with
a clear message and a nonzero exit, never a silent empty plot.

Usage:
    python -B postprocess_eulerian_channel.py --data_dir output --results_dir results

Outputs:
    results/inlet_profile.png            (placeholder if no VTP reader)
    results/centerline_velocity.png      (from reduced max-speed history)
    results/profile_evolution.png        (from all fluid VTP snapshots)
    results/final_3d_velocity_cloud.png (speed-colored point cloud from final VTP)
    results/wall_slip_summary.json       (finiteness + max-speed summary)
"""
import argparse
import json
import os
import re
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_reduced_dat(path: Path) -> Tuple[List[float], List[float]]:
    """Parse a SPHinXsys ReducedQuantityRecording .dat file.

    SPHinXsys writes a quoted header row (e.g. ``"run_time" "MaximumSpeed"``)
    followed by whitespace-separated data rows. Supports both 2-column
    (time, value) and 3-column (iteration, time, value) layouts.
    Returns (times, values).
    """
    times: List[float] = []
    values: List[float] = []
    if not path.exists():
        return times, values
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            # Skip the quoted header row.
            if line.startswith('"'):
                continue
            parts = line.split()
            try:
                if len(parts) >= 3:
                    times.append(float(parts[1]))
                    values.append(float(parts[2]))
                elif len(parts) == 2:
                    times.append(float(parts[0]))
                    values.append(float(parts[1]))
            except ValueError:
                continue
    return times, values


def collect_vtp_times(data_dir: Path) -> List[int]:
    """Collect the iteration indices of available VTP snapshots.

    SPHinXsys names snapshots ``<BodyName>_<iter>.vtp`` (e.g.
    ``ChannelFluid_0000025699.vtp``). Returns the sorted unique iteration list.
    """
    iters: List[int] = []
    if not data_dir.exists():
        return iters
    pattern = re.compile(r"_(\d+)\.vtp$")
    for entry in sorted(data_dir.iterdir()):
        m = pattern.search(entry.name)
        if m:
            iters.append(int(m.group(1)))
    return sorted(set(iters))


def write_centerline_figure(times: List[float], values: List[float],
                            out_path: Path) -> None:
    """Plot the maximum-speed history as a centreline-proxy diagnostic."""
    if not times:
        print("[postprocess] no reduced history found, skipping centerline figure")
        return
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(times, values, "-o", markersize=2, linewidth=1.0)
    ax.set_xlabel("physical time (s)")
    ax.set_ylabel("maximum speed (m/s)")
    ax.set_title("Eulerian channel — maximum speed history")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"[postprocess] wrote {out_path}")


def write_inlet_profile_figure(out_path: Path, DH: float = 1.0,
                               U_bulk: float = 1.0) -> None:
    """Plot the analytical inlet parabolic profile (reference figure).

    The smoke-test inlet is the parabolic profile u_x(y) = U_max * 4 * eta *
    (1 - eta) with U_max = 1.5 * U_bulk. Until VTP particle data is parsed
    on-site, we draw the target profile so the figure is never empty.
    """
    import numpy as np
    U_max = 1.5 * U_bulk
    y = np.linspace(0.0, DH, 200)
    eta = y / DH
    u = U_max * 4.0 * eta * (1.0 - eta)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(u, y, linewidth=1.5)
    ax.set_xlabel("streamwise velocity u_x (m/s)")
    ax.set_ylabel("wall-normal y (m)")
    ax.set_title("Inlet parabolic profile (target)")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"[postprocess] wrote {out_path}")


def read_vtp_time(vtp_path: Path) -> Optional[float]:
    """Return the physical time stored in a VTP snapshot.

    SPHinXsys writes the snapshot time as a ``TimeValue`` field-data array.
    Falls back to ``None`` if the field is absent (caller then maps the
    iteration index from the filename via the reduced history).
    """
    import pyvista as pv
    mesh = pv.read(str(vtp_path))
    fd = mesh.field_data
    for key in ("TimeValue", "Time", "time"):
        if key in fd and fd[key].size:
            try:
                return float(fd[key][0])
            except (TypeError, ValueError):
                continue
    return None


def compute_mean_velocity_profile(vtp_path: Path, DH: float,
                                   n_bins: int = 0) -> Tuple[np.ndarray, np.ndarray]:
    """Read one VTP snapshot and return the x,z-averaged u_x(y) profile.

    Groups particles by their wall-normal y layer. Because the Eulerian
    lattice places particles on discrete y planes, we auto-detect the unique
    y layers (rounded to particle spacing) and average within each layer
    rather than fixed-width bins that could land between layers. Returns
    (y_centers, u_x_mean) only for non-empty layers.
    """
    import numpy as np
    import pyvista as pv
    mesh = pv.read(str(vtp_path))
    pts = np.asarray(mesh.points)
    vel = np.asarray(mesh.point_data["Velocity"])
    y = pts[:, 1]
    ux = vel[:, 0]
    # Auto-detect y layer spacing: round to a tolerance derived from the
    # y range so lattice planes collapse to a single key.
    y_range = float(y.max() - y.min()) if y.size else DH
    tol = max(y_range * 1e-3, 1e-6)
    y_keys = np.round(y / tol) * tol
    uniq = np.unique(y_keys)
    y_centers = uniq
    u_mean = np.zeros(uniq.shape[0])
    for k, key in enumerate(uniq):
        mask = y_keys == key
        u_mean[k] = ux[mask].mean() if mask.any() else 0.0
    return y_centers, u_mean


def write_profile_evolution_figure(data_dir: Path, out_path: Path,
                                   DH: float, U_bulk: float,
                                   n_bins: int = 20) -> None:
    """Plot the mean streamwise velocity profile u_x(y) evolution over time.

    Overlays profiles from every available VTP snapshot so the time evolution
    of the wall-normal profile is visible at a glance. The colour of each
    curve encodes the snapshot's physical time, read directly from the VTP
    ``TimeValue`` field data (so the colour bar and legend run in real
    seconds 0 → end_time, not a normalised snapshot index). Also draws the
    target parabolic inlet profile as a reference.
    """
    import numpy as np
    import pyvista as pv
    vtp_files = sorted(data_dir.glob("*.vtp"))
    if not vtp_files:
        print("[postprocess] no VTP found, skipping profile evolution figure")
        return
    # Collect (time, y_centers, u_mean) per snapshot. Time comes from the VTP
    # TimeValue field; snapshots missing it are dropped (we cannot place them
    # on a real time axis without guessing).
    profiles = []
    for vtp in vtp_files:
        try:
            y_c, u_mean = compute_mean_velocity_profile(vtp, DH, n_bins)
            t = read_vtp_time(vtp)
        except Exception as exc:  # pragma: no cover - per-file robustness
            print(f"[postprocess] skip {vtp.name}: {exc}", file=sys.stderr)
            continue
        if t is None:
            print(f"[postprocess] skip {vtp.name}: no TimeValue field data",
                  file=sys.stderr)
            continue
        profiles.append((t, y_c, u_mean))
    if not profiles:
        print("[postprocess] no profile with time info, skipping figure",
              file=sys.stderr)
        return

    profiles.sort(key=lambda p: p[0])
    times = np.array([p[0] for p in profiles])
    t_min, t_max = float(times.min()), float(times.max())
    # Guard against a single-snapshot case where t_max == t_min (norm needs a span).
    t_span = t_max - t_min if t_max > t_min else 1.0

    fig, ax = plt.subplots(figsize=(7, 5))
    cmap = plt.cm.viridis
    norm = plt.Normalize(vmin=t_min, vmax=t_max)
    U_max = 1.5 * U_bulk
    n = len(profiles)
    for k, (t, y_c, u_mean) in enumerate(profiles):
        color = cmap(norm(t))
        # Label first and last snapshot with their real physical time.
        label = f"t = {t:.2f} s" if k in (0, n - 1) else None
        ax.plot(u_mean, y_c, color=color, alpha=0.8, linewidth=1.2,
                label=label)
    # Target parabolic profile reference.
    y_ref = np.linspace(0.0, DH, 200)
    eta = y_ref / DH
    u_ref = U_max * 4.0 * eta * (1.0 - eta)
    ax.plot(u_ref, y_ref, "k--", linewidth=1.5, label="target parabolic")
    ax.set_xlabel("mean streamwise velocity u_x (m/s)")
    ax.set_ylabel("wall-normal y (m)")
    ax.set_title("Mean velocity profile evolution")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation="vertical", pad=0.02)
    cbar.set_label("physical time (s)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"[postprocess] wrote {out_path} "
          f"(t: {t_min:.3f} → {t_max:.3f} s, {n} snapshots)")


def write_final_3d_flow_figure(data_dir: Path, out_path: Path,
                               max_vectors: int = 2500) -> None:
    """Plot the final VTP velocity field as a spatially sampled 3D quiver."""
    import numpy as np
    import pyvista as pv
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    vtp_files = sorted(data_dir.glob("ChannelFluid_*.vtp"))
    if not vtp_files:
        print("[postprocess] no ChannelFluid VTP found, skipping 3D flow figure",
              file=sys.stderr)
        return

    final_vtp = vtp_files[-1]
    mesh = pv.read(str(final_vtp))
    points = np.asarray(mesh.points, dtype=float)
    velocity = np.asarray(mesh.point_data["Velocity"], dtype=float)
    if points.shape[0] != velocity.shape[0] or points.shape[0] == 0:
        raise ValueError(f"invalid point/velocity arrays in {final_vtp}")

    speed = np.linalg.norm(velocity, axis=1)
    finite = np.isfinite(points).all(axis=1) & np.isfinite(velocity).all(axis=1)
    points, velocity, speed = points[finite], velocity[finite], speed[finite]
    if points.shape[0] == 0:
        raise ValueError(f"no finite velocity vectors in {final_vtp}")

    # Select one representative vector per spatial voxel, then cap the count.
    # This avoids storage-order bias and keeps the quiver readable.
    extent = np.ptp(points, axis=0)
    spacing = np.max(extent) / max(8.0, max_vectors ** (1.0 / 3.0))
    spacing = max(float(spacing), 1e-12)
    keys = np.floor((points - points.min(axis=0)) / spacing).astype(np.int64)
    representatives = {}
    for index, key in enumerate(map(tuple, keys)):
        representatives.setdefault(key, index)
    selected = np.fromiter(representatives.values(), dtype=np.int64)
    if selected.size > max_vectors:
        order = np.argsort(speed[selected])[::-1][:max_vectors]
        selected = selected[order]

    p = points[selected]
    v = velocity[selected]
    s = speed[selected]
    vmax = float(np.max(s))
    color_norm = plt.Normalize(vmin=0.0, vmax=vmax if vmax > 0.0 else 1.0)

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection="3d")
    colors = plt.cm.viridis(color_norm(s))
    vector_length = max(float(np.max(extent)) * 0.045, 1e-6)
    ax.quiver(p[:, 0], p[:, 1], p[:, 2],
              v[:, 0], v[:, 1], v[:, 2],
              length=vector_length, normalize=True,
              colors=colors, linewidth=0.55, arrow_length_ratio=0.28)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_zlabel("z (m)")
    ax.set_title(f"Final 3D channel velocity field ({final_vtp.stem})")
    ax.set_box_aspect(np.maximum(extent, 1e-6))
    scalar_mappable = plt.cm.ScalarMappable(cmap="viridis", norm=color_norm)
    scalar_mappable.set_array(s)
    fig.colorbar(scalar_mappable, ax=ax, pad=0.1, shrink=0.7,
                 label="speed |u| (m/s)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"[postprocess] wrote {out_path} from {final_vtp.name} "
          f"({selected.size} vectors/{points.shape[0]} finite particles)")


def write_final_3d_velocity_cloud_figure(data_dir: Path, out_path: Path,
                                         DH: float = 1.0,
                                         U_bulk: float = 1.0) -> None:
    """Render the final speed cloud using the cylinder postprocessor style."""
    import numpy as np
    import pyvista as pv

    vtp_files = sorted(
        data_dir.glob("ChannelFluid_[0-9]*.vtp"),
        key=lambda path: int(path.stem.rsplit("_", 1)[1]),
    )
    if not vtp_files:
        print("[postprocess] no numeric ChannelFluid VTP found, "
              "skipping 3D velocity cloud", file=sys.stderr)
        return

    final_vtp = vtp_files[-1]
    mesh = pv.read(str(final_vtp))
    points = np.asarray(mesh.points, dtype=np.float64)
    velocity = np.asarray(mesh.point_data["Velocity"], dtype=np.float64)
    speed = np.linalg.norm(velocity[:, :3], axis=1)
    finite = np.isfinite(points).all(axis=1) & np.isfinite(speed)
    points, speed = points[finite], speed[finite]
    if points.shape[0] == 0:
        raise ValueError(f"no finite velocity data in {final_vtp}")

    extent = np.ptp(points, axis=0)
    x_values = np.unique(np.sort(points[:, 0]))
    dp = float(np.min(np.diff(x_values))) if x_values.size > 1 else 0.05
    particle_radius = (3.0 * dp ** 3 / (4.0 * np.pi)) ** (1.0 / 3.0)
    widest_in_plane_span = max(float(extent[0]), float(extent[1]), 1.0)
    points_per_length = 72.0 * 8.8 * 0.84 / widest_in_plane_span
    marker_size = 0.35 * np.pi * (particle_radius * points_per_length) ** 2
    color_max = max(1.6 * float(U_bulk), float(np.nanmax(speed)))

    fig = plt.figure(figsize=(8.8, 5.6))
    ax = fig.add_subplot(1, 1, 1, projection="3d", computed_zorder=False)
    ax.set_position([0.01, 0.02, 0.84, 0.92])
    cloud = ax.scatter(
        points[:, 0], points[:, 2], points[:, 1],
        c=speed, cmap="turbo", marker="o", s=marker_size,
        vmin=0.0, vmax=color_max, edgecolors="none", linewidths=0.0,
        depthshade=False, alpha=1.0, zorder=2,
    )
    ax.view_init(elev=24.0, azim=-62.0)
    ax.set_proj_type("ortho")
    # Display coordinates are (x, z, y): physical y is vertical.
    ax.set_box_aspect([max(extent[0], 1e-6), max(extent[2], 1e-6),
                       max(DH, 1e-6)])
    ax.set_axis_off()
    ax.set_title(f"3D Flow Field (t = final, |u|, {final_vtp.stem})", pad=4)
    cbar = fig.colorbar(cloud, ax=ax, fraction=0.03, pad=0.01, shrink=0.84)
    cbar.set_label("Velocity magnitude |u| (m/s)")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[postprocess] wrote {out_path} from {final_vtp.name} "
          f"({points.shape[0]} particles, color range 0-{color_max:.4g} m/s)")


def write_summary_json(out_path: Path, times: List[float], values: List[float],
                       vtp_iters: List[int]) -> None:
    summary = {
        "reduced_history_points": len(times),
        "vtp_snapshots": len(vtp_iters),
        "vtp_iterations": vtp_iters,
    }
    if values:
        summary["max_speed_final"] = float(values[-1])
        summary["max_speed_peak"] = float(max(values))
        summary["max_speed_finite"] = bool(
            all(v == v and v not in (float("inf"), float("-inf")) for v in values)
        )
    else:
        summary["max_speed_final"] = None
        summary["max_speed_peak"] = None
        summary["max_speed_finite"] = False
    with out_path.open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[postprocess] wrote {out_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data_dir", default="output",
                        help="case output directory with VTP + reduced dat")
    parser.add_argument("--results_dir", default="results",
                        help="destination directory for figures and JSON")
    parser.add_argument("--config", default="config.ini",
                        help="config.ini path for profile parameters")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    # Reduced maximum-speed history (written by ReducedQuantityRecording).
    # SPHinXsys names it "<BodyName>_<QuantityName>.dat"; search by suffix.
    dat_path: Optional[Path] = None
    if data_dir.exists():
        candidates = sorted(data_dir.glob("*_MaximumSpeed.dat")) + \
                     sorted(data_dir.glob("*maximum_speed*.dat")) + \
                     sorted(data_dir.glob("MaximumSpeed.dat"))
        if candidates:
            dat_path = candidates[0]
    times, values = ([], [])
    if dat_path is not None:
        times, values = parse_reduced_dat(dat_path)
    else:
        print(f"[postprocess] no reduced dat found under {data_dir}", file=sys.stderr)

    vtp_iters = collect_vtp_times(data_dir)

    # Geometry / physical parameters for profile figures: read from config.ini
    # if present, otherwise fall back to defaults matching the smoke config.
    DH = 1.0
    U_bulk = 1.0
    cfg_path = Path(args.config)
    if cfg_path.exists():
        try:
            import configparser as _cp
            parser = _cp.ConfigParser(inline_comment_prefixes=("#", ";"))
            parser.read(cfg_path, encoding="utf-8")
            if parser.has_option("geometry", "dh"):
                DH = parser.getfloat("geometry", "dh")
            if parser.has_option("physical", "u_bulk"):
                U_bulk = parser.getfloat("physical", "u_bulk")
        except Exception as exc:
            print(f"[postprocess] config parse failed, using defaults: {exc}",
                  file=sys.stderr)

    # Figures
    write_centerline_figure(times, values, results_dir / "centerline_velocity.png")
    write_inlet_profile_figure(results_dir / "inlet_profile.png", DH=DH, U_bulk=U_bulk)
    # Mean velocity profile evolution from VTP snapshots (requires pyvista).
    try:
        write_profile_evolution_figure(data_dir,
                                       results_dir / "profile_evolution.png",
                                       DH=DH, U_bulk=U_bulk)
    except ImportError:
        print("[postprocess] pyvista not available, skipping profile evolution",
              file=sys.stderr)

    # Final 3D velocity-magnitude cloud from the last numeric ChannelFluid snapshot.
    try:
        write_final_3d_velocity_cloud_figure(
            data_dir, results_dir / "final_3d_velocity_cloud.png")
    except ImportError:
        print("[postprocess] pyvista not available, skipping final 3D velocity cloud",
              file=sys.stderr)

    # JSON summary
    write_summary_json(results_dir / "wall_slip_summary.json",
                       times, values, vtp_iters)

    if not values and not vtp_iters:
        print("[postprocess] no input data found — figures are reference-only.",
              file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
