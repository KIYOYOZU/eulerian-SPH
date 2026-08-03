#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
postprocess_cylinder_cl_cd_curve.py
===================================

Plot Cd(t) / Cl(t) time-series curves for the 3D Eulerian LG cylinder case.

The CSV input is produced by ``postprocess_cl_cd.py`` and contains one row
per logged force timestep with the columns::

    time, fx_pressure, fy_pressure, fx_viscous, fy_viscous,
    fx_total, fy_total, cd, cl

Two plots are produced under ``results/`` by default:

* ``cylinder_cl_cd_curve.png`` -- combined two-panel figure of Cd and Cl vs
  time, with a moving-average overlay (window in seconds, configurable).
* ``cylinder_force_components.png`` -- decomposed pressure vs viscous
  drag/lift contributions (useful when diagnosing near-zero force regimes).

A JSON sidecar with summary statistics (mean Cd/Cl, peak |Cl| for Strouhal
post-processing, RMS fluctuations) is also written.

Example
-------
::

    # Default: read results/cylinder_3d_cl_cd.csv, write into results/
    python postprocess_cylinder_cl_cd_curve.py

    # Custom CSV / output location, no moving-average overlay
    python postprocess_cylinder_cl_cd_curve.py \
        --csv results/cylinder_3d_cl_cd.csv \
        --results_dir results \
        --smooth_window 0
"""

from __future__ import annotations

import argparse
import configparser
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class ForceCurveStats:
    """Aggregated statistics over the Cd/Cl time series."""

    n_samples: int
    time_min: float
    time_max: float
    cd_mean: float
    cd_std: float
    cl_mean: float
    cl_std: float
    cl_rms: float
    cd_min: float
    cd_max: float
    cl_min: float
    cl_max: float


# ---------------------------------------------------------------------------
# CSV loading
# ---------------------------------------------------------------------------

REQUIRED_COLUMNS = (
    "time",
    "fx_pressure",
    "fy_pressure",
    "fx_viscous",
    "fy_viscous",
    "fx_total",
    "fy_total",
    "cd",
    "cl",
)


def load_cl_cd_csv(csv_path: Path) -> Dict[str, np.ndarray]:
    """Load the Cd/Cl CSV produced by ``postprocess_cl_cd.py`` into arrays."""
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    with csv_path.open("r", encoding="utf-8", newline="") as f:
        header_line = f.readline().strip()
    if not header_line:
        raise ValueError(f"CSV is empty: {csv_path}")
    header = [tok.strip() for tok in header_line.split(",")]
    missing = [col for col in REQUIRED_COLUMNS if col not in header]
    if missing:
        raise ValueError(
            f"CSV {csv_path} is missing required columns: {missing}; "
            f"found: {header}"
        )

    data = np.genfromtxt(csv_path, delimiter=",", names=True, dtype=np.float64)
    if data.ndim == 0:  # single row -> 0-d; reshape into 1-row structured
        data = np.atleast_1d(data)
    if data.size == 0:
        raise ValueError(f"CSV has no data rows: {csv_path}")

    out: Dict[str, np.ndarray] = {}
    for col in REQUIRED_COLUMNS:
        out[col] = np.asarray(data[col], dtype=np.float64)

    # sort by time defensively
    order = np.argsort(out["time"], kind="stable")
    for col in REQUIRED_COLUMNS:
        out[col] = out[col][order]
    return out


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def moving_average(values: np.ndarray, times: np.ndarray, window: float) -> np.ndarray:
    """Trailing moving average over a fixed time window (seconds).

    A non-positive window disables smoothing and returns a copy of the
    original series.
    """
    if window <= 0.0 or values.size <= 1:
        return values.astype(np.float64, copy=True)

    smoothed = np.empty_like(values, dtype=np.float64)
    for i in range(values.size):
        t_now = times[i]
        mask = times >= (t_now - window)
        smoothed[i] = float(np.mean(values[mask]))
    return smoothed


def compute_stats(data: Dict[str, np.ndarray]) -> ForceCurveStats:
    cd = data["cd"]
    cl = data["cl"]
    time = data["time"]
    return ForceCurveStats(
        n_samples=int(cd.size),
        time_min=float(time.min()),
        time_max=float(time.max()),
        cd_mean=float(cd.mean()),
        cd_std=float(cd.std(ddof=0)),
        cl_mean=float(cl.mean()),
        cl_std=float(cl.std(ddof=0)),
        cl_rms=float(np.sqrt(np.mean(np.square(cl - cl.mean())))),
        cd_min=float(cd.min()),
        cd_max=float(cd.max()),
        cl_min=float(cl.min()),
        cl_max=float(cl.max()),
    )


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _save_figure(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight")


def plot_cl_cd_combined(
    data: Dict[str, np.ndarray],
    cfg: Optional[Dict[str, float]],
    smooth_window: float,
    output_path: Path,
) -> None:
    """Two-panel Cd(t) / Cl(t) curve."""
    time = data["time"]
    cd_raw = data["cd"]
    cl_raw = data["cl"]
    cd_smooth = moving_average(cd_raw, time, smooth_window)
    cl_smooth = moving_average(cl_raw, time, smooth_window)

    fig, (ax_cd, ax_cl) = plt.subplots(2, 1, figsize=(8.8, 6.4), sharex=True)
    fig.suptitle("3D Eulerian Cylinder -- Cd / Cl time history", y=0.985)

    ax_cd.plot(time, cd_raw, color="#1f77b4", alpha=0.45, linewidth=0.8, label="Cd raw")
    if smooth_window > 0.0:
        ax_cd.plot(time, cd_smooth, color="#1f77b4", linewidth=1.6,
                   label=f"Cd moving avg (window = {smooth_window:.3f}s)")
    ax_cd.set_ylabel(r"$C_d$")
    ax_cd.grid(True, alpha=0.25)
    ax_cd.legend(loc="best", fontsize=9)

    ax_cl.plot(time, cl_raw, color="#d62728", alpha=0.45, linewidth=0.8, label="Cl raw")
    if smooth_window > 0.0:
        ax_cl.plot(time, cl_smooth, color="#d62728", linewidth=1.6,
                   label=f"Cl moving avg (window = {smooth_window:.3f}s)")
    ax_cl.axhline(0.0, color="black", linewidth=0.4, linestyle="--", alpha=0.6)
    ax_cl.set_xlabel("t (s)")
    ax_cl.set_ylabel(r"$C_l$")
    ax_cl.grid(True, alpha=0.25)
    ax_cl.legend(loc="best", fontsize=9)

    if cfg is not None:
        info = (
            f"$Re$ = {cfg.get('re', float('nan')):.0f},  "
            f"$\\rho_0$ = {cfg.get('rho0_f', float('nan')):.3f},  "
            f"$u_\\infty$ = {cfg.get('u_f', float('nan')):.3f},  "
            f"D = {cfg.get('D', float('nan')):.3f}"
        )
        fig.text(0.5, 0.94, info, ha="center", fontsize=9)

    _save_figure(fig, output_path)
    plt.close(fig)


def plot_force_components(
    data: Dict[str, np.ndarray],
    cfg: Optional[Dict[str, float]],
    smooth_window: float,
    output_path: Path,
) -> None:
    """Two-panel decomposed drag/lift split into pressure vs viscous parts."""
    time = data["time"]
    scale_denom = 0.5 * float(cfg.get("rho0_f", 1.0)) * float(cfg.get("u_f", 1.0)) ** 2 * float(cfg.get("D", 1.0)) * float(cfg.get("DW", 1.0))
    if scale_denom <= 0.0:
        scale_denom = 1.0

    cd_p = data["fx_pressure"] / scale_denom
    cd_v = data["fx_viscous"] / scale_denom
    cl_p = data["fy_pressure"] / scale_denom
    cl_v = data["fy_viscous"] / scale_denom

    cd_p_smooth = moving_average(cd_p, time, smooth_window)
    cd_v_smooth = moving_average(cd_v, time, smooth_window)
    cl_p_smooth = moving_average(cl_p, time, smooth_window)
    cl_v_smooth = moving_average(cl_v, time, smooth_window)

    fig, (ax_d, ax_l) = plt.subplots(2, 1, figsize=(8.8, 6.4), sharex=True)
    fig.suptitle("3D Eulerian Cylinder -- Drag / Lift components", y=0.985)

    ax_d.plot(time, cd_p, color="#1f77b4", alpha=0.35, linewidth=0.7, label=r"$C_d$ pressure")
    ax_d.plot(time, cd_v, color="#2ca02c", alpha=0.35, linewidth=0.7, label=r"$C_d$ viscous")
    if smooth_window > 0.0:
        ax_d.plot(time, cd_p_smooth, color="#1f77b4", linewidth=1.4, label=r"$C_d$ pressure (smooth)")
        ax_d.plot(time, cd_v_smooth, color="#2ca02c", linewidth=1.4, label=r"$C_d$ viscous (smooth)")
    ax_d.set_ylabel(r"$C_d$")
    ax_d.grid(True, alpha=0.25)
    ax_d.legend(loc="best", fontsize=9, ncol=2)

    ax_l.plot(time, cl_p, color="#d62728", alpha=0.35, linewidth=0.7, label=r"$C_l$ pressure")
    ax_l.plot(time, cl_v, color="#9467bd", alpha=0.35, linewidth=0.7, label=r"$C_l$ viscous")
    if smooth_window > 0.0:
        ax_l.plot(time, cl_p_smooth, color="#d62728", linewidth=1.4, label=r"$C_l$ pressure (smooth)")
        ax_l.plot(time, cl_v_smooth, color="#9467bd", linewidth=1.4, label=r"$C_l$ viscous (smooth)")
    ax_l.axhline(0.0, color="black", linewidth=0.4, linestyle="--", alpha=0.6)
    ax_l.set_xlabel("t (s)")
    ax_l.set_ylabel(r"$C_l$")
    ax_l.grid(True, alpha=0.25)
    ax_l.legend(loc="best", fontsize=9, ncol=2)

    _save_figure(fig, output_path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def load_case_context(config_path: Path) -> Dict[str, float]:
    """Read the case context for plot annotations (rho0, u_f, D, DW, Re)."""
    if not config_path.exists():
        return {}

    parser = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    parser.read(config_path, encoding="utf-8")
    info: Dict[str, float] = {}
    try:
        info["rho0_f"] = float(parser["physical"]["rho0_f"])
        info["u_f"] = float(parser["physical"]["u_f"])
        info["re"] = float(parser["physical"]["re"])
        radius = float(parser["geometry"]["cylinder_radius"])
        info["D"] = 2.0 * radius
        info["DW"] = float(parser["geometry"]["dw"])
    except (KeyError, ValueError):
        pass
    return info


def load_postprocess_defaults(config_path: Path) -> Dict[str, object]:
    """Read the ``[postprocess]`` section of ``config.ini`` for CLI defaults."""
    defaults: Dict[str, object] = {
        "cl_cd_csv": Path("results") / "cylinder_3d_cl_cd.csv",
        "cd_cl_smooth_window": 0.0,
    }
    if not config_path.exists():
        return defaults
    parser = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    parser.read(config_path, encoding="utf-8")
    if "postprocess" not in parser:
        return defaults
    section = parser["postprocess"]
    if "cl_cd_csv" in section:
        defaults["cl_cd_csv"] = Path(section["cl_cd_csv"].strip())
    if "cd_cl_smooth_window" in section:
        try:
            defaults["cd_cl_smooth_window"] = float(section["cd_cl_smooth_window"])
        except ValueError:
            pass
    return defaults


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Plot Cd/Cl time-series for the 3D Eulerian LG cylinder case.",
    )
    parser.add_argument("--csv", type=Path, default=None,
                        help="Path to cylinder_3d_cl_cd.csv (produced by postprocess_cl_cd.py). "
                             "Defaults to [postprocess] cl_cd_csv in config.ini.")
    parser.add_argument("--results_dir", type=Path, default=Path("results"),
                        help="Output directory for PNG / JSON sidecars.")
    parser.add_argument("--config", type=Path, default=Path("config.ini"),
                        help="Path to case config.ini (used for plot annotations).")
    parser.add_argument(
        "--smooth_window",
        type=float,
        default=None,
        help="Trailing moving-average window (seconds). 0 disables smoothing. "
             "Defaults to [postprocess] cd_cl_smooth_window in config.ini.",
    )
    parser.add_argument("--no_components", action="store_true",
                        help="Skip the pressure/viscous decomposition plot.")
    parser.add_argument("--verbose", action="store_true", help="Print per-plot logs.")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    pp_defaults = load_postprocess_defaults(args.config)

    if args.csv is None:
        args.csv = Path(str(pp_defaults["cl_cd_csv"]))
    if args.smooth_window is None:
        args.smooth_window = float(pp_defaults["cd_cl_smooth_window"])
    if args.smooth_window < 0:
        raise ValueError(f"--smooth_window must be >= 0, got {args.smooth_window}")

    data = load_cl_cd_csv(args.csv)
    cfg_info = load_case_context(args.config)
    stats = compute_stats(data)

    args.results_dir.mkdir(parents=True, exist_ok=True)

    cl_cd_path = args.results_dir / "cylinder_cl_cd_curve.png"
    plot_cl_cd_combined(data, cfg_info, args.smooth_window, cl_cd_path)
    if args.verbose:
        print(f"[OK] Wrote {cl_cd_path}", flush=True)

    if not args.no_components:
        components_path = args.results_dir / "cylinder_force_components.png"
        plot_force_components(data, cfg_info, args.smooth_window, components_path)
        if args.verbose:
            print(f"[OK] Wrote {components_path}", flush=True)

    summary_path = args.results_dir / "cylinder_cl_cd_stats.json"
    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(asdict(stats), f, indent=2, ensure_ascii=False)
    if args.verbose:
        print(f"[OK] Wrote {summary_path}", flush=True)

    print(
        f"[STATS] Cd mean={stats.cd_mean:.4f} ± {stats.cd_std:.4f}; "
        f"Cl mean={stats.cl_mean:.4f} ± {stats.cl_std:.4f} (RMS={stats.cl_rms:.4f}); "
        f"N={stats.n_samples} samples in t=[{stats.time_min:.3f}, {stats.time_max:.3f}]s",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())