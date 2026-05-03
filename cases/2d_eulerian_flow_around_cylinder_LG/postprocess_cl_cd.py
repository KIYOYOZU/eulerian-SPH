#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Post-process cylinder force histories into drag and lift coefficients.

The SPHinXsys output files used here are the total viscous and pressure forces
acting on the cylinder from the fluid. Drag keeps the output x-force sign,
while lift is inverted by default to match the requested convention.
"""

import argparse
import csv
import math
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


RHO0 = 1.0            # reference density
U_REF = 0.06          # m/s, freestream velocity
CYLINDER_DIAMETER = 0.05  # m, D = 2*R = 2*0.025


def read_case_config_reference_values(case_dir):
    config_path = case_dir / "config.ini"
    if not config_path.exists():
        return {}

    values = {}
    section = ""
    with config_path.open("r", encoding="utf-8") as config_file:
        for raw_line in config_file:
            line = raw_line.split("#", 1)[0].split(";", 1)[0].strip()
            if not line:
                continue
            if line.startswith("[") and line.endswith("]"):
                section = line[1:-1].strip().lower()
                continue
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            values[(section, key.strip().lower())] = value.strip()

    reference_values = {}
    if ("fluid", "rho0_f") in values:
        reference_values["rho0"] = float(values[("fluid", "rho0_f")])
    if ("fluid", "u_f") in values:
        reference_values["u_ref"] = float(values[("fluid", "u_f")])
    if ("geometry", "cylinder_radius") in values:
        reference_values["diameter"] = 2.0 * float(values[("geometry", "cylinder_radius")])
    return reference_values


def read_force_file(path):
    if not path.exists():
        raise FileNotFoundError(f"Force file not found: {path}")

    records = []
    with path.open("r", encoding="utf-8") as force_file:
        header = force_file.readline().strip()
        if not header:
            raise ValueError(f"Force file is empty: {path}")

        for line_number, line in enumerate(force_file, start=2):
            stripped = line.strip()
            if not stripped:
                continue

            columns = stripped.split()
            if len(columns) < 3:
                raise ValueError(
                    f"Expected at least 3 columns in {path} line {line_number}, got: {stripped}"
                )

            values = [float(column.replace("-nan(ind)", "nan").replace("nan(ind)", "nan")) for column in columns[:3]]
            if not all(math.isfinite(value) for value in values):
                continue

            records.append((values[0], values[1], values[2]))

    if not records:
        raise ValueError(f"No numeric force records found in: {path}")

    return records


def combine_force_histories(viscous_records, pressure_records, tolerance):
    if len(viscous_records) != len(pressure_records):
        raise ValueError(
            "Viscous and pressure force histories have different lengths: "
            f"{len(viscous_records)} vs {len(pressure_records)}"
        )

    combined = []
    for index, (viscous, pressure) in enumerate(zip(viscous_records, pressure_records), start=1):
        time_v, force_vx, force_vy = viscous
        time_p, force_px, force_py = pressure
        if abs(time_v - time_p) > tolerance:
            raise ValueError(
                f"Time mismatch at row {index}: viscous={time_v}, pressure={time_p}, "
                f"tolerance={tolerance}"
            )

        combined.append(
            {
                "time": 0.5 * (time_v + time_p),
                "force_x_viscous": force_vx,
                "force_y_viscous": force_vy,
                "force_x_pressure": force_px,
                "force_y_pressure": force_py,
                "force_x_total": force_vx + force_px,
                "force_y_total": force_vy + force_py,
            }
        )

    return combined


def read_fluid_wall_force_file(path):
    if not path.exists():
        raise FileNotFoundError(f"Fluid-side wall force file not found: {path}")

    records = []
    with path.open("r", encoding="utf-8") as force_file:
        header = force_file.readline().strip()
        if not header:
            raise ValueError(f"Fluid-side wall force file is empty: {path}")

        for line_number, line in enumerate(force_file, start=2):
            stripped = line.strip()
            if not stripped:
                continue

            columns = stripped.split()
            if len(columns) < 6:
                raise ValueError(
                    f"Expected at least 6 columns in {path} line {line_number}, got: {stripped}"
                )

            values = [float(column.replace("-nan(ind)", "nan").replace("nan(ind)", "nan")) for column in columns[:6]]
            if not all(math.isfinite(value) for value in values):
                continue

            records.append(
                {
                    "time": values[0],
                    "force_x_pressure": values[2],
                    "force_y_pressure": values[3],
                    "force_x_viscous": values[4],
                    "force_y_viscous": values[5],
                    "force_x_total": values[2] + values[4],
                    "force_y_total": values[3] + values[5],
                }
            )

    if not records:
        raise ValueError(f"No numeric fluid-side wall force records found in: {path}")

    return records


def compute_coefficients(records, rho0, u_ref, diameter, invert_lift, coefficient_force_sign=1.0):
    dynamic_force_scale = 0.5 * rho0 * u_ref * u_ref * diameter
    if dynamic_force_scale == 0.0:
        raise ValueError("Coefficient normalization scale is zero.")

    lift_sign = -1.0 if invert_lift else 1.0
    for record in records:
        record["cd"] = coefficient_force_sign * record["force_x_total"] / dynamic_force_scale
        record["cl"] = lift_sign * coefficient_force_sign * record["force_y_total"] / dynamic_force_scale

    return dynamic_force_scale


def write_csv(records, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "time",
        "force_x_viscous",
        "force_y_viscous",
        "force_x_pressure",
        "force_y_pressure",
        "force_x_total",
        "force_y_total",
        "cd",
        "cl",
    ]

    with output_path.open("w", encoding="utf-8", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            writer.writerow({name: f"{record[name]:.12g}" for name in fieldnames})


def compute_tail_summary(records, tail_fraction):
    if not 0.0 < tail_fraction <= 1.0:
        raise ValueError(f"Tail fraction must be in (0, 1], got {tail_fraction}.")

    start_time = records[0]["time"]
    end_time = records[-1]["time"]
    tail_start_time = start_time + (1.0 - tail_fraction) * (end_time - start_time)
    tail_records = [record for record in records if record["time"] >= tail_start_time]
    if not tail_records:
        raise ValueError("No records found in tail averaging window.")

    return {
        "tail_fraction": tail_fraction,
        "tail_start_time": tail_start_time,
        "tail_end_time": end_time,
        "tail_record_count": len(tail_records),
        "cd_tail_mean": sum(record["cd"] for record in tail_records) / len(tail_records),
    }


def compute_tail_mean(records, value_name, tail_fraction):
    if not 0.0 < tail_fraction <= 1.0:
        raise ValueError(f"Tail fraction must be in (0, 1], got {tail_fraction}.")

    start_time = records[0]["time"]
    end_time = records[-1]["time"]
    tail_start_time = start_time + (1.0 - tail_fraction) * (end_time - start_time)
    tail_records = [record for record in records if record["time"] >= tail_start_time]
    if not tail_records:
        raise ValueError("No records found in tail averaging window.")

    return {
        "tail_fraction": tail_fraction,
        "tail_start_time": tail_start_time,
        "tail_end_time": end_time,
        "tail_record_count": len(tail_records),
        "mean": sum(record[value_name] for record in tail_records) / len(tail_records),
    }


def compute_strouhal_summary(records, value_name, tail_fraction, diameter, u_ref):
    if not 0.0 < tail_fraction <= 1.0:
        raise ValueError(f"Tail fraction must be in (0, 1], got {tail_fraction}.")
    if diameter <= 0.0 or u_ref <= 0.0:
        raise ValueError("Strouhal normalization requires positive diameter and u_ref.")

    start_time = records[0]["time"]
    end_time = records[-1]["time"]
    tail_start_time = start_time + (1.0 - tail_fraction) * (end_time - start_time)
    tail_records = [record for record in records if record["time"] >= tail_start_time]
    if len(tail_records) < 3:
        return {
            "st_signal": value_name,
            "st_tail_start_time": tail_start_time,
            "st_tail_end_time": end_time,
            "st_peak_count": 0,
            "st_mean_period": math.nan,
            "st_frequency": math.nan,
            "st": math.nan,
        }

    peak_times = []
    for previous, current, following in zip(tail_records, tail_records[1:], tail_records[2:]):
        if current[value_name] > previous[value_name] and current[value_name] >= following[value_name]:
            peak_times.append(current["time"])

    if len(peak_times) < 2:
        mean_period = math.nan
        frequency = math.nan
        strouhal = math.nan
    else:
        periods = [right - left for left, right in zip(peak_times, peak_times[1:]) if right > left]
        if periods:
            mean_period = sum(periods) / len(periods)
            frequency = 1.0 / mean_period
            strouhal = frequency * diameter / u_ref
        else:
            mean_period = math.nan
            frequency = math.nan
            strouhal = math.nan

    return {
        "st_signal": value_name,
        "st_tail_start_time": tail_start_time,
        "st_tail_end_time": end_time,
        "st_peak_count": len(peak_times),
        "st_mean_period": mean_period,
        "st_frequency": frequency,
        "st": strouhal,
    }


def format_float(value):
    if isinstance(value, float) and math.isnan(value):
        return "nan"
    return f"{value:.12g}"


def write_summary(summary, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "tail_fraction",
        "tail_start_time",
        "tail_end_time",
        "tail_record_count",
        "cd_tail_mean",
    ]

    with output_path.open("w", encoding="utf-8", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(
            {
                "tail_fraction": format_float(summary["tail_fraction"]),
                "tail_start_time": format_float(summary["tail_start_time"]),
                "tail_end_time": format_float(summary["tail_end_time"]),
                "tail_record_count": summary["tail_record_count"],
                "cd_tail_mean": format_float(summary["cd_tail_mean"]),
            }
        )


def write_dual_summary(cd_summary, cl_summary, st_summary, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "tail_fraction",
        "tail_start_time",
        "tail_end_time",
        "tail_record_count",
        "cd_tail_mean",
        "cl_tail_mean",
        "st_signal",
        "st_peak_count",
        "st_mean_period",
        "st_frequency",
        "st",
    ]

    with output_path.open("w", encoding="utf-8", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(
            {
                "tail_fraction": format_float(cd_summary["tail_fraction"]),
                "tail_start_time": format_float(cd_summary["tail_start_time"]),
                "tail_end_time": format_float(cd_summary["tail_end_time"]),
                "tail_record_count": cd_summary["tail_record_count"],
                "cd_tail_mean": format_float(cd_summary["mean"]),
                "cl_tail_mean": format_float(cl_summary["mean"]),
                "st_signal": st_summary["st_signal"],
                "st_peak_count": st_summary["st_peak_count"],
                "st_mean_period": format_float(st_summary["st_mean_period"]),
                "st_frequency": format_float(st_summary["st_frequency"]),
                "st": format_float(st_summary["st"]),
            }
        )


def plot_coefficients(records, output_path, title):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    times = [record["time"] for record in records]
    cd_values = [record["cd"] for record in records]
    cl_values = [record["cl"] for record in records]

    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.grid": True,
            "grid.alpha": 0.3,
            "lines.linewidth": 1.6,
        }
    )

    fig, axes = plt.subplots(2, 1, figsize=(9.0, 6.4), sharex=True)
    axes[0].plot(times, cd_values, color="#1f77b4")
    axes[0].set_ylabel(r"$C_D$")
    axes[0].set_title(title)

    axes[1].plot(times, cl_values, color="#d62728")
    axes[1].set_xlabel("Time")
    axes[1].set_ylabel(r"$C_L$")

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _parse_vtp_data_array(data_array):
    text = data_array.text or ""
    values = np.fromstring(text, sep=" ")
    components = int(data_array.attrib.get("NumberOfComponents", "1"))
    if components > 1:
        if values.size % components != 0:
            raise ValueError(
                f"DataArray {data_array.attrib.get('Name', '<unnamed>')} has "
                f"{values.size} values, not divisible by {components} components."
            )
        values = values.reshape((-1, components))
    return values


def _read_vtp_time(root):
    time_array = root.find(".//FieldData/DataArray[@Name='TimeValue']")
    if time_array is None:
        return None

    values = np.fromstring(time_array.text or "", sep=" ")
    if values.size == 0:
        return None
    return float(values[0])


def _read_vtp_piece_fields(vtp_path):
    root = ET.parse(vtp_path).getroot()
    piece = root.find(".//Piece")
    if piece is None:
        raise ValueError(f"No Piece element found in VTP file: {vtp_path}")

    positions_array = piece.find("./Points/DataArray[@Name='Position']")
    if positions_array is None:
        raise ValueError(f"No Position array found in VTP file: {vtp_path}")

    fields = {"Position": _parse_vtp_data_array(positions_array)}
    point_data = piece.find("./PointData")
    if point_data is not None:
        for data_array in point_data.findall("./DataArray"):
            name = data_array.attrib.get("Name")
            if name:
                fields[name] = _parse_vtp_data_array(data_array)

    fields["TimeValue"] = _read_vtp_time(root)
    fields["PieceName"] = piece.attrib.get("Name", vtp_path.stem).strip()
    return fields


def find_latest_flow_vtp(output_dir, body_name):
    if not output_dir.exists():
        raise FileNotFoundError(f"Output directory not found: {output_dir}")

    if body_name:
        candidates = sorted(output_dir.glob(f"{body_name}_*.vtp"))
    else:
        candidates = sorted(output_dir.glob("*.vtp"))
    if not candidates:
        raise FileNotFoundError(f"No VTP files found in: {output_dir}")

    valid_candidates = []
    for candidate in candidates:
        try:
            fields = _read_vtp_piece_fields(candidate)
        except (ET.ParseError, ValueError):
            continue
        if "Velocity" in fields and "Pressure" in fields:
            valid_candidates.append((candidate, fields))

    if not valid_candidates:
        raise ValueError(
            f"No VTP file with both Velocity and Pressure arrays found in: {output_dir}"
        )

    candidates_with_time = [
        (candidate, fields) for candidate, fields in valid_candidates if fields.get("TimeValue") is not None
    ]
    if candidates_with_time:
        return max(candidates_with_time, key=lambda item: item[1]["TimeValue"])

    return valid_candidates[-1]


def _set_equal_aspect_limits(axis, x, y):
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlim(float(np.min(x)), float(np.max(x)))
    axis.set_ylim(float(np.min(y)), float(np.max(y)))


def _compute_particle_marker_size(point_count):
    if point_count <= 0:
        return 1.0
    return float(np.clip(80000.0 / point_count, 0.15, 2.5))


def _plot_scalar_on_particles(axis, x, y, values, title, colorbar_label, cmap, marker_size):
    scatter = axis.scatter(
        x,
        y,
        c=values,
        s=marker_size,
        marker="o",
        cmap=cmap,
        linewidths=0.0,
        edgecolors="none",
        rasterized=True,
    )
    axis.set_title(title)
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    axis.set_facecolor("#f8f9fb")
    return scatter


def _select_refinement_values(fields, point_count, finite_mask):
    volumetric_measure = fields.get("VolumetricMeasure")
    if volumetric_measure is not None:
        values = np.asarray(volumetric_measure).reshape(-1)
        if len(values) == len(finite_mask):
            values = values[finite_mask]
        if len(values) == point_count and np.all(np.isfinite(values)):
            return np.sqrt(values), "Particle spacing proxy", r"$\sqrt{V}$", "magma_r"

    smoothing_length_ratio = fields.get("SmoothingLengthRatio")
    if smoothing_length_ratio is not None:
        values = np.asarray(smoothing_length_ratio).reshape(-1)
        if len(values) == len(finite_mask):
            values = values[finite_mask]
        if len(values) == point_count and np.all(np.isfinite(values)):
            return values, "Smoothing length ratio", r"$h/h_0$", "magma"

    indicator = fields.get("Indicator")
    if indicator is not None:
        values = np.asarray(indicator).reshape(-1)
        if len(values) == len(finite_mask):
            values = values[finite_mask]
        if len(values) == point_count and np.all(np.isfinite(values)):
            return values, "Particle indicator", "indicator", "magma"

    return None


def compute_particle_vorticity(x, y, velocity):
    if velocity.ndim != 2 or velocity.shape[1] < 2:
        raise ValueError("Velocity array must have at least two components.")

    triangulation = mtri.Triangulation(x, y)
    if triangulation.triangles.size == 0:
        raise ValueError("Not enough non-collinear points to compute vorticity.")

    u_gradient = mtri.LinearTriInterpolator(triangulation, velocity[:, 0]).gradient(x, y)
    v_gradient = mtri.LinearTriInterpolator(triangulation, velocity[:, 1]).gradient(x, y)
    return np.asarray(v_gradient[0]) - np.asarray(u_gradient[1]), triangulation


def plot_latest_flow_fields(output_dir, results_dir, body_name, output_name):
    vtp_path, fields = find_latest_flow_vtp(output_dir, body_name)
    positions = fields["Position"]
    velocity = fields["Velocity"]
    pressure = fields["Pressure"]

    if positions.ndim != 2 or positions.shape[1] < 2:
        raise ValueError(f"Position array in {vtp_path} must have at least two components.")
    if pressure.ndim != 1:
        pressure = np.asarray(pressure).reshape(-1)
    if not (len(positions) == len(velocity) == len(pressure)):
        raise ValueError(
            f"VTP array length mismatch in {vtp_path}: "
            f"Position={len(positions)}, Velocity={len(velocity)}, Pressure={len(pressure)}"
        )

    finite_mask = (
        np.isfinite(positions[:, 0])
        & np.isfinite(positions[:, 1])
        & np.isfinite(velocity[:, 0])
        & np.isfinite(velocity[:, 1])
        & np.isfinite(pressure)
    )
    positions = positions[finite_mask]
    velocity = velocity[finite_mask]
    pressure = pressure[finite_mask]
    if len(positions) < 3:
        raise ValueError(f"Not enough finite fluid points in VTP file: {vtp_path}")

    x = positions[:, 0]
    y = positions[:, 1]
    speed = np.linalg.norm(velocity[:, :2], axis=1)
    vorticity, _ = compute_particle_vorticity(x, y, velocity)
    vort_lo, vort_hi = np.nanpercentile(vorticity, [1, 99])
    vorticity = np.clip(vorticity, vort_lo, vort_hi)
    refinement_panel = _select_refinement_values(fields, len(x), finite_mask)
    marker_size = _compute_particle_marker_size(len(x))

    results_dir.mkdir(parents=True, exist_ok=True)
    figure_path = results_dir / output_name

    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.grid": False,
        }
    )
    panel_count = 4 if refinement_panel is not None else 3
    if panel_count == 4:
        fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.2), constrained_layout=True)
        axes = axes.ravel()
    else:
        fig, axes = plt.subplots(1, 3, figsize=(15.6, 4.8), constrained_layout=True)
    time_value = fields.get("TimeValue")
    time_text = f", t={time_value:.8g}" if time_value is not None else ""
    fig.suptitle(
        f"{fields.get('PieceName', vtp_path.stem)} final particle flow fields{time_text}"
    )

    panels = [
        (speed, "Velocity magnitude", "|U|", "viridis"),
        (pressure, "Pressure", "p", "coolwarm"),
        (vorticity, "Vorticity", r"$\omega_z$", "RdBu_r"),
    ]
    if refinement_panel is not None:
        panels.append(refinement_panel)

    for axis, (values, title, label, cmap) in zip(axes, panels):
        contour = _plot_scalar_on_particles(axis, x, y, values, title, label, cmap, marker_size)
        _set_equal_aspect_limits(axis, x, y)
        fig.colorbar(contour, ax=axis, label=label, shrink=0.88)

    fig.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return vtp_path, figure_path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot drag and lift coefficients for the 2D Eulerian cylinder LG case."
    )
    parser.add_argument(
        "--case-dir",
        default=None,
        help="Case directory containing output/ and results/. Defaults to the current working directory.",
    )
    parser.add_argument("--output-dir", default="output", help="Directory containing SPHinXsys force .dat files.")
    parser.add_argument("--results-dir", default="results", help="Directory for generated CSV and plots.")
    parser.add_argument(
        "--rho0",
        type=float,
        default=None,
        help="Reference density. Defaults to config.ini [fluid] rho0_f, then script fallback.",
    )
    parser.add_argument(
        "--u-ref",
        type=float,
        default=None,
        help="Reference velocity. Defaults to config.ini [fluid] u_f, then script fallback.",
    )
    parser.add_argument(
        "--diameter",
        type=float,
        default=None,
        help="Cylinder diameter. Defaults to 2*config.ini [geometry] cylinder_radius, then script fallback.",
    )
    parser.add_argument(
        "--keep-lift-sign",
        action="store_true",
        help="Do not invert the y-force sign when computing lift coefficient.",
    )
    parser.add_argument(
        "--time-tolerance",
        type=float,
        default=1.0e-8,
        help="Allowed time mismatch between viscous and pressure force files.",
    )
    parser.add_argument(
        "--tail-fraction",
        type=float,
        default=0.4,
        help="Fraction of the final time span used to average Cd/Cl and estimate St.",
    )
    parser.add_argument(
        "--st-signal",
        choices=("cl", "cd"),
        default="cl",
        help="Coefficient signal used for Strouhal estimation from peak-to-peak period.",
    )
    parser.add_argument(
        "--skip-fluid-side",
        action="store_true",
        help="Skip post-processing output/fluid_wall_force.dat.",
    )
    parser.add_argument(
        "--skip-flow-field",
        action="store_true",
        help="Skip plotting the final VTP velocity, pressure and vorticity fields.",
    )
    parser.add_argument(
        "--flow-vtp-body",
        default="WaterBlock",
        help="VTP body/file prefix used for flow-field plotting. Use an empty string to scan all VTP files.",
    )
    parser.add_argument(
        "--flow-field-output",
        default="final_flow_fields.png",
        help="Output figure name for the final VTP flow-field plot.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.case_dir is None:
        case_dir = Path.cwd()
    else:
        case_dir = Path(args.case_dir).expanduser()
        if not case_dir.is_absolute():
            case_dir = Path.cwd() / case_dir
        case_dir = case_dir.resolve()
    output_dir = (case_dir / args.output_dir).resolve()
    results_dir = (case_dir / args.results_dir).resolve()
    config_reference_values = read_case_config_reference_values(case_dir)
    rho0 = args.rho0 if args.rho0 is not None else config_reference_values.get("rho0", RHO0)
    u_ref = args.u_ref if args.u_ref is not None else config_reference_values.get("u_ref", U_REF)
    diameter = (
        args.diameter
        if args.diameter is not None
        else config_reference_values.get("diameter", CYLINDER_DIAMETER)
    )

    viscous_path = output_dir / "Cylinder_TotalViscousForceFromFluid.dat"
    pressure_path = output_dir / "Cylinder_TotalPressureForceFromFluid.dat"

    print(f"[INFO] Reading viscous force: {viscous_path}")
    viscous_records = read_force_file(viscous_path)
    print(f"[INFO] Reading pressure force: {pressure_path}")
    pressure_records = read_force_file(pressure_path)

    records = combine_force_histories(viscous_records, pressure_records, args.time_tolerance)
    force_scale = compute_coefficients(
        records,
        rho0=rho0,
        u_ref=u_ref,
        diameter=diameter,
        invert_lift=not args.keep_lift_sign,
    )

    csv_path = results_dir / "cl_cd_timeseries.csv"
    summary_path = results_dir / "cl_cd_summary.csv"
    figure_path = results_dir / "cl_cd.png"
    cd_summary = compute_tail_mean(records, "cd", args.tail_fraction)
    cl_summary = compute_tail_mean(records, "cl", args.tail_fraction)
    st_summary = compute_strouhal_summary(
        records,
        args.st_signal,
        args.tail_fraction,
        diameter=diameter,
        u_ref=u_ref,
    )
    write_csv(records, csv_path)
    write_dual_summary(cd_summary, cl_summary, st_summary, summary_path)
    plot_coefficients(records, figure_path, "2D Eulerian Flow Around Cylinder LG")

    print(f"[INFO] Normalization scale: 0.5*rho0*U^2*D = {force_scale:.12g}")
    print(f"[INFO] Reference values: rho0={rho0:.12g}, U={u_ref:.12g}, D={diameter:.12g}")
    print(f"[INFO] Lift sign: {'inverted' if not args.keep_lift_sign else 'kept'}")
    print(
        "[INFO] Final "
        f"{100.0 * cd_summary['tail_fraction']:.1f}% time-window Cd mean: "
        f"{cd_summary['mean']:.12g} "
        f"(t={cd_summary['tail_start_time']:.12g} to {cd_summary['tail_end_time']:.12g}, "
        f"n={cd_summary['tail_record_count']})"
    )
    print(f"[INFO] Final Cl mean: {cl_summary['mean']:.12g}")
    print(
        "[INFO] Strouhal estimate "
        f"from {st_summary['st_signal']} peaks: St={format_float(st_summary['st'])}, "
        f"period={format_float(st_summary['st_mean_period'])}, "
        f"peaks={st_summary['st_peak_count']}"
    )
    print(f"[INFO] Wrote CSV: {csv_path}")
    print(f"[INFO] Wrote summary: {summary_path}")
    print(f"[INFO] Wrote figure: {figure_path}")

    fluid_wall_path = output_dir / "fluid_wall_force.dat"
    if not args.skip_fluid_side:
        print(f"[INFO] Reading fluid-side wall force: {fluid_wall_path}")
        fluid_records = read_fluid_wall_force_file(fluid_wall_path)
        compute_coefficients(
            fluid_records,
            rho0=rho0,
            u_ref=u_ref,
            diameter=diameter,
            invert_lift=not args.keep_lift_sign,
            coefficient_force_sign=-1.0,
        )
        fluid_csv_path = results_dir / "cl_cd_fluid_side_timeseries.csv"
        fluid_summary_path = results_dir / "cl_cd_fluid_side_summary.csv"
        fluid_figure_path = results_dir / "cl_cd_fluid_side.png"
        fluid_cd_summary = compute_tail_mean(fluid_records, "cd", args.tail_fraction)
        fluid_cl_summary = compute_tail_mean(fluid_records, "cl", args.tail_fraction)
        fluid_st_summary = compute_strouhal_summary(
            fluid_records,
            args.st_signal,
            args.tail_fraction,
            diameter=diameter,
            u_ref=u_ref,
        )
        write_csv(fluid_records, fluid_csv_path)
        write_dual_summary(fluid_cd_summary, fluid_cl_summary, fluid_st_summary, fluid_summary_path)
        plot_coefficients(fluid_records, fluid_figure_path, "2D Eulerian Flow Around Cylinder LG: Fluid-Side Wall Force")
        print(
            "[INFO] Fluid-side final "
            f"{100.0 * fluid_cd_summary['tail_fraction']:.1f}% time-window Cd mean: "
            f"{fluid_cd_summary['mean']:.12g} "
            f"(t={fluid_cd_summary['tail_start_time']:.12g} to {fluid_cd_summary['tail_end_time']:.12g}, "
            f"n={fluid_cd_summary['tail_record_count']})"
        )
        print(f"[INFO] Fluid-side final Cl mean: {fluid_cl_summary['mean']:.12g}")
        print(
            "[INFO] Fluid-side Strouhal estimate "
            f"from {fluid_st_summary['st_signal']} peaks: St={format_float(fluid_st_summary['st'])}, "
            f"period={format_float(fluid_st_summary['st_mean_period'])}, "
            f"peaks={fluid_st_summary['st_peak_count']}"
        )
        print(f"[INFO] Wrote fluid-side CSV: {fluid_csv_path}")
        print(f"[INFO] Wrote fluid-side summary: {fluid_summary_path}")
        print(f"[INFO] Wrote fluid-side figure: {fluid_figure_path}")
    if not args.skip_flow_field:
        print(f"[INFO] Reading latest flow VTP from: {output_dir}")
        vtp_path, flow_figure_path = plot_latest_flow_fields(
            output_dir,
            results_dir,
            body_name=args.flow_vtp_body.strip(),
            output_name=args.flow_field_output,
        )
        print(f"[INFO] Read flow VTP: {vtp_path}")
        print(f"[INFO] Wrote flow-field figure: {flow_figure_path}")
    print("[INFO] Done.")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"[ERROR] {error}", file=sys.stderr)
        sys.exit(1)
