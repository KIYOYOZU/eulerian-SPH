#!/usr/bin/env python3
"""
postprocess_rotor67.py
Domain-STL post-processing diagnostics for the Rotor 67 smoke test.

Reads the VTP output written by test_3d_eulerian_rotor67 and produces:
  - mass-flux imbalance between inlet (z ~ z_in) and outlet (z ~ z_out) faces,
  - pressure ratio p_out / p_in (trend indicator),
  - a meridional (r-z) slice contour of pressure / Mach to qualitatively
    locate shock structures in the transonic passage.

Usage:
    python postprocess_rotor67.py [--case-dir CASE_DIR] [--output-dir OUTPUT_DIR]

Trend indicators only — this smoke gate does not claim absolute accuracy.
"""
import argparse
import json
import os
import re
import sys
from pathlib import Path

import numpy as np


def parse_config(case_dir: Path) -> dict:
    """Minimal INI reader producing section.key entries."""
    cfg = {}
    ini_path = case_dir / "config.ini"
    if not ini_path.exists():
        return cfg
    section = None
    with open(ini_path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.split("#", 1)[0].split(";", 1)[0].strip()
            if not line:
                continue
            if line.startswith("[") and line.endswith("]"):
                section = line[1:-1].lower()
                continue
            if "=" in line and section is not None:
                key, val = line.split("=", 1)
                cfg[f"{section}.{key.strip().lower()}"] = val.strip()
    return cfg


def to_float(cfg: dict, key: str, default: float = 0.0) -> float:
    try:
        return float(cfg.get(key, default))
    except (TypeError, ValueError):
        return default


def _resolve_case_path(case_dir: Path, raw: str | None) -> Path | None:
    if not raw:
        return None
    path = Path(raw)
    return path if path.is_absolute() else case_dir / path


def _mean_theta(vertices: np.ndarray) -> float:
    theta = np.arctan2(vertices[:, 1], vertices[:, 0])
    return float(np.arctan2(np.sin(theta).sum(), np.cos(theta).sum()))


def _positive_angle_diff(a: float, b: float) -> float:
    diff = a - b
    while diff < 0.0:
        diff += 2.0 * np.pi
    while diff > 2.0 * np.pi:
        diff -= 2.0 * np.pi
    return float(diff)


def load_domain_geometry(case_dir: Path, cfg: dict) -> dict:
    """Read domain-STL inlet/outlet z and periodic angle from named solids."""
    stl_path = _resolve_case_path(case_dir, cfg.get("geometry.domain_stl"))
    geometry = {
        "source": "config_fallback",
        "z_in": to_float(cfg, "geometry.z_in", -0.15),
        "z_out": to_float(cfg, "geometry.z_out", 0.18),
        "periodic_delta_theta": to_float(
            cfg, "rotor.delta_theta",
            2.0 * np.pi / max(int(to_float(cfg, "rotor.n_blades", 22)), 1)),
    }
    if stl_path is None or not stl_path.exists():
        geometry["domain_stl"] = str(stl_path) if stl_path is not None else None
        return geometry

    vertices_by_solid: dict[str, list[tuple[float, float, float]]] = {}
    current = None
    with open(stl_path, "r", encoding="utf-8") as f:
        for raw in f:
            parts = raw.strip().split()
            if not parts:
                continue
            if parts[0] == "solid":
                current = " ".join(parts[1:])
                vertices_by_solid.setdefault(current, [])
                continue
            if parts[0] == "endsolid":
                current = None
                continue
            if current is not None and parts[0] == "vertex" and len(parts) >= 4:
                vertices_by_solid[current].append(
                    (float(parts[1]), float(parts[2]), float(parts[3])))

    inlet = np.asarray(vertices_by_solid.get("rotor67_domain_inlet", []), dtype=float)
    outlet = np.asarray(vertices_by_solid.get("rotor67_domain_outlet", []), dtype=float)
    pmin = np.asarray(vertices_by_solid.get("rotor67_domain_periodic_min", []), dtype=float)
    pmax = np.asarray(vertices_by_solid.get("rotor67_domain_periodic_max", []), dtype=float)
    if len(inlet) and len(outlet):
        geometry["z_in"] = float(inlet[:, 2].min())
        geometry["z_out"] = float(outlet[:, 2].max())
        geometry["source"] = "domain_stl"
    if len(pmin) and len(pmax):
        geometry["periodic_delta_theta"] = _positive_angle_diff(_mean_theta(pmax), _mean_theta(pmin))
        geometry["source"] = "domain_stl"
    geometry["domain_stl"] = str(stl_path)
    return geometry


def load_last_frame(output_dir: Path):
    """Load the last RotorFluid_*.vtp frame."""
    try:
        import pyvista as pv
    except ImportError:
        print("ERROR: pyvista is required for VTP postprocess. Install via `pip install pyvista`.",
              file=sys.stderr)
        sys.exit(1)
    vtp_files = sorted(output_dir.glob("RotorFluid_*.vtp"))
    if not vtp_files:
        raise FileNotFoundError(f"No RotorFluid_*.vtp found in {output_dir}")
    frame_pattern = re.compile(r"^RotorFluid_(\d+)\.vtp$")
    numbered = []
    for path in vtp_files:
        match = frame_pattern.match(path.name)
        if match:
            numbered.append((int(match.group(1)), path))
    if numbered:
        _, latest = max(numbered, key=lambda item: item[0])
    else:
        latest = max(vtp_files, key=lambda path: path.stat().st_mtime)
    return pv.read(str(latest)), latest.name


def run_self_test() -> int:
    pts = np.array([[0.1, 0.0, 0.0], [0.0, 0.1, 0.0]], dtype=float)
    w = np.array([[0.0, -10.0, 100.0], [10.0, 0.0, 100.0]], dtype=float)
    u = absolute_velocity_from_relative(pts, w, 100.0)
    if not np.allclose(u[:, :2], 0.0):
        print("[self-test] absolute velocity conversion failed", file=sys.stderr)
        return 1
    p = np.array([101325.0, 101325.0])
    rho = np.array([1.225, 1.225])
    _, t0, p0 = total_state_from_static(p, rho, np.zeros_like(u), 1.4, 287.05)
    t_static = p / (rho * 287.05)
    if not np.allclose(t0, t_static) or not np.allclose(p0, p):
        print("[self-test] zero-velocity total state failed", file=sys.stderr)
        return 1
    print("[self-test] PASS")
    return 0


def mass_flux_diagnostics(mesh, cfg: dict, geometry: dict) -> dict:
    """Compute inlet/outlet mass flux and pressure ratio (trend indicators)."""
    dp = to_float(cfg, "geometry.global_resolution", 0.003)
    band = to_float(cfg, "boundary.boundary_n_layers", 3) * dp
    z_in = geometry["z_in"]
    z_out = geometry["z_out"]

    pts = np.asarray(mesh.points)
    z = pts[:, 2]
    rho = np.asarray(mesh.point_data["Density"])
    vel = np.asarray(mesh.point_data["Velocity"])
    p = np.asarray(mesh.point_data["Pressure"])
    vol = np.asarray(mesh.point_data["VolumetricMeasure"]) if "VolumetricMeasure" in mesh.point_data else np.full(len(rho), dp ** 3)

    w_z = vel[:, 2]
    A = vol / dp  # cubic-particle face area ~ dp^2

    near_in = np.abs(z - z_in) <= band
    near_out = np.abs(z - z_out) <= band
    m_in = float(np.sum(rho[near_in] * w_z[near_in] * A[near_in]))
    m_out = float(np.sum(rho[near_out] * w_z[near_out] * A[near_out]))
    denom = max(abs(m_in), abs(m_out))
    imbalance = abs(m_in - m_out) / denom if denom > 0 else float("nan")

    p_in = float(np.mean(p[near_in])) if near_in.any() else float("nan")
    p_out = float(np.mean(p[near_out])) if near_out.any() else float("nan")
    ratio = p_out / p_in if abs(p_in) > 0 else float("nan")

    return {
        "m_inlet": m_in,
        "m_outlet": m_out,
        "mass_flux_imbalance": imbalance,
        "p_inlet_avg": p_in,
        "p_outlet_avg": p_out,
        "pressure_ratio": ratio,
    }


def absolute_velocity_from_relative(points: np.ndarray, velocity_relative: np.ndarray,
                                    omega: float) -> np.ndarray:
    """Convert rotating-frame relative velocity w to absolute velocity u."""
    frame = np.column_stack((-omega * points[:, 1],
                             omega * points[:, 0],
                             np.zeros(len(points))))
    return velocity_relative + frame


def total_state_from_static(p: np.ndarray, rho: np.ndarray, velocity_absolute: np.ndarray,
                            gamma: float, gas_constant: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return static temperature, total temperature and total pressure."""
    rho_safe = np.maximum(rho, 1.0e-12)
    p_safe = np.maximum(p, 1.0e-12)
    temperature = p_safe / (rho_safe * gas_constant)
    cp = gamma * gas_constant / (gamma - 1.0)
    total_temperature = temperature + np.sum(velocity_absolute * velocity_absolute, axis=1) / (2.0 * cp)
    ratio = np.maximum(total_temperature / np.maximum(temperature, 1.0e-12), 1.0)
    total_pressure = p_safe * ratio ** (gamma / (gamma - 1.0))
    return temperature, total_temperature, total_pressure


def _mass_weighted_average(values: np.ndarray, weights: np.ndarray) -> float:
    mask = np.isfinite(values) & np.isfinite(weights) & (np.abs(weights) > 0.0)
    if not np.any(mask):
        return float("nan")
    w = np.abs(weights[mask])
    return float(np.sum(values[mask] * w) / np.sum(w))


def literature_metrics(mesh, cfg: dict, geometry: dict) -> dict:
    """Compute one-point Rotor67 literature-style metrics on physical faces."""
    dp = to_float(cfg, "geometry.global_resolution", 0.003)
    band = to_float(cfg, "boundary.boundary_n_layers", 3) * dp
    gamma = to_float(cfg, "physical.gamma", 1.4)
    gas_constant = to_float(cfg, "physical.gas_constant", 287.05)
    rpm = to_float(cfg, "rotor.rpm", 0.0)
    omega = to_float(cfg, "rotor.omega_rad_s", 2.0 * np.pi * rpm / 60.0)

    pts = np.asarray(mesh.points)
    z = pts[:, 2]
    rho = np.asarray(mesh.point_data["Density"])
    p = np.asarray(mesh.point_data["Pressure"])
    vel_rel = np.asarray(mesh.point_data["Velocity"])
    vol = np.asarray(mesh.point_data["VolumetricMeasure"]) if "VolumetricMeasure" in mesh.point_data else np.full(len(rho), dp ** 3)
    u_abs = absolute_velocity_from_relative(pts, vel_rel, omega)
    temperature, total_temperature, total_pressure = total_state_from_static(
        p, rho, u_abs, gamma, gas_constant)
    valid_sound = (p > 0.0) & (rho > 0.0)
    sound = np.full(len(rho), np.nan)
    sound[valid_sound] = np.sqrt(gamma * p[valid_sound] / rho[valid_sound])
    relative_mach = np.linalg.norm(vel_rel, axis=1) / sound

    area = vol / dp
    mass_weight = rho * vel_rel[:, 2] * area
    near_in = np.abs(z - geometry["z_in"]) <= band
    near_out = np.abs(z - geometry["z_out"]) <= band

    p0_in = _mass_weighted_average(total_pressure[near_in], mass_weight[near_in])
    p0_out = _mass_weighted_average(total_pressure[near_out], mass_weight[near_out])
    t0_in = _mass_weighted_average(total_temperature[near_in], mass_weight[near_in])
    t0_out = _mass_weighted_average(total_temperature[near_out], mass_weight[near_out])
    pr_total = p0_out / p0_in if np.isfinite(p0_in) and abs(p0_in) > 0 else float("nan")
    tr_total = t0_out / t0_in if np.isfinite(t0_in) and abs(t0_in) > 0 else float("nan")
    if np.isfinite(pr_total) and np.isfinite(tr_total) and abs(tr_total - 1.0) > 1.0e-12:
        efficiency = (pr_total ** ((gamma - 1.0) / gamma) - 1.0) / (tr_total - 1.0)
    else:
        efficiency = float("nan")

    return {
        "axis_mapping": "code z is Rotor67 axial direction; Velocity is rotor-relative w",
        "mass_flow_physical_inlet": float(np.sum(mass_weight[near_in])),
        "mass_flow_physical_outlet": float(np.sum(mass_weight[near_out])),
        "mass_averaged_total_pressure_in": p0_in,
        "mass_averaged_total_pressure_out": p0_out,
        "mass_averaged_total_temperature_in": t0_in,
        "mass_averaged_total_temperature_out": t0_out,
        "total_pressure_ratio": pr_total,
        "total_temperature_ratio": tr_total,
        "isentropic_efficiency": efficiency,
        "relative_mach_max": float(np.nanmax(relative_mach)) if np.any(np.isfinite(relative_mach)) else float("nan"),
        "relative_mach_mean": float(np.nanmean(relative_mach)) if np.any(np.isfinite(relative_mach)) else float("nan"),
        "physical_inlet_band_count": int(np.sum(near_in)),
        "physical_outlet_band_count": int(np.sum(near_out)),
    }


def claim_control(diag: dict, metrics: dict, geometry_audit: dict | None) -> dict:
    """Assign strict claim labels and missing gates."""
    labels = ["smoke-baseline"]
    missing = []
    if geometry_audit and geometry_audit.get("literature_geometry_status") == "literature-geometry-ready":
        labels.append("literature-geometry-ready")
    else:
        missing.append("geometry audit not literature-ready")
    pr = metrics.get("total_pressure_ratio", float("nan"))
    eff = metrics.get("isentropic_efficiency", float("nan"))
    mach_max = metrics.get("relative_mach_max", float("nan"))
    trend_ok = (
        np.isfinite(pr) and pr > 0.0 and
        np.isfinite(eff) and 0.0 <= eff <= 1.2 and
        np.isfinite(mach_max) and mach_max < 10.0)
    if trend_ok:
        labels.append("clean-laminar-trend")
    else:
        missing.append("physically bounded total metrics unavailable")
    missing.append("turbulence closure/y+ gate unavailable")
    missing.append("resolution independence gate unavailable")
    blocked = ["literature-performance-ready", "rans-comparable"]
    return {
        "claim_labels": labels,
        "blocked_claims": blocked,
        "missing_gates": missing,
        "literature_validation_pass": False,
        "rans_comparable": False,
    }


def load_geometry_audit(output_dir: Path) -> dict | None:
    path = output_dir / "rotor67_geometry_audit.json"
    if not path.exists():
        return None
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_validation_report(output_dir: Path, summary: dict) -> Path:
    report = output_dir / "literature_validation_report.md"
    claims = summary["claim_control"]
    metrics = summary["literature_metrics"]
    with open(report, "w", encoding="utf-8", newline="\n") as f:
        f.write("# Rotor67 Literature Validation Report\n\n")
        f.write("This report is generated from the current SPH smoke/trend workflow. It is not a RANS literature-performance validation.\n\n")
        f.write("## Claim Labels\n")
        for label in claims["claim_labels"]:
            f.write(f"- {label}\n")
        f.write("\n## Blocked Claims\n")
        for label in claims["blocked_claims"]:
            f.write(f"- {label}\n")
        f.write("\n## Missing Gates\n")
        for gate in claims["missing_gates"]:
            f.write(f"- {gate}\n")
        f.write("\n## Metrics\n")
        f.write(f"- total_pressure_ratio: {metrics.get('total_pressure_ratio')}\n")
        f.write(f"- total_temperature_ratio: {metrics.get('total_temperature_ratio')}\n")
        f.write(f"- isentropic_efficiency: {metrics.get('isentropic_efficiency')}\n")
        f.write(f"- relative_mach_max: {metrics.get('relative_mach_max')}\n")
    return report


def json_ready(value):
    if isinstance(value, dict):
        return {k: json_ready(v) for k, v in value.items()}
    if isinstance(value, list):
        return [json_ready(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def periodic_consistency(mesh, cfg: dict, geometry: dict) -> dict:
    """Compare periodic lower/upper bands after rotating vectors to a common side."""
    if "BoundaryKind" not in mesh.point_data:
        return {"available": False, "reason": "BoundaryKind missing"}
    pts = np.asarray(mesh.points)
    vel = np.asarray(mesh.point_data["Velocity"])
    mom = np.asarray(mesh.point_data["Momentum"]) if "Momentum" in mesh.point_data else None
    kind = np.asarray(mesh.point_data["BoundaryKind"])
    lower = kind == 3
    upper = kind == 4
    if not lower.any() or not upper.any():
        return {"available": False, "reason": "periodic band empty"}

    delta = geometry["periodic_delta_theta"]
    c, s = np.cos(delta), np.sin(delta)
    rot_plus = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    lower_vel_to_upper = vel[lower] @ rot_plus.T
    upper_vel = vel[upper]

    out = {
        "available": True,
        "lower_count": int(lower.sum()),
        "upper_count": int(upper.sum()),
        "velocity_mean_vector_diff_lower_to_upper": float(np.linalg.norm(
            lower_vel_to_upper.mean(axis=0) - upper_vel.mean(axis=0))),
        "pressure_mean_abs_diff": float(abs(np.mean(mesh.point_data["Pressure"][lower]) -
                                            np.mean(mesh.point_data["Pressure"][upper]))),
        "density_mean_abs_diff": float(abs(np.mean(mesh.point_data["Density"][lower]) -
                                           np.mean(mesh.point_data["Density"][upper]))),
    }
    if mom is not None:
        lower_mom_to_upper = mom[lower] @ rot_plus.T
        upper_mom = mom[upper]
        out["momentum_mean_vector_diff_lower_to_upper"] = float(np.linalg.norm(
            lower_mom_to_upper.mean(axis=0) - upper_mom.mean(axis=0)))
    return out


def wall_summary(output_dir: Path) -> dict:
    """Read the t=0 combined wall VTP and summarize BodyType counts."""
    wall_files = sorted(output_dir.glob("RotorWalls_*.vtp"))
    if not wall_files:
        return {"available": False, "reason": "RotorWalls_*.vtp missing"}
    try:
        import pyvista as pv
    except ImportError:
        return {"available": False, "reason": "pyvista missing", "frame": wall_files[0].name}
    wall_mesh = pv.read(str(wall_files[0]))
    if "BodyType" not in wall_mesh.point_data:
        return {"available": False, "reason": "BodyType missing", "frame": wall_files[0].name}
    body_type = np.asarray(wall_mesh.point_data["BodyType"]).astype(int)
    labels = {0: "blade", 1: "hub", 2: "shroud"}
    counts = {labels.get(i, f"type_{i}"): int(np.sum(body_type == i))
              for i in sorted(set(body_type.tolist()))}
    return {"available": True, "frame": wall_files[0].name, "counts": counts}


def meridional_slice(mesh, cfg: dict, output_dir: Path):
    """Plot meridional (r-z) pressure and Mach trend views."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pts = np.asarray(mesh.points)
    z = pts[:, 2]
    r = np.sqrt(pts[:, 0] ** 2 + pts[:, 1] ** 2)
    p = np.asarray(mesh.point_data["Pressure"])
    rho = np.asarray(mesh.point_data["Density"])
    vel = np.asarray(mesh.point_data["Velocity"])
    gamma = to_float(cfg, "physical.gamma", 1.4)
    rho_safe = np.maximum(rho, 1.0e-12)
    mach = np.linalg.norm(vel, axis=1) / np.sqrt(np.maximum(gamma * p / rho_safe, 1.0e-12))

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    sc0 = axes[0].scatter(z, r, c=p, s=2, cmap="jet")
    axes[0].set_xlabel("z (axial, m)")
    axes[0].set_ylabel("r (radial, m)")
    axes[0].set_title("Pressure")
    fig.colorbar(sc0, ax=axes[0], label="Pa")
    sc1 = axes[1].scatter(z, r, c=mach, s=2, cmap="viridis")
    axes[1].set_xlabel("z (axial, m)")
    axes[1].set_title("Mach")
    fig.colorbar(sc1, ax=axes[1], label="Mach")
    out = output_dir / "rotor67_meridional_pressure.png"
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def main():
    parser = argparse.ArgumentParser(description="Rotor 67 domain-STL post-processing")
    parser.add_argument("--case-dir", default=str(Path(__file__).parent),
                        help="Case directory containing config.ini")
    parser.add_argument("--output-dir", default=None,
                        help="Output directory with VTP files (default: case-dir/output)")
    parser.add_argument("--latest", action="store_true",
                        help="Read the latest frame (default behavior).")
    parser.add_argument("--self-test", action="store_true",
                        help="Run formula self-tests without reading VTP files.")
    args = parser.parse_args()

    if args.self_test:
        return run_self_test()

    case_dir = Path(args.case_dir)
    output_dir = Path(args.output_dir) if args.output_dir else case_dir / "output"
    cfg = parse_config(case_dir)
    geometry = load_domain_geometry(case_dir, cfg)

    print(f"[postprocess] case_dir={case_dir}")
    print(f"[postprocess] output_dir={output_dir}")
    print(f"[postprocess] geometry_source={geometry['source']} "
          f"z_in={geometry['z_in']} z_out={geometry['z_out']} "
          f"periodic_delta_theta={geometry['periodic_delta_theta']}")

    mesh, frame = load_last_frame(output_dir)
    print(f"[postprocess] loaded frame: {frame}, points={mesh.n_points}")

    diag = mass_flux_diagnostics(mesh, cfg, geometry)
    metrics = literature_metrics(mesh, cfg, geometry)
    periodic = periodic_consistency(mesh, cfg, geometry)
    walls = wall_summary(output_dir)
    geometry_audit = load_geometry_audit(output_dir)
    claims = claim_control(diag, metrics, geometry_audit)
    print("[postprocess] diagnostics (trend indicators):")
    for k, v in diag.items():
        print(f"    {k} = {v}")
    print("[postprocess] literature metrics (claim-gated):")
    for k in ("total_pressure_ratio", "total_temperature_ratio",
              "isentropic_efficiency", "relative_mach_max"):
        print(f"    {k} = {metrics[k]}")

    png = meridional_slice(mesh, cfg, output_dir)
    print(f"[postprocess] meridional pressure slice written to: {png}")

    # Trend-indicator verdicts; this smoke gate does not claim absolute accuracy.
    tol = to_float(cfg, "validation.mass_flux_imbalance_tol", 0.05)
    pr_min = to_float(cfg, "validation.pressure_ratio_min", 1.2)
    pr_max = to_float(cfg, "validation.pressure_ratio_max", 2.0)
    print("[postprocess] trend verdicts:")
    print(f"    mass_flux_imbalance < {tol}: "
          f"{'PASS' if diag['mass_flux_imbalance'] < tol else 'CHECK'} "
          f"({diag['mass_flux_imbalance']:.4f})")
    pr = diag["pressure_ratio"]
    print(f"    pressure_ratio in [{pr_min}, {pr_max}]: "
          f"{'PASS' if pr_min <= pr <= pr_max else 'CHECK'} ({pr:.4f})")

    summary = {
        "case": cfg.get("case.profile", "smoke") + "_" + cfg.get("inlet.case_type", "clean"),
        "config": {
            "case_profile": cfg.get("case.profile", "smoke"),
            "inlet_case_type": cfg.get("inlet.case_type", "clean"),
            "total_pressure_profile": cfg.get("inlet.total_pressure_profile", "uniform"),
            "distortion_intensity": to_float(cfg, "inlet.distortion_intensity", 0.0),
            "swirl_type": cfg.get("inlet.swirl_type", "none"),
            "swirl_angle_deg": to_float(cfg, "inlet.swirl_angle_deg", 0.0),
            "back_pressure": to_float(cfg, "outlet.back_pressure", 0.0),
        },
        "frame": frame,
        "geometry": geometry,
        "diagnostics": diag,
        "literature_metrics": metrics,
        "periodic_consistency": periodic,
        "wall_summary": walls,
        "geometry_audit": geometry_audit,
        "claim_control": claims,
        "trend_verdicts": {
            "mass_flux_imbalance": "PASS" if diag["mass_flux_imbalance"] < tol else "CHECK",
            "pressure_ratio": "PASS" if pr_min <= pr <= pr_max else "CHECK",
        },
    }
    json_path = output_dir / "rotor67_boundary_summary.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(json_ready(summary), f, indent=2, allow_nan=False)
    print(f"[postprocess] boundary summary written to: {json_path}")
    literature_json = output_dir / "rotor67_literature_summary.json"
    with open(literature_json, "w", encoding="utf-8") as f:
        json.dump(json_ready(summary), f, indent=2, allow_nan=False)
    report = write_validation_report(output_dir, summary)
    print(f"[postprocess] literature summary written to: {literature_json}")
    print(f"[postprocess] validation report written to: {report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
