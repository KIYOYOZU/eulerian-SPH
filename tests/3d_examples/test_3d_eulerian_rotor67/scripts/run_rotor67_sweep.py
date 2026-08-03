#!/usr/bin/env python3
"""Generate Rotor67 clean/distortion sweep configs and manifests.

Default behavior is dry-run generation only. Production-size local execution is
blocked by a particle-count preflight; HPC submission can consume the generated
manifest separately. The default matrix is intentionally compact: clean,
pressure distortion, swirl distortion, and one combined case.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from configparser import ConfigParser
from pathlib import Path


LOCAL_PARTICLE_LIMIT = 5_000_000


def read_ini(path: Path) -> ConfigParser:
    parser = ConfigParser()
    parser.optionxform = str
    with open(path, "r", encoding="utf-8") as f:
        parser.read_file(f)
    return parser


def write_ini(path: Path, parser: ConfigParser) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        parser.write(f, space_around_delimiters=True)


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def estimate_particles(parser: ConfigParser) -> int:
    rotor = parser["rotor"]
    geometry = parser["geometry"]
    hub = float(rotor.get("hub_radius", 0.089))
    shroud = float(rotor.get("shroud_radius", 0.257))
    n_blades = int(rotor.get("n_blades", 22))
    z_in = float(geometry.get("z_in", -0.04))
    z_out = float(geometry.get("z_out", 0.1303))
    dp = float(geometry.get("global_resolution", 0.002))
    sponge = float(geometry.get("sponge_width_factor", 5.0)) * dp
    pitch = 2.0 * math.pi / max(n_blades, 1)
    volume = 0.5 * pitch * (shroud * shroud - hub * hub) * ((z_out - z_in) + 2.0 * sponge)
    return int(max(volume / max(dp ** 3, 1.0e-18), 0.0))


def set_inlet_case(parser: ConfigParser, case: dict) -> None:
    inlet = parser["inlet"]
    inlet["case_type"] = str(case["case_type"])
    inlet["total_pressure_profile"] = str(case["total_pressure_profile"])
    inlet["distortion_intensity"] = f"{float(case['distortion_intensity']):.2f}"
    inlet["swirl_type"] = str(case["swirl_type"])
    inlet["swirl_angle_deg"] = f"{float(case['swirl_angle_deg']):.1f}"


def relocate_case_paths(parser: ConfigParser, base_path: Path, target_dir: Path) -> None:
    base_dir = base_path.parent
    for section, key in (("geometry", "domain_stl"),
                         ("meridional", "hub_curve"),
                         ("meridional", "shroud_curve")):
        raw = parser[section].get(key, "")
        path = Path(raw)
        if raw and not path.is_absolute():
            absolute = (base_dir / path).resolve()
            parser[section][key] = os.path.relpath(absolute, target_dir)


def load_matrix(path: Path) -> list[dict]:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    cases: list[dict] = [{"name": "clean", "case_type": "clean",
                          "total_pressure_profile": "uniform",
                          "distortion_intensity": 0.0,
                          "swirl_type": "none", "swirl_angle_deg": 0.0}]
    for key in ("pressure_cases", "swirl_cases", "combined_cases"):
        cases.extend(data.get(key, []))
    return cases


def main() -> int:
    ap = argparse.ArgumentParser(description="Rotor67 sweep config generator")
    ap.add_argument("--case-dir", default=str(Path(__file__).resolve().parents[1]))
    ap.add_argument("--base", default="cases/rotor67_literature/literature_clean.ini")
    ap.add_argument("--matrix", default="cases/rotor67_literature/distortion_matrix.json")
    ap.add_argument("--output", default="cases/rotor67_literature/generated")
    ap.add_argument("--back-pressure", nargs="+", type=float,
                    default=[130000.0])
    ap.add_argument("--dry-run", action="store_true", default=True)
    ap.add_argument("--allow-local-production", action="store_true")
    args = ap.parse_args()

    case_dir = Path(args.case_dir)
    base_path = case_dir / args.base
    matrix_path = case_dir / args.matrix
    out_dir = case_dir / args.output
    cases = load_matrix(matrix_path)
    manifest = {
        "mode": "dry-run" if args.dry_run else "run-requested",
        "case_root": ".",
        "base_config": str(Path(args.base)),
        "base_sha256": file_sha256(base_path),
        "local_particle_limit": LOCAL_PARTICLE_LIMIT,
        "points": [],
        "blocked": [],
    }

    for case in cases:
        for back_pressure in args.back_pressure:
            parser = read_ini(base_path)
            set_inlet_case(parser, case)
            parser["outlet"]["back_pressure"] = f"{back_pressure:.6f}"
            estimate = estimate_particles(parser)
            name = f"{case['name']}_bp{int(back_pressure)}"
            config_path = out_dir / f"{name}.ini"
            relocate_case_paths(parser, base_path, config_path.parent)
            write_ini(config_path, parser)
            point = {
                "name": name,
                "case": case["name"],
                "config": str(config_path.relative_to(case_dir)),
                "config_sha256": file_sha256(config_path),
                "back_pressure": back_pressure,
                "estimated_particles": estimate,
                "local_allowed": estimate <= LOCAL_PARTICLE_LIMIT,
                "claim": "smoke-baseline" if case["name"] == "clean" else "distortion-trend",
            }
            manifest["points"].append(point)
            if estimate > LOCAL_PARTICLE_LIMIT and not args.allow_local_production:
                manifest["blocked"].append({
                    "name": name,
                    "reason": "estimated particle count exceeds local safety limit",
                })

    manifest_path = out_dir / "rotor67_sweep_manifest.json"
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(manifest_path, "w", encoding="utf-8", newline="\n") as f:
        json.dump(manifest, f, indent=2)
    print(f"[sweep] generated {len(manifest['points'])} configs")
    print(f"[sweep] manifest: {manifest_path}")
    if manifest["blocked"]:
        print(f"[sweep] local production blocked for {len(manifest['blocked'])} points")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
