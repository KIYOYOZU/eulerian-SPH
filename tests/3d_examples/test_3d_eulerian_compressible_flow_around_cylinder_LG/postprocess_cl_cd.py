#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute spanwise-normalized Cd/Cl for the 3D Eulerian fully compressible
cylinder case (Ma = 0.3, Re = 100).

Inputs
------
* ``output/Cylinder_TotalPressureForceFromFluid.dat``
* ``output/Cylinder_TotalViscousForceFromFluid.dat``

Outputs
-------
* ``results/cylinder_3d_cl_cd.csv`` with the columns
  ``time, fx_pressure, fy_pressure, fx_viscous, fy_viscous,
  fx_total, fy_total, cd, cl``.

The drag / lift coefficients use the spanwise-normalized reference
``0.5 * rho_inf * u_inf^2 * D * DW`` (D = cylinder diameter, DW = spanwise
spanwise width). By default all four values are derived from the case
``config.ini`` -- ``u_inf = mach_inf * c_inf`` -- so the compressible case can
never be normalised by a stale hard-coded ``u_ref``. The explicit flags remain
available as overrides.

The force files share the simulation time column, so rows are aligned
position-wise and a small mismatch (<= 1e-10 s) is silently snapped to the
minimum of the two.

Example
-------
::

    # Normalisation read from config.ini (recommended)
    python postprocess_cl_cd.py \
        --output-dir output \
        --config config.ini \
        --out results/cylinder_3d_cl_cd.csv
"""

from __future__ import annotations

import argparse
import configparser
import csv
import math
from pathlib import Path
from typing import Iterable, Optional


def _strip_token(token: str) -> str:
    return token.strip().strip('"')


def load_normalization(config_path: Path) -> dict[str, float]:
    """Derive rho_inf / u_inf / D / DW from the case ``config.ini``.

    The freestream speed is derived the same way the C++ config does
    (``u_inf = mach_inf * c_inf``) so the two can never disagree.
    """
    if not config_path.exists():
        raise FileNotFoundError(f"config.ini not found: {config_path}")
    parser = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    parser.read(config_path, encoding="utf-8")

    def need(section: str, key: str) -> float:
        if section not in parser or key not in parser[section]:
            raise KeyError(f"missing key '{key}' in [{section}] of {config_path}")
        return float(parser[section][key].strip())

    rho_inf = need("physical", "rho_inf")
    c_inf = need("physical", "c_inf")
    mach_inf = need("physical", "mach_inf")
    radius = need("geometry", "cylinder_radius")
    span_width = need("geometry", "dw")
    return {
        "rho0": rho_inf,
        "u_ref": mach_inf * c_inf,
        "diameter": 2.0 * radius,
        "span_width": span_width,
    }


def read_force_file(path: Path) -> tuple[list[str], list[list[float]]]:
    if not path.exists():
        raise FileNotFoundError(f"force file not found: {path}")
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"force file is empty: {path}")

    header = [_strip_token(token) for token in lines[0].split()]
    rows: list[list[float]] = []
    for line_no, line in enumerate(lines[1:], start=2):
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"{path}:{line_no}: expected time and 3 force components")
        try:
            rows.append([float(value) for value in parts[:4]])
        except ValueError as exc:
            raise ValueError(f"{path}:{line_no}: cannot parse numeric force row") from exc
    if not rows:
        raise ValueError(f"force file has header but no data rows: {path}")
    return header, rows


def combine_rows(pressure_rows: list[list[float]], viscous_rows: list[list[float]]) -> Iterable[tuple[float, float, float, float, float, float, float]]:
    n = min(len(pressure_rows), len(viscous_rows))
    if n == 0:
        raise ValueError("no overlapping force rows")
    for i in range(n):
        tp, fpx, fpy, fpz = pressure_rows[i]
        tv, fvx, fvy, fvz = viscous_rows[i]
        time = tp if abs(tp - tv) <= 1.0e-10 else min(tp, tv)
        yield time, fpx, fpy, fpz, fvx, fvy, fvz


def compute(output_dir: Path, rho0: float, u_ref: float, diameter: float, span_width: float, out_csv: Path) -> None:
    pressure_path = output_dir / "Cylinder_TotalPressureForceFromFluid.dat"
    viscous_path = output_dir / "Cylinder_TotalViscousForceFromFluid.dat"
    _, pressure_rows = read_force_file(pressure_path)
    _, viscous_rows = read_force_file(viscous_path)

    if min(rho0, u_ref, diameter, span_width) <= 0.0:
        raise ValueError(
            f"normalization inputs must all be positive, got rho0={rho0}, u_ref={u_ref}, "
            f"diameter={diameter}, span_width={span_width}"
        )
    scale = 0.5 * rho0 * u_ref * u_ref * diameter * span_width
    if not math.isfinite(scale) or scale <= 0.0:
        raise ValueError("invalid normalization scale; check rho0/u_ref/diameter/span_width")
    print(
        f"[cl_cd] normalisation: rho_inf={rho0} u_inf={u_ref} D={diameter} DW={span_width} "
        f"scale={scale:.9g}"
    )

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "time",
            "fx_pressure",
            "fy_pressure",
            "fx_viscous",
            "fy_viscous",
            "fx_total",
            "fy_total",
            "cd",
            "cl",
        ])
        for time, fpx, fpy, _fpz, fvx, fvy, _fvz in combine_rows(pressure_rows, viscous_rows):
            fx = fpx + fvx
            fy = fpy + fvy
            cd = fx / scale
            cl = fy / scale
            if not all(math.isfinite(value) for value in (fx, fy, cd, cl)):
                raise ValueError("non-finite Cd/Cl encountered")
            writer.writerow([time, fpx, fpy, fvx, fvy, fx, fy, cd, cl])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("output"), help="SPHinXsys output directory")
    parser.add_argument("--config", type=Path, default=Path("config.ini"),
                        help="case config.ini used to derive the normalisation (default: config.ini)")
    parser.add_argument("--rho0", type=float, default=None, help="override reference density rho_inf")
    parser.add_argument("--u-ref", type=float, default=None,
                        help="override freestream velocity u_inf (config gives mach_inf * c_inf)")
    parser.add_argument("--diameter", type=float, default=None, help="override cylinder diameter D")
    parser.add_argument("--span-width", type=float, default=None, help="override spanwise width DW")
    parser.add_argument("--out", type=Path, default=Path("results") / "cylinder_3d_cl_cd.csv", help="output CSV path")
    args = parser.parse_args()

    overrides = {
        "rho0": args.rho0,
        "u_ref": args.u_ref,
        "diameter": args.diameter,
        "span_width": args.span_width,
    }
    if all(value is not None for value in overrides.values()):
        norm = {key: float(value) for key, value in overrides.items()}
    else:
        norm = load_normalization(args.config)
        for key, value in overrides.items():
            if value is not None:
                norm[key] = float(value)

    compute(args.output_dir, norm["rho0"], norm["u_ref"], norm["diameter"], norm["span_width"], args.out)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
