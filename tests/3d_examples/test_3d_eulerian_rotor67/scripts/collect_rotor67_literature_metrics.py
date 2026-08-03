#!/usr/bin/env python3
"""Aggregate Rotor67 literature summary JSON files into CSV."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def load_summary(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    metrics = data.get("literature_metrics", {})
    diag = data.get("diagnostics", {})
    claims = data.get("claim_control", {})
    config = data.get("config", {})
    return {
        "summary": str(path),
        "case": data.get("case", path.parent.name),
        "case_profile": config.get("case_profile"),
        "inlet_case_type": config.get("inlet_case_type"),
        "total_pressure_profile": config.get("total_pressure_profile"),
        "distortion_intensity": config.get("distortion_intensity"),
        "swirl_type": config.get("swirl_type"),
        "swirl_angle_deg": config.get("swirl_angle_deg"),
        "back_pressure": config.get("back_pressure"),
        "mass_flow_outlet": metrics.get("mass_flow_physical_outlet"),
        "total_pressure_ratio": metrics.get("total_pressure_ratio"),
        "total_temperature_ratio": metrics.get("total_temperature_ratio"),
        "isentropic_efficiency": metrics.get("isentropic_efficiency"),
        "relative_mach_max": metrics.get("relative_mach_max"),
        "static_pressure_ratio_trend": diag.get("pressure_ratio"),
        "mass_flux_imbalance": diag.get("mass_flux_imbalance"),
        "claim_labels": "|".join(claims.get("claim_labels", [])),
        "blocked_claims": "|".join(claims.get("blocked_claims", [])),
        "missing_gates": "|".join(claims.get("missing_gates", [])),
        "literature_validation_pass": claims.get("literature_validation_pass"),
        "rans_comparable": claims.get("rans_comparable"),
    }


def add_clean_deltas(rows: list[dict]) -> None:
    clean = next((r for r in rows if r.get("inlet_case_type") == "clean"), None)
    if clean is None:
        for row in rows:
            row["delta_mass_flow_vs_clean"] = ""
            row["delta_total_pressure_ratio_vs_clean"] = ""
            row["delta_stability_range_vs_clean"] = "unavailable"
        return
    try:
        clean_m = float(clean["mass_flow_outlet"])
        clean_pr = float(clean["total_pressure_ratio"])
    except (TypeError, ValueError):
        clean_m = clean_pr = None
    for row in rows:
        try:
            row["delta_mass_flow_vs_clean"] = (
                float(row["mass_flow_outlet"]) - clean_m if clean_m is not None else "")
            row["delta_total_pressure_ratio_vs_clean"] = (
                float(row["total_pressure_ratio"]) - clean_pr if clean_pr is not None else "")
        except (TypeError, ValueError):
            row["delta_mass_flow_vs_clean"] = ""
            row["delta_total_pressure_ratio_vs_clean"] = ""
        row["delta_stability_range_vs_clean"] = "unavailable"


def main() -> int:
    ap = argparse.ArgumentParser(description="Collect Rotor67 literature metrics")
    ap.add_argument("summaries", nargs="*", help="rotor67_literature_summary.json files")
    ap.add_argument("--glob", default=None, help="Glob pattern for summaries")
    ap.add_argument("--output", default="output/rotor67_literature_metrics.csv")
    args = ap.parse_args()

    paths = [Path(p) for p in args.summaries]
    if args.glob:
        paths.extend(Path().glob(args.glob))
    if not paths:
        default = Path("output/rotor67_literature_summary.json")
        if default.exists():
            paths.append(default)
    rows = [load_summary(p) for p in paths if p.exists()]
    add_clean_deltas(rows)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "summary", "case", "case_profile", "inlet_case_type",
        "total_pressure_profile", "distortion_intensity", "swirl_type",
        "swirl_angle_deg", "back_pressure", "mass_flow_outlet", "total_pressure_ratio",
        "total_temperature_ratio", "isentropic_efficiency", "relative_mach_max",
        "static_pressure_ratio_trend", "mass_flux_imbalance",
        "delta_mass_flow_vs_clean", "delta_total_pressure_ratio_vs_clean",
        "delta_stability_range_vs_clean", "claim_labels", "blocked_claims",
        "missing_gates", "literature_validation_pass", "rans_comparable",
    ]
    with open(out, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"[collect] wrote {len(rows)} rows to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
