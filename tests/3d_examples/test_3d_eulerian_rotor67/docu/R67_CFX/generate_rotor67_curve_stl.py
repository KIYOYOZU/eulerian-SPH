#!/usr/bin/env python3
"""Generate STL previews from the Rotor67 CFX-BladeGen curve export.

The CFX files use the Z axis as the rotor axis:
  - Rot_Hub.curve / Rot_Shd.curve: (theta, r, z), with theta currently 0.
  - Rot_Profile.curve: blade section samples as (r, circumferential_arc, z).

The profile points are mapped back to Cartesian coordinates with
theta = theta_blade + circumferential_arc / r.
"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


@dataclass
class Mesh:
    name: str
    vertices: np.ndarray
    faces: np.ndarray


PROFILE_HEADER = re.compile(r"^#\s*Profile\s+(\d+)\s+at\s+([0-9.]+)\s*%")


def read_inf(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("!"):
            continue
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        values[key.strip().lower()] = value.strip()
    return values


def read_curve_points(path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 3:
            raise ValueError(f"{path}: expected at least 3 columns, got {line!r}")
        rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
    if not rows:
        raise ValueError(f"{path}: no data rows found")
    return np.asarray(rows, dtype=float)


def read_profiles(path: Path) -> tuple[list[float], list[np.ndarray]]:
    spans: list[float] = []
    profiles: list[np.ndarray] = []
    current: list[list[float]] = []

    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        header = PROFILE_HEADER.match(line)
        if header:
            if current:
                profiles.append(np.asarray(current, dtype=float))
                current = []
            spans.append(float(header.group(2)))
            continue
        parts = line.split()
        if len(parts) < 3:
            raise ValueError(f"{path}: expected profile data row, got {line!r}")
        current.append([float(parts[0]), float(parts[1]), float(parts[2])])

    if current:
        profiles.append(np.asarray(current, dtype=float))
    if len(spans) != len(profiles):
        raise ValueError(
            f"{path}: found {len(spans)} headers but {len(profiles)} profiles"
        )
    return spans, profiles


def profile_to_cartesian(points: np.ndarray, theta_blade: float) -> np.ndarray:
    radius = points[:, 0]
    arc = points[:, 1]
    z = points[:, 2]
    if np.any(radius <= 0.0):
        raise ValueError("Profile radius must be positive")
    theta = theta_blade + arc / radius
    return np.column_stack((radius * np.cos(theta), radius * np.sin(theta), z))


def resample_closed_loop(points: np.ndarray, count: int) -> np.ndarray:
    loop = np.asarray(points, dtype=float)
    if np.linalg.norm(loop[0] - loop[-1]) > 1.0e-10:
        loop = np.vstack((loop, loop[0]))
    segment_lengths = np.linalg.norm(np.diff(loop, axis=0), axis=1)
    cumulative = np.r_[0.0, np.cumsum(segment_lengths)]
    total = cumulative[-1]
    if total <= 0.0:
        raise ValueError("Degenerate closed loop")
    samples = np.linspace(0.0, total, count + 1)[:-1]
    out = np.empty((count, 3), dtype=float)
    for axis in range(3):
        out[:, axis] = np.interp(samples, cumulative, loop[:, axis])
    return out


def cyclically_align_rings(rings: list[np.ndarray]) -> list[np.ndarray]:
    aligned = [rings[0]]
    for ring in rings[1:]:
        previous = aligned[-1]
        best_shift = 0
        best_score = math.inf
        for shift in range(len(ring)):
            rolled = np.roll(ring, -shift, axis=0)
            score = float(np.mean(np.linalg.norm(rolled - previous, axis=1)))
            if score < best_score:
                best_score = score
                best_shift = shift
        aligned.append(np.roll(ring, -best_shift, axis=0))
    return aligned


def signed_volume(vertices: np.ndarray, faces: np.ndarray) -> float:
    tri = vertices[faces]
    return float(np.einsum("ij,ij->i", tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum() / 6.0)


def build_blade_mesh(
    profiles: list[np.ndarray],
    theta_blade: float,
    samples_per_section: int,
) -> Mesh:
    rings = [
        resample_closed_loop(profile_to_cartesian(profile, theta_blade), samples_per_section)
        for profile in profiles
    ]
    rings = cyclically_align_rings(rings)

    section_count = len(rings)
    vertices = np.vstack(rings)
    faces: list[list[int]] = []
    n = samples_per_section

    for section in range(section_count - 1):
        base = section * n
        nxt = (section + 1) * n
        for i in range(n):
            a = base + i
            b = base + (i + 1) % n
            c = nxt + (i + 1) % n
            d = nxt + i
            faces.append([a, b, c])
            faces.append([a, c, d])

    root_center = vertices[:n].mean(axis=0)
    tip_start = (section_count - 1) * n
    tip_center = vertices[tip_start : tip_start + n].mean(axis=0)
    root_center_id = len(vertices)
    tip_center_id = len(vertices) + 1
    vertices = np.vstack((vertices, root_center, tip_center))

    for i in range(n):
        a = i
        b = (i + 1) % n
        faces.append([root_center_id, b, a])

    for i in range(n):
        a = tip_start + i
        b = tip_start + (i + 1) % n
        faces.append([tip_center_id, a, b])

    faces_array = np.asarray(faces, dtype=np.int64)
    if signed_volume(vertices, faces_array) < 0.0:
        faces_array = faces_array[:, [0, 2, 1]]

    return Mesh("rotor67_blade_profile", vertices, faces_array)


def blade_bounds(profiles: list[np.ndarray], theta_blade: float) -> dict[str, float]:
    data = np.vstack(profiles)
    radius = data[:, 0]
    theta = theta_blade + data[:, 1] / radius
    z = data[:, 2]
    return {
        "r_min": float(radius.min()),
        "r_max": float(radius.max()),
        "theta_min": float(theta.min()),
        "theta_max": float(theta.max()),
        "z_min": float(z.min()),
        "z_max": float(z.max()),
    }


def interp_radius_at_z(curve: np.ndarray, z: np.ndarray | float) -> np.ndarray | float:
    order = np.argsort(curve[:, 2])
    return np.interp(z, curve[order, 2], curve[order, 1])


def z_samples_for_window(
    z_min: float,
    z_max: float,
    segments: int,
) -> np.ndarray:
    if segments < 1:
        raise ValueError("Axial segments must be positive")
    return np.linspace(z_min, z_max, segments + 1)


def build_theta_z_surface(
    name: str,
    curve: np.ndarray,
    z_values: np.ndarray,
    theta_values: np.ndarray,
    flip: bool,
) -> Mesh:
    vertices = np.empty((len(z_values) * len(theta_values), 3), dtype=float)
    radius_values = interp_radius_at_z(curve, z_values)
    for j, (z, radius) in enumerate(zip(z_values, radius_values)):
        base = j * len(theta_values)
        vertices[base : base + len(theta_values), 0] = radius * np.cos(theta_values)
        vertices[base : base + len(theta_values), 1] = radius * np.sin(theta_values)
        vertices[base : base + len(theta_values), 2] = z

    faces: list[list[int]] = []
    nt = len(theta_values)
    for j in range(len(z_values) - 1):
        base = j * nt
        nxt = (j + 1) * nt
        for i in range(nt - 1):
            a = base + i
            b = base + i + 1
            c = nxt + i + 1
            d = nxt + i
            if flip:
                faces.append([a, c, b])
                faces.append([a, d, c])
            else:
                faces.append([a, b, c])
                faces.append([a, c, d])
    return Mesh(name, vertices, np.asarray(faces, dtype=np.int64))


def build_z_plane(
    name: str,
    hub: np.ndarray,
    shroud: np.ndarray,
    z: float,
    theta_values: np.ndarray,
    radial_segments: int,
    flip: bool,
) -> Mesh:
    r_hub = float(interp_radius_at_z(hub, z))
    r_shroud = float(interp_radius_at_z(shroud, z))
    radius_values = np.linspace(r_hub, r_shroud, radial_segments + 1)
    vertices = np.empty((len(radius_values) * len(theta_values), 3), dtype=float)

    for j, radius in enumerate(radius_values):
        base = j * len(theta_values)
        vertices[base : base + len(theta_values), 0] = radius * np.cos(theta_values)
        vertices[base : base + len(theta_values), 1] = radius * np.sin(theta_values)
        vertices[base : base + len(theta_values), 2] = z

    faces: list[list[int]] = []
    nt = len(theta_values)
    for j in range(len(radius_values) - 1):
        base = j * nt
        nxt = (j + 1) * nt
        for i in range(nt - 1):
            a = base + i
            b = base + i + 1
            c = nxt + i + 1
            d = nxt + i
            if flip:
                faces.append([a, c, b])
                faces.append([a, d, c])
            else:
                faces.append([a, b, c])
                faces.append([a, c, d])
    return Mesh(name, vertices, np.asarray(faces, dtype=np.int64))


def build_theta_plane(
    name: str,
    hub: np.ndarray,
    shroud: np.ndarray,
    z_values: np.ndarray,
    theta: float,
    radial_segments: int,
    flip: bool,
) -> Mesh:
    vertices = np.empty((len(z_values) * (radial_segments + 1), 3), dtype=float)
    for j, z in enumerate(z_values):
        r_hub = float(interp_radius_at_z(hub, z))
        r_shroud = float(interp_radius_at_z(shroud, z))
        radius_values = np.linspace(r_hub, r_shroud, radial_segments + 1)
        base = j * (radial_segments + 1)
        vertices[base : base + radial_segments + 1, 0] = radius_values * math.cos(theta)
        vertices[base : base + radial_segments + 1, 1] = radius_values * math.sin(theta)
        vertices[base : base + radial_segments + 1, 2] = z

    faces: list[list[int]] = []
    nr = radial_segments + 1
    for j in range(len(z_values) - 1):
        base = j * nr
        nxt = (j + 1) * nr
        for i in range(radial_segments):
            a = base + i
            b = base + i + 1
            c = nxt + i + 1
            d = nxt + i
            if flip:
                faces.append([a, c, b])
                faces.append([a, d, c])
            else:
                faces.append([a, b, c])
                faces.append([a, c, d])
    return Mesh(name, vertices, np.asarray(faces, dtype=np.int64))


def build_local_domain_meshes(
    hub: np.ndarray,
    shroud: np.ndarray,
    profiles: list[np.ndarray],
    theta_blade: float,
    profile_samples: int,
    theta_segments: int,
    axial_segments: int,
    radial_segments: int,
    axial_margin: float,
    theta_margin: float,
) -> tuple[list[Mesh], dict[str, float]]:
    bounds = blade_bounds(profiles, theta_blade)
    hub_z_min = float(hub[:, 2].min())
    hub_z_max = float(hub[:, 2].max())
    shroud_z_min = float(shroud[:, 2].min())
    shroud_z_max = float(shroud[:, 2].max())
    z_min = max(bounds["z_min"] - axial_margin, hub_z_min, shroud_z_min)
    z_max = min(bounds["z_max"] + axial_margin, hub_z_max, shroud_z_max)
    if z_min >= z_max:
        raise ValueError("Invalid local domain z window")

    theta_min = bounds["theta_min"] - theta_margin
    theta_max = bounds["theta_max"] + theta_margin
    if theta_min >= theta_max:
        raise ValueError("Invalid local domain theta window")

    z_values = z_samples_for_window(z_min, z_max, axial_segments)
    theta_values = np.linspace(theta_min, theta_max, theta_segments + 1)

    meshes = [
        build_theta_z_surface("rotor67_domain_hub", hub, z_values, theta_values, flip=True),
        build_theta_z_surface("rotor67_domain_shroud", shroud, z_values, theta_values, flip=False),
        build_z_plane(
            "rotor67_domain_inlet",
            hub,
            shroud,
            z_min,
            theta_values,
            radial_segments,
            flip=True,
        ),
        build_z_plane(
            "rotor67_domain_outlet",
            hub,
            shroud,
            z_max,
            theta_values,
            radial_segments,
            flip=False,
        ),
        build_theta_plane(
            "rotor67_domain_periodic_min",
            hub,
            shroud,
            z_values,
            theta_min,
            radial_segments,
            flip=True,
        ),
        build_theta_plane(
            "rotor67_domain_periodic_max",
            hub,
            shroud,
            z_values,
            theta_max,
            radial_segments,
            flip=False,
        ),
        build_blade_mesh(profiles, theta_blade, profile_samples),
    ]
    window = {
        **bounds,
        "domain_z_min": z_min,
        "domain_z_max": z_max,
        "domain_theta_min": theta_min,
        "domain_theta_max": theta_max,
    }
    return meshes, window


def stl_normal(triangle: np.ndarray) -> np.ndarray:
    normal = np.cross(triangle[1] - triangle[0], triangle[2] - triangle[0])
    length = float(np.linalg.norm(normal))
    if length <= 0.0:
        return np.zeros(3)
    return normal / length


def write_ascii_solid(out, mesh: Mesh) -> None:
    out.write(f"solid {mesh.name}\n")
    for face in mesh.faces:
        triangle = mesh.vertices[face]
        normal = stl_normal(triangle)
        out.write(f"  facet normal {normal[0]:.9e} {normal[1]:.9e} {normal[2]:.9e}\n")
        out.write("    outer loop\n")
        for vertex in triangle:
            out.write(f"      vertex {vertex[0]:.9e} {vertex[1]:.9e} {vertex[2]:.9e}\n")
        out.write("    endloop\n")
        out.write("  endfacet\n")
    out.write(f"endsolid {mesh.name}\n")


def write_multi_solid_stl(path: Path, meshes: list[Mesh]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as out:
        for mesh in meshes:
            write_ascii_solid(out, mesh)


def edge_report(mesh: Mesh) -> tuple[dict[int, int], int]:
    counts: dict[tuple[int, int], int] = {}
    for face in mesh.faces:
        for i in range(3):
            a = int(face[i])
            b = int(face[(i + 1) % 3])
            edge = (min(a, b), max(a, b))
            counts[edge] = counts.get(edge, 0) + 1
    histogram: dict[int, int] = {}
    for count in counts.values():
        histogram[count] = histogram.get(count, 0) + 1
    boundary_edges = sum(1 for count in counts.values() if count == 1)
    return histogram, boundary_edges


def write_domain_preview(
    path: Path,
    meshes: list[Mesh],
    window: dict[str, float],
) -> None:
    fig = plt.figure(figsize=(10, 8), constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    colors = {
        "rotor67_domain_hub": "#4477aa",
        "rotor67_domain_shroud": "#cc6677",
        "rotor67_domain_inlet": "#66ccee",
        "rotor67_domain_outlet": "#aa3377",
        "rotor67_domain_periodic_min": "#228833",
        "rotor67_domain_periodic_max": "#ee7733",
        "rotor67_blade_profile": "#222222",
    }
    for mesh in meshes:
        face_stride = max(1, len(mesh.faces) // 1200)
        tris = mesh.vertices[mesh.faces[::face_stride]]
        ax.plot_trisurf(
            tris[:, :, 0].ravel(),
            tris[:, :, 1].ravel(),
            tris[:, :, 2].ravel(),
            triangles=np.arange(tris.size // 3).reshape(-1, 3),
            color=colors.get(mesh.name, "#999999"),
            alpha=0.25 if mesh.name != "rotor67_blade_profile" else 0.65,
            linewidth=0.15,
            edgecolor="#333333",
        )
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_zlabel("z (m)")
    ax.set_title(
        "Rotor67 local computational domain\n"
        f"z=[{window['domain_z_min']:.3f}, {window['domain_z_max']:.3f}], "
        f"theta=[{math.degrees(window['domain_theta_min']):.1f}deg, "
        f"{math.degrees(window['domain_theta_max']):.1f}deg]"
    )
    all_points = np.vstack([mesh.vertices for mesh in meshes])
    center = all_points.mean(axis=0)
    radius = float((all_points.max(axis=0) - all_points.min(axis=0)).max() * 0.55)
    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)
    ax.view_init(elev=23, azim=-58)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--theta-blade", type=float, default=0.0, help="Blade reference angle in radians")
    parser.add_argument("--profile-samples", type=int, default=192, help="Resampled points per blade profile")
    parser.add_argument("--theta-segments", type=int, default=72, help="Angular segments for local domain walls")
    parser.add_argument("--axial-segments", type=int, default=96, help="Axial segments for local domain walls")
    parser.add_argument("--radial-segments", type=int, default=24, help="Radial segments for inlet/outlet/periodic faces")
    parser.add_argument("--axial-margin", type=float, default=0.04, help="Axial margin around blade profile z-range")
    parser.add_argument("--theta-margin", type=float, default=0.04, help="Theta margin around blade profile theta-range")
    parser.add_argument("--output-prefix", default="rotor67_domain", help="Output file prefix")
    args = parser.parse_args()

    base = Path(__file__).resolve().parent
    inf = read_inf(base / "Rot.inf")
    axis = inf.get("axis of rotation", "Z").upper()
    if axis != "Z":
        raise ValueError(f"Only Z-axis Rotor67 export is supported, got axis {axis!r}")

    hub_file = inf.get("hub data file", "Rot_Hub.curve")
    shroud_file = inf.get("shroud data file", "Rot_Shd.curve")
    profile_file = inf.get("profile data file", "Rot_Profile.curve")
    blade_sets = inf.get("number of blade sets", "unknown")
    units = inf.get("geometry units", "unknown")

    hub = read_curve_points(base / hub_file)
    shroud = read_curve_points(base / shroud_file)
    spans, profiles = read_profiles(base / profile_file)

    meshes, window = build_local_domain_meshes(
        hub,
        shroud,
        profiles,
        args.theta_blade,
        args.profile_samples,
        args.theta_segments,
        args.axial_segments,
        args.radial_segments,
        args.axial_margin,
        args.theta_margin,
    )
    blade_mesh = meshes[-1]
    output_stl = base / f"{args.output_prefix}.stl"
    write_multi_solid_stl(output_stl, meshes)

    preview = base / f"{args.output_prefix}_preview.png"
    write_domain_preview(preview, meshes, window)

    blade_edges, blade_boundary_edges = edge_report(blade_mesh)
    print("Input summary:")
    print(f"  axis: {axis}, blade sets: {blade_sets}, units: {units}")
    print(f"  hub points: {len(hub)}")
    print(f"  shroud points: {len(shroud)}")
    print(f"  profiles: {len(profiles)}, spans: {spans[0]:.1f}%..{spans[-1]:.1f}%")
    print("Blade STL check:")
    print(f"  vertices: {len(blade_mesh.vertices)}, faces: {len(blade_mesh.faces)}")
    print(f"  signed volume: {signed_volume(blade_mesh.vertices, blade_mesh.faces):.9e} m^3")
    print(f"  edge count histogram: {blade_edges}")
    print(f"  boundary edges: {blade_boundary_edges}")
    print("Local domain window:")
    print(f"  z: [{window['domain_z_min']:.9e}, {window['domain_z_max']:.9e}] m")
    print(
        "  theta: "
        f"[{window['domain_theta_min']:.9e}, {window['domain_theta_max']:.9e}] rad "
        f"([{math.degrees(window['domain_theta_min']):.3f}, "
        f"{math.degrees(window['domain_theta_max']):.3f}] deg)"
    )
    print("Solid facet counts:")
    for mesh in meshes:
        print(f"  {mesh.name}: {len(mesh.faces)}")
    print("Wrote:")
    print(f"  {output_stl}")
    print(f"  {preview}")


if __name__ == "__main__":
    main()
