#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Post-process for the 2D Eulerian SPH shock-bubble interaction case (Kapila
five-equation model, stiffened gas EOS).

The Eulerian particles sit on a fixed lattice, so each particles_<step>.csv
(columns x, y, rho, p, u, v, alpha) is pivoted onto a 2D grid and rendered in
the reference window [0, L] x [0, H] (the extended right reservoir is cropped).

Outputs (in results/):
  * evolution_alpha.png    - montage of volume-fraction snapshots (bubble
                             deformation sequence),
  * evolution_density.png  - montage of density snapshots (shock + collapse),
  * evolution_pressure.png - montage of pressure snapshots,
  * evolution_alpha.gif    - animation of the alpha field (if Pillow exists),
  * evolution_summary.txt  - per-frame t, alpha bounds, mass drift, and the
                             early shock-front position/speed (expected
                             ~2788 m/s from the Rankine-Hugoniot jump).

Usage:
    python postprocess_multiphase_shock_bubble.py
"""

import os
import sys
import glob
import math

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# ---------------- reference setup (defaults mirror config.ini) -----------
L_WINDOW = 0.027          # reference window length (figure)
H_BOX = 0.012
BUBBLE_X, BUBBLE_Y, BUBBLE_R = 0.012, 0.006, 0.003
RHO_PS, P_PS, U_PS = 1323.65, 1.9e9, -681.58
RHO_0, P_0 = 1000.0, 1.0e5


def load_config():
    """Override the defaults above from the case config.ini (same keys as the
    C++ side) so editing the config cannot silently desync post-processing."""
    vals = {}
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config.ini")
    if os.path.exists(path):
        section = None
        with open(path, encoding="utf-8") as fh:
            for raw in fh:
                line = raw.split("#", 1)[0].split(";", 1)[0].strip()
                if not line:
                    continue
                if line.startswith("[") and line.endswith("]"):
                    section = line[1:-1].strip().lower()
                    continue
                if "=" in line and section:
                    key, val = (s.strip().lower() for s in line.split("=", 1))
                    vals[f"{section}.{key}"] = val

    def num(key, default):
        try:
            return float(vals[key])
        except (KeyError, ValueError):
            return default

    global L_WINDOW, H_BOX, BUBBLE_X, BUBBLE_Y, BUBBLE_R, RHO_PS, P_PS, U_PS, RHO_0, P_0
    L_WINDOW = num("geometry.l", L_WINDOW)
    H_BOX = num("geometry.h", H_BOX)
    BUBBLE_X = num("geometry.bubble_x", BUBBLE_X)
    BUBBLE_Y = num("geometry.bubble_y", BUBBLE_Y)
    BUBBLE_R = num("geometry.bubble_r", BUBBLE_R)
    RHO_PS = num("ic.rho_postshock", RHO_PS)
    P_PS = num("ic.p_postshock", P_PS)
    U_PS = num("ic.u_postshock", U_PS)
    RHO_0 = num("ic.rho_preshock", RHO_0)
    P_0 = num("ic.p_preshock", P_0)


def expected_shock_speed():
    """Left-running shock speed from the Rankine-Hugoniot jump conditions."""
    m = (P_PS - P_0) / (0.0 - U_PS)      # mass flux through the shock
    return m / RHO_0                      # |D|, m/s


def load_frame(csv_path):
    return np.loadtxt(csv_path, delimiter=",", skiprows=1)


def build_grid(data0):
    """Fixed-lattice pivot axes from the first frame."""
    ux = np.unique(data0[:, 0])
    uy = np.unique(data0[:, 1])
    return ux, uy


def pivot(data, ux, uy):
    ix = np.searchsorted(ux, data[:, 0])
    iy = np.searchsorted(uy, data[:, 1])
    nx, ny = ux.size, uy.size
    fields = {}
    for name, col in (("rho", 2), ("p", 3), ("u", 4), ("v", 5), ("alpha", 6)):
        f = np.full((ny, nx), np.nan)
        f[iy, ix] = data[:, col]
        fields[name] = f
    return fields


def shock_front_x(fields, ux, uy, dp):
    """Shock position on the centerline while it is right of the bubble.

    Scans the row nearest y = H/2 from just right of the bubble edge towards
    the window end; the front is the first x whose pressure exceeds the
    mid-jump value. Returns None once the shock has reached the bubble.
    """
    j = np.argmin(np.abs(uy - 0.5 * H_BOX))
    p_row = fields["p"][j, :]
    x_lo = BUBBLE_X + BUBBLE_R + 4.0 * dp
    p_mid = 0.5 * (P_0 + P_PS)
    sel = (ux > x_lo) & (ux < L_WINDOW)
    if not np.any(sel):
        return None
    xs_sel = ux[sel]
    p_sel = p_row[sel]
    hit = np.flatnonzero(p_sel > p_mid)
    if hit.size == 0:
        return None
    x_front = xs_sel[hit[0]]
    # only meaningful while the front is still right of the bubble
    return x_front if x_front > x_lo + 2.0 * dp else None


def montage(frames, times, key, cmap, vmin, vmax, log_scale, out_path, title):
    n = len(frames)
    ncol = min(3, n)
    nrow = int(math.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 3.4 * nrow),
                             squeeze=False)
    im = None
    for r in range(nrow):
        for cidx in range(ncol):
            ax = axes[r][cidx]
            k = r * ncol + cidx
            if k < n:
                f = frames[k][key]
                if log_scale:
                    f = np.log10(np.clip(f, 1.0, None))
                im = ax.pcolormesh(frames[k]["_ux"], frames[k]["_uy"], f,
                                   vmin=vmin, vmax=vmax, cmap=cmap,
                                   shading="nearest")
                ax.add_patch(plt.Circle((BUBBLE_X, BUBBLE_Y), BUBBLE_R,
                                        fill=False, ec="w", ls="--", lw=0.8))
                ax.set_title(f"t={times[k]:.2f} us", fontsize=9)
                ax.set_aspect("equal")
                ax.tick_params(labelsize=6)
            else:
                ax.axis("off")
    fig.subplots_adjust(right=0.90)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    fig.suptitle(title, fontsize=11)
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    print(f"Saved {out_path}")


def main():
    case_dir = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(case_dir, "output")
    res_dir = os.path.join(case_dir, "results")
    os.makedirs(res_dir, exist_ok=True)
    load_config()
    print(f"[config] window L={L_WINDOW}, H={H_BOX}, bubble=({BUBBLE_X},{BUBBLE_Y},R={BUBBLE_R}), "
          f"p_ps={P_PS:.3e}, p_0={P_0:.3e}")

    csvs = sorted(glob.glob(os.path.join(out_dir, "particles_*.csv")),
                  key=lambda p: int(os.path.basename(p)[len("particles_"):-len(".csv")]))
    if not csvs:
        print(f"No particles_*.csv found in {out_dir}")
        return 1

    dp = None
    frames, times, masses = [], [], []
    t_sec = []
    alpha_bounds = []
    shock_track = []  # (t, x_front)

    for csv_path in csvs:
        step = os.path.basename(csv_path)[len("particles_"):-len(".csv")]
        data = load_frame(csv_path)
        time_path = os.path.join(out_dir, f"time_{step}.txt")
        if os.path.exists(time_path):
            with open(time_path) as fh:
                t = float(fh.read().strip())
        else:
            t = 0.0

        if dp is None:
            ux0, uy0 = build_grid(data)
            dp = float(np.median(np.diff(ux0)))
            ux_w = ux0[ux0 <= L_WINDOW + 0.5 * dp]

        # crop to the reference window (drop the extended reservoir)
        keep = data[:, 0] <= L_WINDOW + 0.5 * dp
        data_w = data[keep]
        fields = pivot(data_w, ux_w, uy0)
        fields["_ux"] = ux_w
        fields["_uy"] = uy0
        t_sec.append(t)

        a = data_w[:, 6]
        alpha_bounds.append((float(a.min()), float(a.max())))
        # mass is conserved on the full domain (walls exchange no mass);
        # the windowed crop would show physical inflow/outflow instead.
        masses.append(float(data[:, 2].sum()) * dp * dp)
        x_front = shock_front_x(fields, ux_w, uy0, dp)
        if x_front is not None:
            shock_track.append((t, x_front))
        frames.append(fields)
        times.append(t * 1e6)  # microseconds for display

    n_frames = len(frames)
    print(f"Loaded {n_frames} frames, dp={dp:.3e}, t=[{times[0]:.3f}, {times[-1]:.3f}] us")

    # stale files from an older, longer run would make time non-monotone and
    # silently corrupt the shock-speed fit and the animation
    if any(b <= a for a, b in zip(t_sec, t_sec[1:])):
        print("ERROR: frame times not strictly increasing; stale particles_*.csv "
              "in output/. Clear output/ (or remove leftover files) and rerun.")
        return 2

    # ---------------- physical checks ----------------
    d_exp = expected_shock_speed()
    lines = []
    lines.append(f"shock-bubble interaction: {n_frames} frames, dp={dp:.3e}")
    lines.append(f"alpha bounds per frame: min={min(b[0] for b in alpha_bounds):.6f} "
                 f"max={max(b[1] for b in alpha_bounds):.6f}")
    m0 = masses[0]
    drift = max(abs(m - m0) / m0 for m in masses)
    lines.append(f"mass drift (full domain) max |dm|/m0 = {drift:.3e}")
    if len(shock_track) >= 2:
        ts = np.array([s[0] for s in shock_track])
        xs = np.array([s[1] for s in shock_track])
        # keep only the approach branch: segment speeds stay close to the
        # first segment's until impact; post-impact samples stall (the front
        # detector then tracks the bubble compression, not the shock).
        seg = -np.diff(xs) / np.diff(ts)
        n_seg = 0
        while n_seg < len(seg) and seg[n_seg] > 0.5 * seg[0]:
            n_seg += 1
        lines.append("shock front track (t [s], x [m]):")
        for tt_, xx_ in shock_track:
            lines.append(f"    {tt_:.6e}  {xx_:.6f}")
        if n_seg >= 1:
            speed = -float(np.polyfit(ts[: n_seg + 1], xs[: n_seg + 1], 1)[0])
            lines.append(f"measured |D| = {speed:.1f} m/s over t=[0, {ts[n_seg]:.3e}] s "
                         f"(expected {d_exp:.1f} m/s, RH)")
        else:
            lines.append("shock front track: no monotone approach branch")
    else:
        lines.append("shock front track: fewer than 2 samples (shock already at bubble)")
    summary = "\n".join(lines) + "\n"
    print(summary)
    with open(os.path.join(res_dir, "evolution_summary.txt"), "w") as f:
        f.write(summary)

    # ---------------- montages ----------------
    n_show = min(9, n_frames)
    idx = np.unique(np.linspace(0, n_frames - 1, n_show).astype(int))
    sel_frames = [frames[i] for i in idx]
    sel_times = [times[i] for i in idx]

    montage(sel_frames, sel_times, "alpha", "viridis", 0.0, 1.0, False,
            os.path.join(res_dir, "evolution_alpha.png"),
            "Shock-bubble interaction: volume fraction (gas)")
    montage(sel_frames, sel_times, "rho", "inferno", 0.0, math.log10(1500.0), True,
            os.path.join(res_dir, "evolution_density.png"),
            "Shock-bubble interaction: log10(density)")
    montage(sel_frames, sel_times, "p", "magma", 0.0, 2.0e9, False,
            os.path.join(res_dir, "evolution_pressure.png"),
            "Shock-bubble interaction: pressure")

    # ---------------- gif ----------------
    try:
        fig, ax = plt.subplots(figsize=(6.4, 3.6))
        im = ax.pcolormesh(frames[0]["_ux"], frames[0]["_uy"], frames[0]["alpha"],
                           vmin=0.0, vmax=1.0, cmap="viridis", shading="nearest")
        ax.add_patch(plt.Circle((BUBBLE_X, BUBBLE_Y), BUBBLE_R, fill=False,
                                ec="w", ls="--", lw=0.8))
        ax.set_aspect("equal")
        ttl = ax.set_title(f"t={times[0]:.2f} us")
        fig.colorbar(im, ax=ax, label="alpha (gas)")

        def animate(k):
            # fixed lattice: update the mappable in place instead of
            # rebuilding the QuadMesh (clearing collections breaks colorbar)
            im.set_array(np.asarray(frames[k]["alpha"], dtype=float).ravel())
            ttl.set_text(f"t={times[k]:.2f} us")
            return im

        anim = FuncAnimation(fig, animate, frames=n_frames)
        gif_path = os.path.join(res_dir, "evolution_alpha.gif")
        anim.save(gif_path, writer=PillowWriter(fps=max(2, n_frames // 6)))
        plt.close(fig)
        print(f"Saved {gif_path}")
    except Exception as exc:  # Pillow missing etc. - montages already saved
        print(f"[warn] gif skipped: {exc}")

    bounded = all(b[0] >= -1e-8 and b[1] <= 1 + 1e-8 for b in alpha_bounds)
    return 0 if bounded else 2


if __name__ == "__main__":
    sys.exit(main())
