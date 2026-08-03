#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Post-process for the 2D Eulerian SPH shock tube (Lax problem).

Reads particles_<step>.csv dumped by the C++ probe, bins particles along x
to get the numerical profile, overlays the exact Riemann solution (Toro,
Riemann Solvers and Numerical Methods for Fluid Dynamics, Ch. 4), and
reports L1 errors.

Usage:
    python postprocess_shock_tube.py [step]
    step defaults to the latest particles_*.csv found in ./output.
"""

import os
import sys
import glob
import math

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------
# Problem parameters (must match shock_tube.h)
# ---------------------------------------------------------------------
GAMMA = 1.4
RHO_L, P_L, U_L = 1.0, 0.425, 0.0
RHO_R, P_R, U_R = 0.125, 0.1, 0.0
X_MEMBRANE = 0.5
DOMAIN_X0, DOMAIN_X1 = 0.0, 1.0


def sound_speed(p, rho, gamma=GAMMA):
    return math.sqrt(gamma * p / rho)


def _pressure_functions(p, dk, rho_k, p_k, gamma):
    """Toro eq. (4.5)/(4.6): f_L, f_R and their derivatives."""
    if p <= 0.0:
        # Non-physical; return large value to push Newton iteration back.
        return 1.0e18, 0.0
    c_k = sound_speed(p_k, rho_k, gamma)
    a_k = 2.0 / ((gamma + 1.0) * rho_k)
    b_k = (gamma - 1.0) / (gamma + 1.0) * p_k
    if p <= p_k:  # rarefaction (rarefaction fan)
        f = 2.0 * c_k / (gamma - 1.0) * ((p / p_k) ** ((gamma - 1.0) / (2.0 * gamma)) - 1.0)
        df = (1.0 / (rho_k * c_k)) * (p / p_k) ** (-(gamma + 1.0) / (2.0 * gamma))
    else:  # shock
        sq = math.sqrt(a_k / (b_k + p))
        f = (p - p_k) * sq
        df = sq * (1.0 - 0.5 * (p - p_k) / (b_k + p))
    return f, df


def _solve_star_pressure(p_guess=0.5 * (P_L + P_R), tol=1.0e-11, max_iter=50):
    """Newton-Raphson for the star-region pressure (Toro Sec. 4.3)."""
    p = max(p_guess, 1.0e-9)
    for _ in range(max_iter):
        fL, dfL = _pressure_functions(p, "L", RHO_L, P_L, GAMMA)
        fR, dfR = _pressure_functions(p, "R", RHO_R, P_R, GAMMA)
        f = fL + fR + (U_R - U_L)
        df = dfL + dfR
        dp = -f / df if df != 0.0 else -f
        p_new = p + dp
        if p_new <= 0.0:
            p_new = 1.0e-9
        if abs(p_new - p) / (0.5 * (p_new + p)) < tol:
            return p_new
        p = p_new
    return p


def _star_velocity(p_star):
    fL, _ = _pressure_functions(p_star, "L", RHO_L, P_L, GAMMA)
    fR, _ = _pressure_functions(p_star, "R", RHO_R, P_R, GAMMA)
    return 0.5 * (U_L + U_R + fR - fL)


def _rho_star_shock(p_star, rho_k, p_k):
    """Post-shock density (Rankine-Hugoniot), Toro eq. (4.50)."""
    return rho_k * ((p_star / p_k + (GAMMA - 1.0) / (GAMMA + 1.0)) /
                    ((GAMMA - 1.0) / (GAMMA + 1.0) * p_star / p_k + 1.0))


def _rho_star_rarefaction(p_star, rho_k, p_k):
    """Isentropic density behind rarefaction, Toro eq. (4.47)."""
    return rho_k * (p_star / p_k) ** (1.0 / GAMMA)


def exact_solution(x, t, x0=X_MEMBRANE):
    """Return (rho, p, u) arrays along x at time t."""
    x = np.asarray(x, dtype=float)
    rho = np.empty_like(x)
    p = np.empty_like(x)
    u = np.empty_like(x)

    if t <= 0.0:
        left = x < x0
        rho[left] = RHO_L; p[left] = P_L; u[left] = U_L
        rho[~left] = RHO_R; p[~left] = P_R; u[~left] = U_R
        return rho, p, u

    p_star = _solve_star_pressure()
    u_star = _star_velocity(p_star)
    cL = sound_speed(P_L, RHO_L)
    cR = sound_speed(P_R, RHO_R)

    # Toro shock speed (eq. 4.50/4.52): S = u_k +/- c_k * sqrt((gamma+1)/(2gamma) * (p*/p_k - 1) + 1)
    def shock_speed(u_k, c_k, p_k, p_star, sign):
        return u_k + sign * c_k * math.sqrt((GAMMA + 1.0) / (2.0 * GAMMA) * (p_star / p_k - 1.0) + 1.0)

    # Left wave
    if p_star <= P_L:  # left rarefaction (left-going)
        rho_star_L = _rho_star_rarefaction(p_star, RHO_L, P_L)
        c_star_L = cL * (p_star / P_L) ** ((GAMMA - 1.0) / (2.0 * GAMMA))
        S_HL = x0 + (U_L - cL) * t          # fan head (left-going, fastest)
        S_TL = x0 + (u_star - c_star_L) * t  # fan tail (left-going, slowest)
    else:  # left shock (left-going)
        rho_star_L = _rho_star_shock(p_star, RHO_L, P_L)
        S_L = x0 + shock_speed(U_L, cL, P_L, p_star, -1.0) * t
        S_HL = S_TL = None  # no fan

    # Right wave
    if p_star <= P_R:  # right rarefaction (right-going)
        rho_star_R = _rho_star_rarefaction(p_star, RHO_R, P_R)
        c_star_R = cR * (p_star / P_R) ** ((GAMMA - 1.0) / (2.0 * GAMMA))
        S_HR = x0 + (U_R + cR) * t          # fan head (right-going, fastest)
        S_TR = x0 + (u_star + c_star_R) * t  # fan tail (right-going, slowest)
    else:  # right shock (right-going)
        rho_star_R = _rho_star_shock(p_star, RHO_R, P_R)
        S_R = x0 + shock_speed(U_R, cR, P_R, p_star, +1.0) * t
        S_HR = S_TR = None

    # Assemble by region
    for i, xi in enumerate(x):
        s = (xi - x0) / t  # similarity variable
        if s < u_star:  # left of contact
            if p_star <= P_L:  # left rarefaction
                if xi < S_HL:  # left state (pre-fan)
                    rho[i], p[i], u[i] = RHO_L, P_L, U_L
                elif xi > S_TL:  # post-fan star state
                    rho[i], p[i], u[i] = rho_star_L, p_star, u_star
                else:  # inside fan (isentropic)
                    f = 2.0 / (GAMMA + 1.0) + (GAMMA - 1.0) / ((GAMMA + 1.0) * cL) * (U_L - s)
                    p_f = P_L * f ** (2.0 * GAMMA / (GAMMA - 1.0))
                    rho_f = RHO_L * f ** (2.0 / (GAMMA - 1.0))
                    u_f = 2.0 / (GAMMA + 1.0) * (cL + (GAMMA - 1.0) / 2.0 * U_L + s)
                    rho[i], p[i], u[i] = rho_f, p_f, u_f
            else:  # left shock: left state vs star
                if xi < S_L:
                    rho[i], p[i], u[i] = RHO_L, P_L, U_L
                else:
                    rho[i], p[i], u[i] = rho_star_L, p_star, u_star
        else:  # right of contact
            if p_star <= P_R:  # right rarefaction
                if xi > S_HR:
                    rho[i], p[i], u[i] = RHO_R, P_R, U_R
                elif xi < S_TR:
                    rho[i], p[i], u[i] = rho_star_R, p_star, u_star
                else:
                    f = 2.0 / (GAMMA + 1.0) - (GAMMA - 1.0) / ((GAMMA + 1.0) * cR) * (U_R - s)
                    p_f = P_R * f ** (2.0 * GAMMA / (GAMMA - 1.0))
                    rho_f = RHO_R * f ** (2.0 / (GAMMA - 1.0))
                    u_f = 2.0 / (GAMMA + 1.0) * (-cR + (GAMMA - 1.0) / 2.0 * U_R + s)
                    rho[i], p[i], u[i] = rho_f, p_f, u_f
            else:  # right shock
                if xi > S_R:
                    rho[i], p[i], u[i] = RHO_R, P_R, U_R
                else:
                    rho[i], p[i], u[i] = rho_star_R, p_star, u_star
    return rho, p, u


# ---------------------------------------------------------------------
# Numerical profile: bin particles along x
# ---------------------------------------------------------------------
def load_particles(csv_path):
    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    return data  # columns: x, y, rho, p, u


def bin_profile(data, n_bins=200):
    x = data[:, 0]
    x_edges = np.linspace(DOMAIN_X0, DOMAIN_X1, n_bins + 1)
    centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    rho = np.zeros(n_bins)
    p = np.zeros(n_bins)
    u = np.zeros(n_bins)
    cnt = np.zeros(n_bins, dtype=int)
    for row in data:
        xi = row[0]
        idx = int((xi - DOMAIN_X0) / (DOMAIN_X1 - DOMAIN_X0) * n_bins)
        idx = min(max(idx, 0), n_bins - 1)
        rho[idx] += row[2]
        p[idx] += row[3]
        u[idx] += row[4]
        cnt[idx] += 1
    mask = cnt > 0
    rho[mask] /= cnt[mask]
    p[mask] /= cnt[mask]
    u[mask] /= cnt[mask]
    return centers, rho, p, u, cnt


def main():
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
    res_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(res_dir, exist_ok=True)
    if len(sys.argv) > 1:
        step = sys.argv[1]
        csv_path = os.path.join(out_dir, f"particles_{step}.csv")
    else:
        csvs = sorted(glob.glob(os.path.join(out_dir, "particles_*.csv")),
                      key=lambda p: int(os.path.basename(p)[len("particles_"):-len(".csv")]))
        if not csvs:
            print(f"No particles_*.csv found in {out_dir}")
            sys.exit(1)
        csv_path = csvs[-1]
        step = os.path.basename(csv_path)[len("particles_"):-len(".csv")]

    # Physical time at this step
    time_path = os.path.join(out_dir, f"time_{step}.txt")
    t = float(open(time_path).read().strip()) if os.path.exists(time_path) else 0.2

    data = load_particles(csv_path)
    xc, rho_n, p_n, u_n, cnt = bin_profile(data, n_bins=400)

    rho_e, p_e, u_e = exact_solution(xc, t)

    # L1 error (only over bins that have particles)
    mask = cnt > 0
    err_rho = np.mean(np.abs(rho_n[mask] - rho_e[mask]))
    err_p = np.mean(np.abs(p_n[mask] - p_e[mask]))
    err_u = np.mean(np.abs(u_n[mask] - u_e[mask]))
    summary = (f"Lax shock tube comparison at t={t:.6f}\n"
               f"L1(rho) = {err_rho:.6e}\n"
               f"L1(p)   = {err_p:.6e}\n"
               f"L1(u)   = {err_u:.6e}\n")
    print(summary)
    with open(os.path.join(res_dir, "error_summary.txt"), "w") as f:
        f.write(summary)

    fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    titles = ["Density", "Pressure", "Velocity (x)"]
    num = [rho_n, p_n, u_n]
    exa = [rho_e, p_e, u_e]
    for ax, ttl, nn, ee in zip(axes, titles, num, exa):
        ax.plot(xc, ee, "k-", lw=1.5, label="exact")
        ax.plot(xc, nn, "r.", ms=3, label="SPH")
        ax.set_ylabel(ttl)
        ax.legend(loc="best")
        ax.grid(True, alpha=0.3)
    axes[-1].set_xlabel("x")
    fig.suptitle(f"Lax shock tube (t={t:.4f})")
    fig.tight_layout()
    fig_path = os.path.join(res_dir, f"comparison_{step}.png")
    fig.savefig(fig_path, dpi=120)
    print(f"Saved {fig_path}")


if __name__ == "__main__":
    main()
