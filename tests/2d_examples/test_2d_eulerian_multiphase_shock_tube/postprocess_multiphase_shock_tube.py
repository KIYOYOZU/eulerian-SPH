#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Post-process for the 2D Eulerian SPH multiphase shock tube (Kapila
five-equation model, stiffened gas EOS).

Two modes, auto-detected from the data:

  * Sod mode (TEST_CASE==2): two ideal-gas materials, alpha=1 | alpha=0, Sod
    states. The mixture reduces to the ideal-gas Euler equations, so rho/p/u
    are compared against the exact Sod (Toro) solution and alpha against the
    material contact riding on the Sod contact. This reproduces the classic
    rarefaction-fan + contact + shock curves.

  * Gas-water mode (TEST_CASE==1): high-pressure gas | water. Compared against
    the exact two-phase Riemann solution (gas isentropic rarefaction + water
    stiffened-gas shock Hugoniot).

Reads particles_<step>.csv (columns x, y, rho, p, u, alpha), bins along x,
checks alpha in [0,1], overlays the exact solution and saves a figure.

Usage:
    python postprocess_multiphase_shock_tube.py [step]
"""

import os
import sys
import glob
import math

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

DOMAIN_X0, DOMAIN_X1 = 0.0, 1.0

# ---------------- Sod-like two-gas (TEST_CASE==2) parameters -------------
SOD_GAMMA_L, SOD_PINF_L = 1.4, 0.0
SOD_GAMMA_R, SOD_PINF_R = 1.6, 0.0
SOD_RHO_L, SOD_P_L, SOD_U_L = 1.0, 0.425, 0.0
SOD_RHO_R, SOD_P_R, SOD_U_R = 0.125, 0.1, 0.0
SOD_X0 = 0.5

# ---------------- Gas-water (TEST_CASE==1) parameters ---------------------
GW_GAMMA_L, GW_PINF_L = 1.4, 0.0
GW_GAMMA_R, GW_PINF_R = 4.4, 6.0e8
GW_RHO_L, GW_P_L, GW_U_L = 1.4, 1.0e6, 0.0
GW_RHO_R, GW_P_R, GW_U_R = 1000.0, 1.0e5, 0.0
GW_X0 = 0.7


# =====================================================================
# General two-phase exact Riemann solver (stiffened gas on each side).
# Uses shifted pressure P = p + p_inf so each side is ideal-gas-like.
# Handles rarefaction or shock on each side; contact carries the alpha jump.
# =====================================================================
def _c(gamma, pinf, p, rho):
    return math.sqrt(gamma * (p + pinf) / rho)


def _wave_f(p, gamma, pinf, rho_k, p_k):
    """Stiffened-gas wave curve f(p) and derivative (Toro 4.5/4.6, shifted)."""
    P = p + pinf
    Pk = p_k + pinf
    if P <= 0.0:
        return 1.0e18, 0.0
    c_k = math.sqrt(gamma * Pk / rho_k)
    if p <= p_k:  # rarefaction
        f = 2.0 * c_k / (gamma - 1.0) * ((P / Pk) ** ((gamma - 1.0) / (2.0 * gamma)) - 1.0)
        df = (1.0 / (rho_k * c_k)) * (P / Pk) ** (-(gamma + 1.0) / (2.0 * gamma))
    else:  # shock
        a_k = 2.0 / ((gamma + 1.0) * rho_k)
        b_k = (gamma - 1.0) / (gamma + 1.0) * Pk
        sq = math.sqrt(a_k / (b_k + P))
        f = (P - Pk) * sq
        df = sq * (1.0 - 0.5 * (P - Pk) / (b_k + P))
    return f, df


def tp_exact_solution(x, t, cfg):
    gL, piL, rL, pL, uL = cfg["gL"], cfg["piL"], cfg["rL"], cfg["pL"], cfg["uL"]
    gR, piR, rR, pR, uR = cfg["gR"], cfg["piR"], cfg["rR"], cfg["pR"], cfg["uR"]
    x0 = cfg["x0"]
    x = np.asarray(x, dtype=float)
    rho = np.empty_like(x)
    p = np.empty_like(x)
    u = np.empty_like(x)
    alpha = np.empty_like(x)

    # Newton iteration for star pressure (p continuous at contact).
    p_star = 0.5 * (pL + pR)
    for _ in range(80):
        fL, dfL = _wave_f(p_star, gL, piL, rL, pL)
        fR, dfR = _wave_f(p_star, gR, piR, rR, pR)
        f = fL + fR + (uR - uL)
        df = dfL + dfR
        dp = -f / df if df != 0.0 else -f
        p_new = max(p_star + dp, 1.0e-9)
        if abs(p_new - p_star) / (0.5 * (p_new + p_star)) < 1e-11:
            p_star = p_new
            break
        p_star = p_new
    fL, _ = _wave_f(p_star, gL, piL, rL, pL)
    fR, _ = _wave_f(p_star, gR, piR, rR, pR)
    u_star = 0.5 * (uL + uR + fR - fL)

    cL = _c(gL, piL, pL, rL)
    cR = _c(gR, piR, pR, rR)
    PL, PR = pL + piL, pR + piR
    PstarL, PstarR = p_star + piL, p_star + piR

    def rho_star_shock(g, rk, Pk, Pstar):
        return rk * ((Pstar / Pk + (g - 1) / (g + 1)) / ((g - 1) / (g + 1) * Pstar / Pk + 1.0))

    # Left wave (material L)
    if p_star <= pL:  # left rarefaction
        left_raref = True
        rho_sL = rL * (PstarL / PL) ** (1.0 / gL)
        c_sL = cL * (PstarL / PL) ** ((gL - 1.0) / (2.0 * gL))
        S_HL = x0 + (uL - cL) * t
        S_TL = x0 + (u_star - c_sL) * t
    else:  # left shock
        left_raref = False
        rho_sL = rho_star_shock(gL, rL, PL, PstarL)
        S_L = x0 + (uL - cL * math.sqrt((gL + 1) / (2 * gL) * (PstarL / PL - 1) + 1)) * t
    # Right wave (material R)
    if p_star <= pR:  # right rarefaction
        right_raref = True
        rho_sR = rR * (PstarR / PR) ** (1.0 / gR)
        c_sR = cR * (PstarR / PR) ** ((gR - 1.0) / (2.0 * gR))
        S_HR = x0 + (uR + cR) * t
        S_TR = x0 + (u_star + c_sR) * t
    else:  # right shock
        right_raref = False
        rho_sR = rho_star_shock(gR, rR, PR, PstarR)
        S_R = x0 + (uR + cR * math.sqrt((gR + 1) / (2 * gR) * (PstarR / PR - 1) + 1)) * t

    x_contact = x0 + u_star * t
    for i, xi in enumerate(x):
        s = (xi - x0) / t
        if s < u_star:  # left of contact (material L)
            if left_raref:
                if xi < S_HL:
                    rho[i], p[i], u[i] = rL, pL, uL
                elif xi > S_TL:
                    rho[i], p[i], u[i] = rho_sL, p_star, u_star
                else:
                    f = 2.0 / (gL + 1) + (gL - 1) / ((gL + 1) * cL) * (uL - s)
                    P_f = PL * f ** (2 * gL / (gL - 1))
                    p[i] = P_f - piL
                    rho[i] = rL * f ** (2 / (gL - 1))
                    u[i] = 2.0 / (gL + 1) * (cL + (gL - 1) / 2 * uL + s)
            else:
                if xi < S_L:
                    rho[i], p[i], u[i] = rL, pL, uL
                else:
                    rho[i], p[i], u[i] = rho_sL, p_star, u_star
        else:  # right of contact (material R)
            if right_raref:
                if xi > S_HR:
                    rho[i], p[i], u[i] = rR, pR, uR
                elif xi < S_TR:
                    rho[i], p[i], u[i] = rho_sR, p_star, u_star
                else:
                    f = 2.0 / (gR + 1) - (gR - 1) / ((gR + 1) * cR) * (uR - s)
                    P_f = PR * f ** (2 * gR / (gR - 1))
                    p[i] = P_f - piR
                    rho[i] = rR * f ** (2 / (gR - 1))
                    u[i] = 2.0 / (gR + 1) * (-cR + (gR - 1) / 2 * uR + s)
            else:
                if xi > S_R:
                    rho[i], p[i], u[i] = rR, pR, uR
                else:
                    rho[i], p[i], u[i] = rho_sR, p_star, u_star
        alpha[i] = 1.0 if xi < x_contact else 0.0
    return rho, p, u, alpha, dict(p_star=p_star, u_star=u_star, x_contact=x_contact)


def make_cfg(sod_mode):
    if sod_mode:
        return dict(gL=SOD_GAMMA_L, piL=SOD_PINF_L, rL=SOD_RHO_L, pL=SOD_P_L, uL=SOD_U_L,
                    gR=SOD_GAMMA_R, piR=SOD_PINF_R, rR=SOD_RHO_R, pR=SOD_P_R, uR=SOD_U_R, x0=SOD_X0)
    return dict(gL=GW_GAMMA_L, piL=GW_PINF_L, rL=GW_RHO_L, pL=GW_P_L, uL=GW_U_L,
                gR=GW_GAMMA_R, piR=GW_PINF_R, rR=GW_RHO_R, pR=GW_P_R, uR=GW_U_R, x0=GW_X0)


# =====================================================================
def load_particles(csv_path):
    return np.loadtxt(csv_path, delimiter=",", skiprows=1)


def bin_profile(data, n_bins):
    x_edges = np.linspace(DOMAIN_X0, DOMAIN_X1, n_bins + 1)
    centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    rho = np.zeros(n_bins)
    p = np.zeros(n_bins)
    u = np.zeros(n_bins)
    alpha = np.zeros(n_bins)
    cnt = np.zeros(n_bins, dtype=int)
    for row in data:
        idx = int((row[0] - DOMAIN_X0) / (DOMAIN_X1 - DOMAIN_X0) * n_bins)
        idx = min(max(idx, 0), n_bins - 1)
        rho[idx] += row[2]
        p[idx] += row[3]
        u[idx] += row[4]
        alpha[idx] += row[5]
        cnt[idx] += 1
    mask = cnt > 0
    rho[mask] /= cnt[mask]
    p[mask] /= cnt[mask]
    u[mask] /= cnt[mask]
    alpha[mask] /= cnt[mask]
    return centers, rho, p, u, alpha, cnt


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

    time_path = os.path.join(out_dir, f"time_{step}.txt")
    t = float(open(time_path).read().strip()) if os.path.exists(time_path) else 0.0

    data = load_particles(csv_path)
    sod_mode = data[:, 2].max() < 10.0  # Sod densities <= ~1; gas-water reaches 1000
    n_bins = 400 if sod_mode else 200
    centers, rho_n, p_n, u_n, alpha_n, cnt = bin_profile(data, n_bins)

    alpha_raw = data[:, 5]
    a_min, a_max = float(alpha_raw.min()), float(alpha_raw.max())
    bounded = (a_min >= -1e-8) and (a_max <= 1 + 1e-8)

    tt = max(t, 1e-12)
    cfg = make_cfg(sod_mode)
    rho_e, p_e, u_e, alpha_e, meta = tp_exact_solution(centers, tt, cfg)
    mode = "two-gas Sod (g=1.4|1.6)" if sod_mode else "gas-water"

    mask = cnt > 0
    dx = (DOMAIN_X1 - DOMAIN_X0) / len(centers)
    l1_rho = float(np.sum(np.abs(rho_n[mask] - rho_e[mask])) * dx)
    l1_p = float(np.sum(np.abs(p_n[mask] - p_e[mask])) * dx)
    l1_u = float(np.sum(np.abs(u_n[mask] - u_e[mask])) * dx)

    summary = (
        f"[{mode}] multiphase shock tube at t={t:.6e} (step={step})\n"
        f"  particles  : {data.shape[0]}\n"
        f"  alpha range: [{a_min:.6f}, {a_max:.6f}]  bounded={bounded}\n"
        f"  exact      : p*={meta['p_star']:.6e}  u*={meta['u_star']:.4f}  contact x={meta['x_contact']:.5f}\n"
        f"  L1(rho)={l1_rho:.4e}  L1(p)={l1_p:.4e}  L1(u)={l1_u:.4e}\n"
    )
    print(summary)
    with open(os.path.join(res_dir, "profile_summary.txt"), "w") as f:
        f.write(summary)

    x_ref = np.linspace(DOMAIN_X0, DOMAIN_X1, 2000)
    rE, pE, uE, aE, _ = tp_exact_solution(x_ref, tt, cfg)

    fig, axes = plt.subplots(4, 1, figsize=(8, 11), sharex=True)
    titles = ["Density", "Pressure", "Velocity (x)", r"Volume fraction $\alpha_1$"]
    sph = [rho_n, p_n, u_n, alpha_n]
    ext = [rE, pE, uE, aE]
    for ax, ttl, ss, ee in zip(axes, titles, sph, ext):
        ax.plot(x_ref, ee, "k-", lw=1.5, label="exact")
        ax.plot(centers[mask], ss[mask], "r.", ms=3, label="SPH")
        ax.axvline(GW_X0 if not sod_mode else SOD_X0, color="b", ls="--", lw=0.8, alpha=0.5)
        ax.set_ylabel(ttl)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=8)
    axes[-1].set_xlabel("x")
    fig.suptitle(f"Multiphase shock tube [{mode}] (t={t:.4f})")
    fig.tight_layout()
    fig_path = os.path.join(res_dir, f"profiles_{step}.png")
    fig.savefig(fig_path, dpi=120)
    print(f"Saved {fig_path}")
    return 0 if bounded else 2


if __name__ == "__main__":
    sys.exit(main())
