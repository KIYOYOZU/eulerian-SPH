#include "eulerian_multiphase_riemann_solver.h"

namespace SPH
{
//=================================================================================================//
MultiphaseHLLCRiemannSolver::MultiphaseHLLCRiemannSolver(MultiphaseMixture &mixture)
    : mixture_(mixture) {};
//=================================================================================================//
MultiphaseFluidStarState MultiphaseHLLCRiemannSolver::
    getInterfaceState(const MultiphaseFluidState &state_i,
                      const MultiphaseFluidState &state_j,
                      const Vecd &e_ij)
{
    // Normal velocities projected onto -e_ij (from i towards j).
    Real ul = -e_ij.dot(state_i.vel_);
    Real ur = -e_ij.dot(state_j.vel_);

    // Wood mixture sound speeds for wave speed estimation.
    Real c_l = mixture_.MixtureSoundSpeed(state_i.alpha_, state_i.rho_, state_i.p_);
    Real c_r = mixture_.MixtureSoundSpeed(state_j.alpha_, state_j.rho_, state_j.p_);

    // Wave speed estimates (Davis bounds).
    Real s_l = SMIN(ul - c_l, ur - c_r);
    Real s_r = SMAX(ul + c_l, ur + c_r);

    // Contact wave speed from pressure continuity.
    // s* = (p_r - p_l + rho_l*ul*(s_l - ul) - rho_r*ur*(s_r - ur))
    //      / (rho_l*(s_l - ul) - rho_r*(s_r - ur))
    Real denominator = state_i.rho_ * (s_l - ul) - state_j.rho_ * (s_r - ur);
    Real s_star;
    if (std::abs(denominator) < 1e-14)
    {
        s_star = 0.5 * (ul + ur);
    }
    else
    {
        s_star = (state_j.p_ - state_i.p_
                  + state_i.rho_ * ul * (s_l - ul)
                  - state_j.rho_ * ur * (s_r - ur)) / denominator;
        // Robustness guard (as in the single-phase reference): a contact speed far
        // outside the wave-speed range signals a degenerate solve; fall back to the
        // velocity average rather than poisoning the fluxes.
        if (s_star < s_l - 1000.0 || s_star > s_r + 1000.0)
        {
            s_star = 0.5 * (ul + ur);
        }
    }

    // Star-region pressure, estimated from each side (the two coincide at the
    // contact; averaging is numerically robust). Only the two star branches use it.
    Real p_star_l = state_i.p_ + state_i.rho_ * (s_l - ul) * (s_star - ul);
    Real p_star_r = state_j.p_ + state_j.rho_ * (s_r - ur) * (s_star - ur);

    Real p_star = 0.0;
    Real rho_star = 0.0;
    Vecd v_star = Vecd::Zero();
    Real energy_star = 0.0;
    Real alpha_star = 0.0;

    if (0.0 < s_l)
    {
        // Supersonic from the left: the interface state is the left state itself,
        // so its pressure is p_l (not the star pressure).
        p_star = state_i.p_;
        rho_star = state_i.rho_;
        v_star = state_i.vel_;
        energy_star = state_i.E_;
        alpha_star = state_i.alpha_;
    }
    else if (s_l <= 0.0 && 0.0 <= s_star)
    {
        // Left star region. The star velocity has normal component s_star along
        // the left-to-right normal n = -e_ij (e_ij points from j to i):
        // v_star = vel_i + n * (s_star - ul) = vel_i - e_ij * (s_star - ul).
        p_star = 0.5 * (p_star_l + p_star_r);
        Real factor = state_i.rho_ * (s_l - ul) / (s_l - s_star);
        rho_star = factor;
        v_star = state_i.vel_ - e_ij * (s_star - ul);
        energy_star = ((s_l - ul) * state_i.E_ - state_i.p_ * ul + p_star * s_star) / (s_l - s_star);
        alpha_star = state_i.alpha_;
    }
    else if (s_star <= 0.0 && 0.0 <= s_r)
    {
        // Right star region (same normal-direction convention as above).
        p_star = 0.5 * (p_star_l + p_star_r);
        Real factor = state_j.rho_ * (s_r - ur) / (s_r - s_star);
        rho_star = factor;
        v_star = state_j.vel_ - e_ij * (s_star - ur);
        energy_star = ((s_r - ur) * state_j.E_ - state_j.p_ * ur + p_star * s_star) / (s_r - s_star);
        alpha_star = state_j.alpha_;
    }
    else
    {
        // Supersonic from the right: the interface state is the right state itself,
        // so its pressure is p_r (not the star pressure).
        p_star = state_j.p_;
        rho_star = state_j.rho_;
        v_star = state_j.vel_;
        energy_star = state_j.E_;
        alpha_star = state_j.alpha_;
    }

    // Clamp volume fraction to [0, 1].
    alpha_star = SMAX(0.0, SMIN(1.0, alpha_star));

    return MultiphaseFluidStarState(rho_star, v_star, p_star, energy_star, alpha_star);
}
//=================================================================================================//
MultiphaseMUSCLBridge::MultiphaseMUSCLBridge(MultiphaseMixture &mixture,
                                                 const fluid_dynamics::SecondOrderConfig &cfg)
    : mixture_(mixture), cfg_(cfg), riemann_solver_(mixture) {}
//=================================================================================================//
#if SPH_NDIM == 2
MultiphaseFluidStarState MultiphaseMUSCLBridge::getInterfaceState(
    const MultiphaseFluidState &state_i, const MultiphaseFluidState &state_j,
    const Vecd &xi, const Vecd &xj, const Vecd &xf, const Vecd &e_ij,
    const Vecd &grad_rho_i, const Vecd &grad_rho_j,
    const Vecd &grad_u_i, const Vecd &grad_u_j,
    const Vecd &grad_v_i, const Vecd &grad_v_j,
    const Vecd &grad_p_i, const Vecd &grad_p_j,
    const Vecd &grad_a_i, const Vecd &grad_a_j)
#else
MultiphaseFluidStarState MultiphaseMUSCLBridge::getInterfaceState(
    const MultiphaseFluidState &state_i, const MultiphaseFluidState &state_j,
    const Vecd &xi, const Vecd &xj, const Vecd &xf, const Vecd &e_ij,
    const Vecd &grad_rho_i, const Vecd &grad_rho_j,
    const Vecd &grad_u_i, const Vecd &grad_u_j,
    const Vecd &grad_v_i, const Vecd &grad_v_j,
    const Vecd &grad_w_i, const Vecd &grad_w_j,
    const Vecd &grad_p_i, const Vecd &grad_p_j,
    const Vecd &grad_a_i, const Vecd &grad_a_j)
#endif
{
    using fluid_dynamics::LR;
    using fluid_dynamics::Primitives;
    using fluid_dynamics::reconstruct_primitives_muscl;
    using fluid_dynamics::reconstruct_scalar_muscl;

    Primitives Pi{state_i.rho_, state_i.vel_, state_i.p_, state_i.E_};
    Primitives Pj{state_j.rho_, state_j.vel_, state_j.p_, state_j.E_};

    // Slope suppression masks.
    // (1) piecewise_rho_alpha (config): keep rho/alpha piecewise constant.
    // (2) across_interface: pairs straddling a jump of the stiffened-gas
    //     invariant B(alpha) (e.g. gas-water, where B jumps by O(gamma*p_inf))
    //     fall back to piecewise-constant (first-order) interface states for
    //     ALL fields; reconstructing across that energy-scale jump drives
    //     velocity/pressure spikes at the contact. Contacts with equal B on
    //     both sides (e.g. two-gas Sod, B=0) keep full MUSCL accuracy.
    const Vecd zero_slope = Vecd::Zero();
    const Real B_i = mixture_.MixtureB(state_i.alpha_);
    const Real B_j = mixture_.MixtureB(state_j.alpha_);
    const bool across_interface =
        std::abs(B_i - B_j) > (Real)1e-9 * SMAX(std::abs(B_i), std::abs(B_j), (Real)1);
    const bool rho_alpha_pw = cfg_.piecewise_rho_alpha || across_interface;
    const Vecd &gr_i = rho_alpha_pw ? zero_slope : grad_rho_i;
    const Vecd &gr_j = rho_alpha_pw ? zero_slope : grad_rho_j;
    const Vecd &ga_i = rho_alpha_pw ? zero_slope : grad_a_i;
    const Vecd &ga_j = rho_alpha_pw ? zero_slope : grad_a_j;
    const Vecd &gu_i = across_interface ? zero_slope : grad_u_i;
    const Vecd &gu_j = across_interface ? zero_slope : grad_u_j;
    const Vecd &gv_i = across_interface ? zero_slope : grad_v_i;
    const Vecd &gv_j = across_interface ? zero_slope : grad_v_j;
#if SPH_NDIM == 3
    const Vecd &gw_i = across_interface ? zero_slope : grad_w_i;
    const Vecd &gw_j = across_interface ? zero_slope : grad_w_j;
#endif
    const Vecd &gp_i = across_interface ? zero_slope : grad_p_i;
    const Vecd &gp_j = across_interface ? zero_slope : grad_p_j;

    LR lr = reconstruct_primitives_muscl(
        Pi, Pj,
        gr_i, gr_j,
        gu_i, gu_j,
        gv_i, gv_j,
#if SPH_NDIM == 3
        gw_i, gw_j,
#endif
        gp_i, gp_j,
        xi, xj, xf, cfg_);

    // Volume fraction reconstructed and clamped to [0, 1].
    auto a_lr = reconstruct_scalar_muscl(state_i.alpha_, ga_i,
                                         state_j.alpha_, ga_j,
                                         xi, xj, xf, cfg_);
    Real alpha_L = SMAX(0.0, SMIN(1.0, a_lr.first));
    Real alpha_R = SMAX(0.0, SMIN(1.0, a_lr.second));

    Real rho_L = lr.L.rho, rho_R = lr.R.rho;
    Real p_L = lr.L.p, p_R = lr.R.p;
    Vecd v_L = lr.L.vel, v_R = lr.R.vel;

    // EOS-consistent mixture total energy (per volume) from reconstructed state.
    Real E_L = mixture_.MixtureInternalEnergyPerVolume(alpha_L, p_L)
             + 0.5 * rho_L * v_L.squaredNorm();
    Real E_R = mixture_.MixtureInternalEnergyPerVolume(alpha_R, p_R)
             + 0.5 * rho_R * v_R.squaredNorm();

    Real rhoL_l = rho_L, pL_l = p_L, EL_l = E_L, aL_l = alpha_L;
    Real rhoR_l = rho_R, pR_l = p_R, ER_l = E_R, aR_l = alpha_R;
    Vecd vL_l = v_L, vR_l = v_R;

    MultiphaseFluidState Ls(rhoL_l, vL_l, pL_l, EL_l, aL_l);
    MultiphaseFluidState Rs(rhoR_l, vR_l, pR_l, ER_l, aR_l);

    return riemann_solver_.getInterfaceState(Ls, Rs, e_ij);
}
//=================================================================================================//
} // namespace SPH
