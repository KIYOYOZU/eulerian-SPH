/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file    eulerian_multiphase_riemann_solver.h
 * @brief   Interface-state computation for the Kapila five-equation two-phase
 *          model: the first-order HLLC solver and the second-order MUSCL-HLLC
 *          bridge (slope-limited reconstruction + mixture-EOS-consistent
 *          energy + the same HLLC solve). State includes volume fraction
 *          alpha1 in addition to the standard compressible variables.
 * @author  KIYOYOZU
 */

#ifndef EULERIAN_MULTIPHASE_RIEMANN_SOLVER_H
#define EULERIAN_MULTIPHASE_RIEMANN_SOLVER_H

#include "multiphase_mixture.h"
#include "muscl_reconstruction.hpp"
#include "riemann_solver.h"

namespace SPH
{
/**
 * @struct MultiphaseFluidState
 * @brief Input state for five-equation Riemann solver.
 *        Extends FluidStateIn with total energy per unit volume
 *        and volume fraction of phase 1.
 */
struct MultiphaseFluidState : FluidStateIn
{
    Real &E_;      /**< total energy per unit volume */
    Real &alpha_;  /**< volume fraction of phase 1 */
    MultiphaseFluidState(Real &rho, Vecd &vel, Real &p, Real &E, Real &alpha)
        : FluidStateIn(rho, vel, p), E_(E), alpha_(alpha) {};
};

/**
 * @struct MultiphaseFluidStarState
 * @brief Output star state from five-equation Riemann solver.
 */
struct MultiphaseFluidStarState : FluidStateOut
{
    Real E_;      /**< total energy per unit volume */
    Real alpha_;  /**< volume fraction of phase 1 */
    MultiphaseFluidStarState(Real rho, Vecd vel, Real p, Real E, Real alpha)
        : FluidStateOut(rho, vel, p), E_(E), alpha_(alpha) {};
};

/**
 * @class MultiphaseHLLCRiemannSolver
 * @brief HLLC Riemann solver for the Kapila five-equation model.
 *        Uses Wood mixture sound speed for wave speed estimation.
 *        Volume fraction is transported upwind.
 */
class MultiphaseHLLCRiemannSolver
{
    // Stored by value: the solver is built from a mixture passed through
    // DynamicsArgs, which copies it into a temporary tuple. Holding a
    // reference to that tuple element would dangle once construction ends.
    // MultiphaseMixture is a small stateless object (two phase references),
    // so copying it is safe and keeps the phase references valid.
    MultiphaseMixture mixture_;

  public:
    MultiphaseHLLCRiemannSolver(MultiphaseMixture &mixture);

    /**
     * @brief Compute the interface star state.
     * @param state_i Left state (particle i).
     * @param state_j Right state (particle j).
     * @param e_ij Unit vector from j to i (SPH convention: e_ij = (x_i - x_j)/|x_i - x_j|).
     * @return Star state at the interface.
     */
    MultiphaseFluidStarState getInterfaceState(
        const MultiphaseFluidState &state_i,
        const MultiphaseFluidState &state_j,
        const Vecd &e_ij);
};

/**
 * @class MultiphaseMUSCLBridge
 * @brief Second-order interface states: reconstruct left/right primitives
 *        (rho, vel, p, alpha) with a slope limiter and rebuild the mixture
 *        total energy EOS-consistently, then call the five-equation HLLC
 *        solver. Pairs across a jump of the stiffened invariant B(alpha)
 *        fall back to piecewise-constant (first-order) states, which is
 *        required for stability at stiff material interfaces (gas-water).
 *        With SecondOrderConfig::piecewise_rho_alpha, rho/alpha stay
 *        piecewise constant everywhere (vel+p reconstruction only).
 */
class MultiphaseMUSCLBridge
{
  public:
    MultiphaseMUSCLBridge(MultiphaseMixture &mixture, const fluid_dynamics::SecondOrderConfig &cfg);

#if SPH_NDIM == 2
    MultiphaseFluidStarState getInterfaceState(
        const MultiphaseFluidState &state_i, const MultiphaseFluidState &state_j,
        const Vecd &xi, const Vecd &xj, const Vecd &xf, const Vecd &e_ij,
        const Vecd &grad_rho_i, const Vecd &grad_rho_j,
        const Vecd &grad_u_i, const Vecd &grad_u_j,
        const Vecd &grad_v_i, const Vecd &grad_v_j,
        const Vecd &grad_p_i, const Vecd &grad_p_j,
        const Vecd &grad_a_i, const Vecd &grad_a_j);
#else
    MultiphaseFluidStarState getInterfaceState(
        const MultiphaseFluidState &state_i, const MultiphaseFluidState &state_j,
        const Vecd &xi, const Vecd &xj, const Vecd &xf, const Vecd &e_ij,
        const Vecd &grad_rho_i, const Vecd &grad_rho_j,
        const Vecd &grad_u_i, const Vecd &grad_u_j,
        const Vecd &grad_v_i, const Vecd &grad_v_j,
        const Vecd &grad_w_i, const Vecd &grad_w_j,
        const Vecd &grad_p_i, const Vecd &grad_p_j,
        const Vecd &grad_a_i, const Vecd &grad_a_j);
#endif

  protected:
    MultiphaseMixture mixture_;
    fluid_dynamics::SecondOrderConfig cfg_;
    MultiphaseHLLCRiemannSolver riemann_solver_;
};

} // namespace SPH

#endif // EULERIAN_MULTIPHASE_RIEMANN_SOLVER_H
