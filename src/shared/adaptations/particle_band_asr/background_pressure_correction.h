/**
 * @file 	background_pressure_correction.h
 * @brief 	Well-balanced consistent-flux correction for the Eulerian
 * 			multiphase Godunov scheme on graded band lattices.
 * @author 	KIYOYOZU
 */

#ifndef BACKGROUND_PRESSURE_CORRECTION_H
#define BACKGROUND_PRESSURE_CORRECTION_H

#include "base_body_relation.h"
#include "base_particles.h"
#include "multiphase_mixture.h"

namespace SPH
{
/**
 * @class MultiphaseBackgroundPressureCorrection
 * @brief Cancels the spurious force that ANY locally uniform flow state exerts
 * on particles whose kernel stencil is asymmetric -- band boundaries of a
 * graded lattice, where the kernel first moment M_i = sum_j V_j dW_ij e_ij is
 * nonzero. For a uniform state the exact pair fluxes reduce to the frozen-state
 * fluxes, and the discrete divergence collapses to F . M_i for mass
 * (F = rho u), momentum (F = rho u(x)u + p I) and energy (F = (E + p) u).
 * Each channel's frozen contribution is subtracted per pair, using the
 * symmetric average of the two endpoint states (with an optional advective
 * upwind endpoint selection at strong density jumps, see
 * execMomentumBetweenHalves).
 *
 * Two-stage timing (matches where each channel's spurious flux is injected):
 * - execMomentumBetweenHalves: after the 1st half (momentum) and BEFORE the
 *   2nd half. The 2nd half's mass/energy/alpha fluxes are built from the
 *   current velocity; leaving the band-boundary pseudo-velocity in place
 *   feeds it into the alpha upwind term at material interfaces (where
 *   alpha* != alpha_i) and into the EOS pressure recovery, which then
 *   corrupts the next step's momentum flux -- a destabilizing feedback loop.
 * - execMassEnergyAfterHalves: after the 2nd half, which is where the frozen
 *   mass/energy flux divergences are actually injected. Recovers rho, vel and
 *   the EOS pressure from the corrected conserved state.
 * The volume fraction needs no correction: its update is the upwind form
 * u*_n (alpha_i - alpha*), which vanishes identically for any uniform state.
 *
 * Properties:
 * - exact well-balancedness: any uniform state stays untouched (the Riemann
 *   interface state of identical inputs is the input state itself);
 * - global conservation: inner-pair corrections are pairwise antisymmetric
 *   (symmetric cached dW, e_ji = -e_ij, symmetric endpoint average);
 * - identity on uniform lattices: M_i = 0 makes every pair sum vanish;
 * - it subsumes the constant background-pressure gauge transformation
 *   p I -> (p - p0) I of the momentum flux proven by the stationary gates.
 *
 * Wall pairs enter the momentum stage only, with the frozen wall flux p_i e
 * (the reflective-ghost identity p* = p_i, exact at rest). The exact star
 * pressure under wall-normal flow must NOT be used here: it carries a
 * rho c u_n term whose sign makes the correction a positive feedback loop on
 * any residual wall-normal velocity. Wall mass/energy fluxes vanish
 * identically because the mirror star state has u* . e = 0 exactly.
 *
 * energy_correlation = false reduces the correction to the momentum stage
 * (the pressure gauge proven by the stationary tests); the mass/energy
 * stage is skipped entirely so that no thermodynamically inconsistent partial
 * recovery can occur.
 */
class MultiphaseBackgroundPressureCorrection
{
  public:
	MultiphaseBackgroundPressureCorrection(BaseInnerRelation &inner_relation,
										   BaseContactRelation &wall_contact_relation,
										   MultiphaseMixture &mixture,
										   bool energy_correlation = true);

	/** capture the state the two half steps will consume as flux input.
	 *  Call once per step, before the 1st half. The 1st half overwrites
	 *  vel/mom with the raw (uncorrected) band-boundary kick before the
	 *  momentum stage runs; building the frozen fluxes from that polluted
	 *  state injects a systematic O(rho |kick|^2) error -- 4e-4 m/s/step on
	 *  the gas side of a quiescent run, saturating at ~1e-2 m/s. The snapshot restores
	 *  bit-level cancellation (star vs endpoint-average residual ~1e-15). */
	void snapshotState();
	/** apply between the two half steps: momentum channel + velocity recovery.
	 *  upwind_select: replace the symmetric endpoint average by the advective
	 *  upwind endpoint (chosen by the sign of the normal relative velocity) for
	 *  pairs whose density jumps by more than interface_density_ratio. At a
	 *  pressure-equilibrated contact with continuous velocity the exact star
	 *  flux IS the upwind endpoint flux, while the arithmetic average misses it
	 *  by (rho* - rho_bar) u(x)u -- a spurious ~10 m/s/step kick on the light
	 *  side of a gas-water interface at u=100 (measured on a translating gas-water interface). The selection is
	 *  swap-invariant (the same physical upwind particle from either side) and
	 *  the advective flux is linear in the state, so exact pairwise
	 *  antisymmetry -- and thus global conservation -- is preserved. */
	void execMomentumBetweenHalves(Real dt, bool upwind_select = false,
								   Real interface_density_ratio = 2.0);
	/** apply after the 2nd-half update: mass/energy channels + full recovery.
	 *  upwind_select as in execMomentumBetweenHalves: at an equilibrated
	 *  contact the star mass flux rho*u and energy flux (E*+p*)u* equal the
	 *  upwind endpoint's exactly, while the arithmetic average misses them by
	 *  +-0.5 (rho_j - rho_i) u -- a per-particle mass/energy injection at
	 *  interface band-boundary particles (global conservation is kept either
	 *  way by pairwise antisymmetry). */
	void execMassEnergyAfterHalves(Real dt, bool upwind_select = false,
								   Real interface_density_ratio = 2.0);

  protected:
	BaseParticles &particles_;
	ParticleConfiguration &inner_configuration_;
	BaseContactRelation &wall_contact_;
	MultiphaseMixture &mixture_;
	bool energy_correlation_;
	Real *Vol_, *mass_, *rho_, *p_, *E_, *alpha_;
	Vecd *vel_, *mom_;
	// step-start state consumed by the frozen fluxes (see snapshotState)
	StdVec<Vecd> vel_snap_;
	StdVec<Real> rho_snap_, p_snap_, E_snap_;
};

} // namespace SPH

#endif // BACKGROUND_PRESSURE_CORRECTION_H
