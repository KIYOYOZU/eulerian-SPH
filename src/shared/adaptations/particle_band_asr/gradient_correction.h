/**
 * @file 	gradient_correction.h
 * @brief 	Zeroth-order consistent kernel-gradient correction for the
 * 			particle-band ASR Eulerian Godunov scheme.
 * @author 	KIYOYOZU
 */

#ifndef GRADIENT_CORRECTION_H
#define GRADIENT_CORRECTION_H

#include "base_body_relation.h"
#include "base_particles.h"

namespace SPH
{
/**
 * @class ComputeGradientCorrection
 * @brief Per-particle gradient-bias vector c_i that closes the kernel first
 * 		moment on a graded band lattice.
 * @details The Eulerian Godunov flux of a locally uniform state collapses to
 * 		-2 V_i F . M_i with the kernel first moment
 * 		M_i = sum_j V_j dW_ij e_ij. M_i vanishes to machine precision on a
 * 		uniform lattice but is O(1/h) at band boundaries and at the wall (the
 * 		mirror wall stencil is cut off asymmetrically at the support radius),
 * 		which is the spurious-force source behind the refined-band
 * 		pseudo-velocity. Unlike MultiphaseBackgroundPressureCorrection -- which
 * 		subtracts a per-pair frozen flux and thereby removes the physical
 * 		divergence of any smooth non-uniform flow -- this correction only
 * 		removes the constant-offset part F_i . M_i and leaves the gradient of a
 * 		non-uniform state intact, so it is safe for dynamic (shock-tube) flows.
 *
 * 		The flux loops replace dW_ij e_ij V_j by (dW_ij e_ij - c_i) V_j for the
 * 		mass/momentum/energy channels; then
 * 			sum_j V_j (dW_ij e_ij - c_i) = M_i - c_i VolSum_i = 0
 * 		with c_i = M_i / VolSum_i and VolSum_i = sum_j V_j, so a uniform state
 * 		exerts exactly zero net flux. The alpha upwind term is a difference form
 * 		(alpha_i - alpha*) and is already immune to M_i, so it is left untouched.
 *
 * 		c_i aggregates the inner AND wall-contact stencils so that the combined
 * 		 stencil closes at the wall as well. The variable is registered
 * 		zero-initialized; on a uniform lattice M_i ~ 0 makes c_i ~ 0 and the
 * 		correction is an exact no-op, so non-ASR cases (which never run this
 * 		dynamics) see bit-identical behavior.
 */
class ComputeGradientCorrection
{
  public:
	/** wall_contact_relation may be nullptr for a wall-less domain. */
	ComputeGradientCorrection(BaseInnerRelation &inner_relation,
							  BaseContactRelation *wall_contact_relation);
	void exec();

  protected:
	BaseParticles &particles_;
	ParticleConfiguration &inner_configuration_;
	BaseContactRelation *wall_contact_;
	Real *Vol_;
	Vecd *grad_corr_;
};
} // namespace SPH
#endif // GRADIENT_CORRECTION_H
