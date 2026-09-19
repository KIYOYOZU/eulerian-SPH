/**
 * @file 	shepard_density_filter.h
 * @brief 	Shepard density filter of the particle-band ASR method,
 * 			Eq. (25) of Yang, Kong & Liu, Phys. Rev. E 104, 055308 (2021).
 * @author 	KIYOYOZU
 */

#ifndef SHEPARD_DENSITY_FILTER_H
#define SHEPARD_DENSITY_FILTER_H

#include "multiphase_mixture.h"

#include "base_body_relation.h"

namespace SPH
{
/**
 * @class ShepardDensityFilter
 * @brief rho_hat_i = sum_j m_j W_ij / sum_j V_j W_ij over same-phase
 * neighbors only (alpha on the same side of 0.5), followed by the pressure
 * recovery p_hat = EOS(alpha_i, E_i/V_i - rho_hat u_i^2 / 2).
 *
 * Override-style insertion: run after the 2nd half step (which recomputes
 * rho = mass/Vol and p from the energy equation), so the next step's fluxes
 * see the filtered state. The conserved bookkeeping variables mass, momentum
 * and total energy are never touched, hence global conservation is exact by
 * construction. Raw values are backed up to "DensityRaw"/"PressureRaw".
 */
class ShepardDensityFilter
{
  public:
	ShepardDensityFilter(RealBody &real_body, BaseInnerRelation &inner_relation,
						 MultiphaseMixture &mixture);
	void exec();

  protected:
	BaseParticles &particles_;
	ParticleConfiguration &inner_configuration_;
	MultiphaseMixture &mixture_;
	Real *Vol_, *mass_, *rho_, *p_, *E_, *alpha_, *rho_raw_, *p_raw_;
	Vecd *vel_;
};

} // namespace SPH

#endif // SHEPARD_DENSITY_FILTER_H
