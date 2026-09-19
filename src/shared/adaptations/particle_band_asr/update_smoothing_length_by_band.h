/**
 * @file 	update_smoothing_length_by_band.h
 * @brief 	Smoothing length evolution of the particle-band ASR method,
 * 			Eqs. (6)-(9) of Yang, Kong & Liu, Phys. Rev. E 104, 055308 (2021).
 * @author 	KIYOYOZU
 */

#ifndef UPDATE_SMOOTHING_LENGTH_BY_BAND_H
#define UPDATE_SMOOTHING_LENGTH_BY_BAND_H

#include "particle_band_adaptation.h"

#include "base_body_relation.h"
#include "base_particle_dynamics.h"

namespace SPH
{
/**
 * @class UpdateSmoothingLengthByBand
 * @brief Evolves the local smoothing length towards the local particle number.
 * @details h~ = h/2 [1 + (N_r/N_i)^{1/2}] (Eq. 6), h_bar = neighbor average of
 * 			the previous smoothing lengths (Eq. 7), h_new = (h~ + h_bar)/2
 * 			(Eq. 8), all clamped to [0.5 h_r, min(1.5 h_r, h_ref)] with
 * 			h_r = h_spacing_ratio * V^{1/d} (Eq. 9, upper bound tightened from
 * 			the paper's 2 h_r so that h_ratio >= 1 always holds, which the
 * 			multi-level cell linked list requires). The update is two-pass:
 * 			all neighbor statistics use the previous smoothing lengths.
 */
class UpdateSmoothingLengthByBand
{
  public:
    UpdateSmoothingLengthByBand(RealBody &real_body, BaseInnerRelation &inner_relation);
    void exec();

  protected:
    ParticleBandAdaptation &adaptation_;
    BaseParticles &particles_;
    ParticleConfiguration &inner_configuration_;
    Real h_ref_, n_r_, h_spacing_ratio_;
    Real *h_ratio_, *Vol_;
    StdVec<Real> h_new_; /**< temporary buffer for the two-pass update */
};
} // namespace SPH
#endif // UPDATE_SMOOTHING_LENGTH_BY_BAND_H
