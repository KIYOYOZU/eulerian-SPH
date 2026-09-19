/**
 * @file 	update_smoothing_length_by_band.cpp
 * @author 	KIYOYOZU
 */

#include "update_smoothing_length_by_band.h"

#include "base_particles.hpp"
#include "particle_iterators.h"

namespace SPH
{
//=================================================================================================//
UpdateSmoothingLengthByBand::UpdateSmoothingLengthByBand(
    RealBody &real_body, BaseInnerRelation &inner_relation)
    : adaptation_(DynamicCast<ParticleBandAdaptation>(this, real_body.getSPHAdaptation())),
      particles_(real_body.getBaseParticles()),
      inner_configuration_(inner_relation.getRelation().inner_configuration_),
      h_ref_(adaptation_.ReferenceSmoothingLength()),
      n_r_(adaptation_.ReferenceNeighborNumber()),
      h_spacing_ratio_(adaptation_.SmoothingLengthSpacingRatio()),
      h_ratio_(particles_.getVariableDataByName<Real>("SmoothingLengthRatio")),
      Vol_(particles_.getVariableDataByName<Real>("VolumetricMeasure")) {}
//=================================================================================================//
void UpdateSmoothingLengthByBand::exec()
{
    size_t total_real_particles = particles_.TotalRealParticles();
    h_new_.resize(total_real_particles);

    // pass 1: statistics from the previous smoothing lengths
    particle_for(execution::ParallelPolicy(), IndexRange(0, total_real_particles),
                 [&](size_t index_i)
                 {
                     Real h_old = h_ref_ / h_ratio_[index_i];
                     Real h_r = h_spacing_ratio_ * pow(Vol_[index_i], 1.0 / Real(Dimensions));
                     Real h_lower = 0.5 * h_r;
                     Real h_upper = SMIN(1.5 * h_r, h_ref_);

                     Neighborhood &neighborhood = inner_configuration_[index_i];
                     int n_i = 0;
                     Real sum_h = 0.0;
                     for (size_t n = 0; n != neighborhood.current_size_; ++n)
                     {
                         // W_ij is cached with i-side truncation, i.e. it is
                         // positive exactly within the support of h_i
                         if (neighborhood.W_ij_[n] > 0.0)
                         {
                             n_i++;
                             sum_h += h_ref_ / h_ratio_[neighborhood.j_[n]];
                         }
                     }

                     Real h_new = h_old;
                     if (n_i > 0)
                     {
                         Real h_tilde = 0.5 * h_old * (1.0 + sqrt(n_r_ / Real(n_i)));
                         h_tilde = SMAX(h_lower, SMIN(h_tilde, h_upper));
                         Real h_bar = sum_h / Real(n_i);
                         h_bar = SMAX(h_lower, SMIN(h_bar, h_upper));
                         h_new = 0.5 * (h_tilde + h_bar);
                     }
                     h_new_[index_i] = SMAX(h_lower, SMIN(h_new, h_upper));
                 });

    // pass 2: commit h_ratio >= 1 (h_new <= h_ref by the clamp above)
    particle_for(execution::ParallelPolicy(), IndexRange(0, total_real_particles),
                 [&](size_t index_i)
                 {
                     h_ratio_[index_i] = h_ref_ / h_new_[index_i];
                 });
}
//=================================================================================================//
} // namespace SPH
