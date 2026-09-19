/**
 * @file 	shepard_density_filter.cpp
 * @author 	KIYOYOZU
 */

#include "shepard_density_filter.h"

#include "base_body.h"
#include "base_particles.hpp"
#include "particle_iterators.h"

namespace SPH
{
//=================================================================================================//
ShepardDensityFilter::ShepardDensityFilter(
	RealBody &real_body, BaseInnerRelation &inner_relation, MultiphaseMixture &mixture)
	: particles_(real_body.getBaseParticles()),
	  inner_configuration_(inner_relation.getRelation().inner_configuration_),
	  mixture_(mixture),
	  Vol_(particles_.getVariableDataByName<Real>("VolumetricMeasure")),
	  mass_(particles_.getVariableDataByName<Real>("Mass")),
	  rho_(particles_.getVariableDataByName<Real>("Density")),
	  p_(particles_.getVariableDataByName<Real>("Pressure")),
	  E_(particles_.getVariableDataByName<Real>("TotalEnergy")),
	  alpha_(particles_.getVariableDataByName<Real>("VolumeFraction")),
	  vel_(particles_.getVariableDataByName<Vecd>("Velocity"))
{
	particles_.registerStateVariableData<Real>("DensityRaw");
	particles_.registerStateVariableData<Real>("PressureRaw");
	rho_raw_ = particles_.getVariableDataByName<Real>("DensityRaw");
	p_raw_ = particles_.getVariableDataByName<Real>("PressureRaw");
}
//=================================================================================================//
void ShepardDensityFilter::exec()
{
	size_t total_real_particles = particles_.TotalRealParticles();
	particle_for(execution::ParallelPolicy(), IndexRange(0, total_real_particles),
				 [&](size_t index_i)
				 {
					 rho_raw_[index_i] = rho_[index_i];
					 p_raw_[index_i] = p_[index_i];

					 // same-phase subset: exclude neighbors across the
					 // alpha = 0.5 material split; W_ij is cached with i-side
					 // truncation, so W_ij > 0 selects the support of h_i
					 bool alpha_side_i = alpha_[index_i] >= 0.5;
					 Real sum_mass_W = 0.0;
					 Real sum_vol_W = 0.0;
					 Neighborhood &neighborhood = inner_configuration_[index_i];
					 for (size_t n = 0; n != neighborhood.current_size_; ++n)
					 {
						 Real W_ij = neighborhood.W_ij_[n];
						 if (W_ij <= 0.0)
							 continue;
						 size_t index_j = neighborhood.j_[n];
						 if ((alpha_[index_j] >= 0.5) != alpha_side_i)
							 continue;
						 sum_mass_W += mass_[index_j] * W_ij;
						 sum_vol_W += Vol_[index_j] * W_ij;
					 }

					 if (sum_vol_W > TinyReal)
					 {
						 Real rho_hat = sum_mass_W / sum_vol_W;
						 rho_[index_i] = rho_hat;
						 Real rho_e = E_[index_i] / Vol_[index_i] -
									  0.5 * rho_hat * vel_[index_i].squaredNorm();
						 p_[index_i] = mixture_.MixturePressure(alpha_[index_i], rho_e);
					 }
				 });
}
//=================================================================================================//
} // namespace SPH
