/**
 * @file 	gradient_correction.cpp
 * @author 	KIYOYOZU
 */

#include "gradient_correction.h"

#include "base_body.h"
#include "base_particles.hpp"
#include "particle_iterators.h"

namespace SPH
{
//=================================================================================================//
ComputeGradientCorrection::ComputeGradientCorrection(
	BaseInnerRelation &inner_relation, BaseContactRelation *wall_contact_relation)
	: particles_(inner_relation.getRelation().real_body_->getBaseParticles()),
	  inner_configuration_(inner_relation.getRelation().inner_configuration_),
	  wall_contact_(wall_contact_relation),
	  Vol_(particles_.getVariableDataByName<Real>("VolumetricMeasure")),
	  grad_corr_(particles_.registerStateVariableData<Vecd>("GradientCorrection")) {}
//=================================================================================================//
void ComputeGradientCorrection::exec()
{
	size_t total_real_particles = particles_.TotalRealParticles();

	// wall contact configuration is optional; when present its particles carry
	// only VolumetricMeasure (the mirror stencil), so only Vol is gathered
	StdVec<Real *> wall_Vol;
	StdVec<ParticleConfiguration> *contact_configuration = nullptr;
	if (wall_contact_ != nullptr)
	{
		StdVec<BaseParticles *> contact_particles = wall_contact_->getContactParticles();
		contact_configuration = &wall_contact_->contact_configuration_;
		for (size_t k = 0; k != contact_particles.size(); ++k)
			wall_Vol.push_back(contact_particles[k]->getVariableDataByName<Real>("VolumetricMeasure"));
	}

	particle_for(execution::ParallelPolicy(), IndexRange(0, total_real_particles),
				 [&](size_t index_i)
				 {
					 Vecd moment = Vecd::Zero();
					 Real vol_sum = 0.0;
					 Neighborhood &inner_neighborhood = inner_configuration_[index_i];
					 for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
					 {
						 Real V_j = Vol_[inner_neighborhood.j_[n]];
						 moment += V_j * inner_neighborhood.dW_ij_[n] * inner_neighborhood.e_ij_[n];
						 vol_sum += V_j;
					 }
					 if (contact_configuration != nullptr)
					 {
						 for (size_t k = 0; k != contact_configuration->size(); ++k)
						 {
							 Neighborhood &wall_neighborhood = (*contact_configuration)[k][index_i];
							 for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
							 {
								 Real V_j = wall_Vol[k][wall_neighborhood.j_[n]];
								 moment += V_j * wall_neighborhood.dW_ij_[n] * wall_neighborhood.e_ij_[n];
								 vol_sum += V_j;
							 }
						 }
					 }
					 if (vol_sum > TinyReal)
						 grad_corr_[index_i] = moment / vol_sum;
					 else
						 grad_corr_[index_i] = Vecd::Zero();
				 });
}
//=================================================================================================//
} // namespace SPH
