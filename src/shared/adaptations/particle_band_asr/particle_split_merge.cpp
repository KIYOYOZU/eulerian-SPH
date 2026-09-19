/**
 * @file 	particle_split_merge.cpp
 * @author 	KIYOYOZU
 */

#include "particle_split_merge.h"

#include "particle_operation.hpp"

#include "base_body.h"
#include "base_particles.hpp"
#include "particle_iterators.h"

#include <algorithm>
#include <iostream>
#include <limits>

namespace SPH
{
//=================================================================================================//
ParticleSplittingByBand::ParticleSplittingByBand(
	RealBody &real_body, BaseInnerRelation &inner_relation,
	Real split_threshold, Real offset_factor, Real periodic_height)
	: adaptation_(DynamicCast<ParticleBandAdaptation>(this, real_body.getSPHAdaptation())),
	  particles_(real_body.getBaseParticles()),
	  inner_configuration_(inner_relation.getRelation().inner_configuration_),
	  gamma_split_(split_threshold), lambda_(offset_factor),
	  periodic_height_(periodic_height),
	  pos_(particles_.getVariableDataByName<Vecd>("Position")),
	  vel_(particles_.getVariableDataByName<Vecd>("Velocity")),
	  mom_(particles_.getVariableDataByName<Vecd>("Momentum")),
	  Vol_(particles_.getVariableDataByName<Real>("VolumetricMeasure")),
	  mass_(particles_.getVariableDataByName<Real>("Mass")),
	  rho_(particles_.getVariableDataByName<Real>("Density")),
	  p_(particles_.getVariableDataByName<Real>("Pressure")),
	  E_(particles_.getVariableDataByName<Real>("TotalEnergy")),
	  alpha_(particles_.getVariableDataByName<Real>("VolumeFraction")),
	  band_(particles_.getVariableDataByName<int>("ParticleBand")),
	  h_ratio_(particles_.getVariableDataByName<Real>("SmoothingLengthRatio")),
	  ref_spacing_(particles_.getVariableDataByName<Real>("ReferenceSpacing")),
	  spawn_(&particles_) {}
//=================================================================================================//
size_t ParticleSplittingByBand::findNearestNeighbor(size_t index_i) const
{
	const Neighborhood &neighborhood = inner_configuration_[index_i];
	size_t nearest = std::numeric_limits<size_t>::max();
	Real r_min = MaxReal;
	for (size_t n = 0; n != neighborhood.current_size_; ++n)
	{
		if (neighborhood.r_ij_[n] < r_min)
		{
			r_min = neighborhood.r_ij_[n];
			nearest = neighborhood.j_[n];
		}
	}
	return nearest;
}
//=================================================================================================//
size_t ParticleSplittingByBand::exec()
{
	size_t total_real_particles = particles_.TotalRealParticles();
	UnsignedInt particles_bound = particles_.ParticlesBound();

	// frozen candidate list: spawn only appends at the end, so mother and
	// nearest-neighbor indices stay valid throughout the pass
	StdVec<std::pair<size_t, size_t>> candidates;
	for (size_t index_i = 0; index_i != total_real_particles; ++index_i)
	{
		Real ds_band = adaptation_.BandSpacing(band_[index_i]);
		Real gamma = Vol_[index_i] / (ds_band * ds_band);
		if (gamma > gamma_split_)
		{
			size_t index_j = findNearestNeighbor(index_i);
			if (index_j != std::numeric_limits<size_t>::max())
				candidates.push_back({index_i, index_j});
		}
	}

	// capacity pre-check: SpawnRealParticle increments first, checks after
	if (total_real_particles + candidates.size() > particles_bound)
	{
		size_t allowed = particles_bound > total_real_particles
							 ? particles_bound - total_real_particles
							 : 0;
		std::cout << "\n Warning: particle buffer exhausted, "
				  << candidates.size() - allowed << " of " << candidates.size()
				  << " split candidates skipped." << std::endl;
		candidates.resize(allowed);
	}

	SpawnRealParticle::ComputingKernel spawn_kernel(execution::SequencedPolicy(), spawn_);
	size_t count = 0;
	for (auto &candidate : candidates)
	{
		size_t index_i = candidate.first;
		size_t index_j = candidate.second;
		UnsignedInt new_index_i = spawn_kernel(UnsignedInt(index_i));

		// offset perpendicular to the nearest-neighbor line, Eqs. (45)-(47);
		// the 2D perpendicular construction (3D splitting is out of scope)
		Vecd displacement = pos_[index_j] - pos_[index_i];
		Real norm = displacement.norm() + TinyReal;
		Vecd e(-displacement[1] / norm, displacement[0] / norm);
		Vecd offset = 0.5 * lambda_ * pow(Vol_[index_i], 1.0 / Real(Dimensions)) * e;

		pos_[new_index_i] = pos_[index_i] + offset;
		pos_[index_i] = pos_[index_i] - offset;
		wrapPeriodicY(pos_[new_index_i]);
		wrapPeriodicY(pos_[index_i]);

		// halve the extensive state (spawn copied the evolving state already);
		// intensive state is inherited: halving leaves rho, p and vel unchanged
		Real vol = 0.5 * Vol_[index_i];
		Vol_[index_i] = vol;
		Vol_[new_index_i] = vol;
		Real mass = 0.5 * mass_[index_i];
		mass_[index_i] = mass;
		mass_[new_index_i] = mass;
		Vecd mom = 0.5 * mom_[index_i];
		mom_[index_i] = mom;
		mom_[new_index_i] = mom;
		Real energy = 0.5 * E_[index_i];
		E_[index_i] = energy;
		E_[new_index_i] = energy;

		vel_[new_index_i] = vel_[index_i];
		rho_[new_index_i] = rho_[index_i];
		p_[new_index_i] = p_[index_i];
		alpha_[new_index_i] = alpha_[index_i];
		count++;
	}
	return count;
}
//=================================================================================================//
ParticleMergingByBand::ParticleMergingByBand(
	RealBody &real_body, BaseInnerRelation &inner_relation,
	MultiphaseMixture &mixture, Real merge_threshold)
	: adaptation_(DynamicCast<ParticleBandAdaptation>(this, real_body.getSPHAdaptation())),
	  particles_(real_body.getBaseParticles()),
	  inner_configuration_(inner_relation.getRelation().inner_configuration_),
	  mixture_(mixture), gamma_merge_(merge_threshold),
	  pos_(particles_.getVariableDataByName<Vecd>("Position")),
	  vel_(particles_.getVariableDataByName<Vecd>("Velocity")),
	  mom_(particles_.getVariableDataByName<Vecd>("Momentum")),
	  Vol_(particles_.getVariableDataByName<Real>("VolumetricMeasure")),
	  mass_(particles_.getVariableDataByName<Real>("Mass")),
	  rho_(particles_.getVariableDataByName<Real>("Density")),
	  p_(particles_.getVariableDataByName<Real>("Pressure")),
	  E_(particles_.getVariableDataByName<Real>("TotalEnergy")),
	  alpha_(particles_.getVariableDataByName<Real>("VolumeFraction")),
	  band_(particles_.getVariableDataByName<int>("ParticleBand")),
	  h_ratio_(particles_.getVariableDataByName<Real>("SmoothingLengthRatio")),
	  ref_spacing_(particles_.getVariableDataByName<Real>("ReferenceSpacing")),
	  remover_(&particles_),
	  life_status_(particles_.ParticlesBound(), 0) {}
//=================================================================================================//
size_t ParticleMergingByBand::exec()
{
	size_t total_real_particles = particles_.TotalRealParticles();
	std::fill(life_status_.begin(), life_status_.end(), 0);
	partner_.assign(total_real_particles, -1);

	// candidate scan on the frozen configuration: nearest neighbor must be
	// on the same side of the alpha = 0.5 material split and must not be a
	// band-0 (finest) interface seed
	for (size_t index_i = 0; index_i != total_real_particles; ++index_i)
	{
		Real ds_band = adaptation_.BandSpacing(band_[index_i]);
		if (Vol_[index_i] / (ds_band * ds_band) >= gamma_merge_)
			continue;

		const Neighborhood &neighborhood = inner_configuration_[index_i];
		bool alpha_side_i = alpha_[index_i] >= 0.5;
		Real r_min = MaxReal;
		int nearest = -1;
		for (size_t n = 0; n != neighborhood.current_size_; ++n)
		{
			size_t index_j = neighborhood.j_[n];
			if (band_[index_j] == 0)
				continue;
			if ((alpha_[index_j] >= 0.5) != alpha_side_i)
				continue;
			if (neighborhood.r_ij_[n] < r_min)
			{
				r_min = neighborhood.r_ij_[n];
				nearest = int(index_j);
			}
		}
		partner_[index_i] = nearest;
	}

	// mutual pairs merge once; the lower index wins and accumulates the
	// conserved state of the loser
	size_t count = 0;
	for (size_t index_i = 0; index_i != total_real_particles; ++index_i)
	{
		int partner = partner_[index_i];
		if (partner < 0 || size_t(partner) <= index_i ||
			partner_[size_t(partner)] != int(index_i))
			continue;
		size_t index_j = size_t(partner);

		Real vol = Vol_[index_i] + Vol_[index_j];
		Real alpha = (alpha_[index_i] * Vol_[index_i] + alpha_[index_j] * Vol_[index_j]) / vol;
		Real mass = mass_[index_i] + mass_[index_j];
		Vecd mom = mom_[index_i] + mom_[index_j];
		Real energy = E_[index_i] + E_[index_j];
		Vecd pos = (mass_[index_i] * pos_[index_i] + mass_[index_j] * pos_[index_j]) / mass;
		Vecd vel = mom / mass;
		Real rho = mass / vol;
		Real rho_e = energy / vol - 0.5 * rho * vel.squaredNorm();

		pos_[index_i] = pos;
		Vol_[index_i] = vol;
		mass_[index_i] = mass;
		mom_[index_i] = mom;
		E_[index_i] = energy;
		alpha_[index_i] = alpha;
		vel_[index_i] = vel;
		rho_[index_i] = rho;
		p_[index_i] = mixture_.MixturePressure(alpha, rho_e);
		ref_spacing_[index_i] = adaptation_.BandSpacing(band_[index_i]);

		life_status_[index_j] = 1;
		count++;
	}

	if (count == 0)
		return 0;

	// single compression pass: ascending order, the kernel consumes marked
	// tail slots itself and resets their life status, so no marked particle
	// is ever removed twice
	RemoveRealParticle::ComputingKernel remove_kernel(execution::SequencedPolicy(), remover_);
	for (size_t index_i = 0; index_i != total_real_particles; ++index_i)
	{
		if (life_status_[index_i] == 1)
			remove_kernel(UnsignedInt(index_i), life_status_.data());
	}

	// the swap-back relocated survivors carry stale non-evolving state
	recoverState(particles_.TotalRealParticles());
	return count;
}
//=================================================================================================//
void ParticleMergingByBand::recoverState(size_t total_particles)
{
	particle_for(execution::ParallelPolicy(), IndexRange(0, total_particles),
				 [&](size_t index_i)
				 {
					 vel_[index_i] = mom_[index_i] / mass_[index_i];
					 rho_[index_i] = mass_[index_i] / Vol_[index_i];
					 Real rho_e = E_[index_i] / Vol_[index_i] -
								  0.5 * rho_[index_i] * vel_[index_i].squaredNorm();
					 p_[index_i] = mixture_.MixturePressure(alpha_[index_i], rho_e);
				 });
}
//=================================================================================================//
} // namespace SPH
