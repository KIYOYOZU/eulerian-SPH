/**
 * @file 	update_particle_bands.cpp
 * @author 	KIYOYOZU
 */

#include "update_particle_bands.h"

#include "base_body.h"
#include "base_particles.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <queue>

namespace SPH
{
//=================================================================================================//
UpdateParticleBands::UpdateParticleBands(
	RealBody &real_body, BaseInnerRelation &inner_relation,
	Real interface_alpha_tol, Real band_hysteresis, bool force_dijkstra,
	bool shock_band, Real shock_rel_jump)
	: adaptation_(DynamicCast<ParticleBandAdaptation>(this, real_body.getSPHAdaptation())),
	  particles_(real_body.getBaseParticles()),
	  inner_configuration_(inner_relation.getRelation().inner_configuration_),
	  alpha_tol_(interface_alpha_tol), band_hysteresis_(band_hysteresis),
	  force_dijkstra_(force_dijkstra),
	  shock_band_(shock_band), shock_rel_jump_(shock_rel_jump),
	  detected_interface_(std::numeric_limits<Real>::quiet_NaN()),
	  detected_shock_(std::numeric_limits<Real>::quiet_NaN()),
	  alpha_(particles_.getVariableDataByName<Real>("VolumeFraction")),
	  Vol_(particles_.getVariableDataByName<Real>("VolumetricMeasure")),
	  ref_spacing_(particles_.getVariableDataByName<Real>("ReferenceSpacing")),
	  p_(particles_.getVariableDataByName<Real>("Pressure")),
	  band_(particles_.getVariableDataByName<int>("ParticleBand")),
	  pos_(particles_.getVariableDataByName<Vecd>("Position")),
	  vel_(particles_.getVariableDataByName<Vecd>("Velocity")) {}
//=================================================================================================//
size_t UpdateParticleBands::exec()
{
	size_t total_real_particles = particles_.TotalRealParticles();
	detected_interface_ = std::numeric_limits<Real>::quiet_NaN();
	detected_shock_ = std::numeric_limits<Real>::quiet_NaN();
	if (!force_dijkstra_)
		detected_interface_ = detectPlanarInterface(total_real_particles);

	// shock tracking: one multi-source Dijkstra from the union of interface
	// and shock seeds, so every particle is banded by the graph distance to
	// the nearest feature -- geometry-free and 3D-ready, and no feature can
	// starve another (unlike an argmax over column jumps)
	if (shock_band_)
	{
		markShockSeeds();
		return assignByDijkstra(total_real_particles);
	}

	size_t changes = 0;
	if (detected_interface_ == detected_interface_) // planar path (not NaN)
	{
		for (size_t index_i = 0; index_i != total_real_particles; ++index_i)
			if (assignBand(index_i, std::abs(pos_[index_i][0] - detected_interface_)))
				changes++;
	}
	else
	{
		changes = assignByDijkstra(total_real_particles);
	}
	return changes;
}
//=================================================================================================//
Real UpdateParticleBands::detectPlanarInterface(size_t total_real_particles) const
{
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	if (total_real_particles == 0)
		return nan;

	Real x_min = MaxReal, x_max = -MaxReal;
	for (size_t index_i = 0; index_i != total_real_particles; ++index_i)
	{
		x_min = SMIN(x_min, pos_[index_i][0]);
		x_max = SMAX(x_max, pos_[index_i][0]);
	}

	// column sums of the phase indicator weighted by EQUAL phase volumes
	// (alpha_i V_i vs (1-alpha_i) V_i): this counts the gas and water
	// volumes per column directly, so the crossing tracks the alpha = 0.5
	// contour. A (2 alpha - 1) V weighting instead balances at
	// alpha ~ 0.92 for a 1.4/1000 density pair -- the smeared profile never
	// crosses it and the detected interface stalls while the true one moves.
	Real bin_width = 0.5 * adaptation_.FinestSpacing();
	int n_bins = SMAX(2, (int)std::ceil((x_max - x_min) / bin_width));
	StdVec<Real> bin_gas(n_bins, 0.0), bin_water(n_bins, 0.0);
	for (size_t index_i = 0; index_i != total_real_particles; ++index_i)
	{
		int b = (int)((pos_[index_i][0] - x_min) / bin_width);
		b = SMAX(0, SMIN(n_bins - 1, b));
		bin_gas[b] += alpha_[index_i] * Vol_[index_i];
		bin_water[b] += (1.0 - alpha_[index_i]) * Vol_[index_i];
	}

	int crossings = 0;
	int prev_bin = -1;
	Real bin_sum_prev = 0.0;
	Real x_if = nan;
	for (int b = 0; b != n_bins; ++b)
	{
		// S = gas volume - water volume per column: >0 on the gas side, <0 on
		// the water side; empty or exactly balanced bins carry no side
		// information and must not register as crossings (columns in the
		// coarse bands are wider than bin_width, leaving empty bins between)
		Real s_b = bin_gas[b] - bin_water[b];
		if (s_b == 0.0)
			continue;
		if (prev_bin >= 0 && (s_b > 0.0) != (bin_sum_prev > 0.0))
		{
			crossings++;
			Real x_a = x_min + (Real(prev_bin) + 0.5) * bin_width;
			Real x_b = x_min + (Real(b) + 0.5) * bin_width;
			x_if = x_a + (x_b - x_a) * bin_sum_prev / (bin_sum_prev - s_b);
		}
		bin_sum_prev = s_b;
		prev_bin = b;
	}
	return crossings == 1 ? x_if : nan;
}
//=================================================================================================//
size_t UpdateParticleBands::markShockSeeds()
{
	const size_t n = particles_.TotalRealParticles();
	shock_seed_.assign(n, 0);
	detected_shock_ = std::numeric_limits<Real>::quiet_NaN();

	size_t seed_count = 0;
	Real seed_x_sum = 0.0;
	for (size_t index_i = 0; index_i != n; ++index_i)
	{
		// gate 1: the anchor particle must be in a pure phase, otherwise the
		// pressure difference to its neighbors probes the material interface
		bool gas_i = alpha_[index_i] >= 1.0 - alpha_tol_;
		bool wat_i = alpha_[index_i] <= alpha_tol_;
		if (!gas_i && !wat_i)
			continue;

		Neighborhood &neighborhood = inner_configuration_[index_i];
		Real div_u = 0.0; // kernel-gradient velocity divergence
		Real max_rel_jump = 0.0;
		for (size_t k = 0; k != neighborhood.current_size_; ++k)
		{
			size_t index_j = neighborhood.j_[k];
			if (index_j >= n)
				continue; // wall ghosts carry mirrored states, not features
			// same-pure-phase pairs only: pressure is continuous across the
			// material interface, so cross-phase pairs carry no shock signal
			bool same_pure = (gas_i && alpha_[index_j] >= 1.0 - alpha_tol_) ||
							 (wat_i && alpha_[index_j] <= alpha_tol_);
			if (!same_pure)
				continue;

			Real w = Vol_[index_j] * neighborhood.dW_ij_[k];
			Vecd du = vel_[index_j] - vel_[index_i];
			div_u += w * du.dot(neighborhood.e_ij_[k]);
			Real dp = p_[index_j] - p_[index_i];
			max_rel_jump = SMAX(max_rel_jump,
								std::abs(dp) / SMAX(SMAX(std::abs(p_[index_i]), std::abs(p_[index_j])), TinyReal));
		}
		// gate 2 (compression): div u < 0 keeps shocks, excludes rarefactions
		// (the difference form is exactly zero for a uniform velocity field);
		// gate 3 (magnitude): a pair pressure jump above shock_rel_jump of the
		// local level excludes pure-noise fields. A Jameson monotonicity gate
		// was tried and removed: it rejects smeared shocks (signed sums cancel
		// on a monotone ramp) while passing oscillation peaks (one-sided), so
		// it discriminated exactly backwards here. Ringing particles can pass
		// the remaining gates; that only refines a noisy region and is benign,
		// since seeds never compete with each other.
		bool is_seed = div_u < 0.0 && max_rel_jump >= shock_rel_jump_;
		if (is_seed)
		{
			shock_seed_[index_i] = 1;
			seed_count++;
			seed_x_sum += pos_[index_i][0];
		}
	}
	if (seed_count > 0)
		detected_shock_ = seed_x_sum / Real(seed_count); // log/monitor only
	return seed_count;
}
//=================================================================================================//
size_t UpdateParticleBands::assignByDijkstra(size_t total_real_particles)
{
	StdVec<Real> distance(total_real_particles, MaxReal);
	typedef std::pair<Real, size_t> Entry;
	std::priority_queue<Entry, StdVec<Entry>, std::greater<Entry>> queue;

	// interface seeds: alpha within tol of the mixed state, or a neighbor
	// alpha jump across the material split; shock seeds come from
	// markShockSeeds (empty when shock_band is off)
	for (size_t index_i = 0; index_i != total_real_particles; ++index_i)
	{
		bool is_seed = alpha_[index_i] > alpha_tol_ &&
					   alpha_[index_i] < 1.0 - alpha_tol_;
		if (shock_band_ && index_i < shock_seed_.size() && shock_seed_[index_i])
			is_seed = true;
		if (!is_seed)
		{
			Neighborhood &neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != neighborhood.current_size_; ++n)
			{
				size_t index_j = neighborhood.j_[n];
				if (index_j < total_real_particles &&
					std::abs(alpha_[index_i] - alpha_[index_j]) > 0.5)
				{
					is_seed = true;
					break;
				}
			}
		}
		if (is_seed)
		{
			distance[index_i] = 0.0;
			queue.push(Entry(0.0, index_i));
		}
	}

	if (queue.empty())
	{
		// warn once: single-phase runs legitimately have no interface and this
		// path would otherwise print every adapt_interval for the whole run
		if (!warned_no_seeds_)
		{
			std::cout << "\n Warning: no interface/shock seeds found, particle bands unchanged."
					  << std::endl;
			warned_no_seeds_ = true;
		}
		return 0;
	}

	while (!queue.empty())
	{
		Entry top = queue.top();
		queue.pop();
		Real d = top.first;
		size_t index_i = top.second;
		if (d > distance[index_i])
			continue; // stale queue entry
		Neighborhood &neighborhood = inner_configuration_[index_i];
		for (size_t n = 0; n != neighborhood.current_size_; ++n)
		{
			size_t index_j = neighborhood.j_[n];
			if (index_j >= total_real_particles)
				continue;
			Real candidate = d + neighborhood.r_ij_[n];
			if (candidate < distance[index_j])
			{
				distance[index_j] = candidate;
				queue.push(Entry(candidate, index_j));
			}
		}
	}

	size_t changes = 0;
	for (size_t index_i = 0; index_i != total_real_particles; ++index_i)
		if (distance[index_i] < MaxReal)
			if (assignBand(index_i, distance[index_i]))
				changes++;
	return changes;
}
//=================================================================================================//
bool UpdateParticleBands::assignBand(size_t index_i, Real distance)
{
	int target = adaptation_.BandOfDistance(distance);
	if (target == band_[index_i])
		return false;

	// hysteresis: only commit a band change once the distance has crossed the
	// relevant band boundary by a margin, so boundary jitter does not flip the
	// spacing target (and with it gamma = Vol / ds_band^2) every update
	Real margin = band_hysteresis_ * adaptation_.BandSpacing(target);
	if (target > band_[index_i])
	{
		if (distance < adaptation_.BandBoundary(target) + margin)
			return false;
	}
	else
	{
		if (distance > adaptation_.BandBoundary(band_[index_i]) - margin)
			return false;
	}

	band_[index_i] = target;
	ref_spacing_[index_i] = adaptation_.BandSpacing(target);
	return true;
}
//=================================================================================================//
} // namespace SPH
