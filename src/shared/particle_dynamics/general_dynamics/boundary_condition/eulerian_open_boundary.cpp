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
 * @file 	eulerian_open_boundary.cpp
 * @brief 	Implementation of the weakly compressible Eulerian open boundary
 * 			state synchronization helpers.
 * @author 	KIYOYOZU, Xiangyu Hu
 */
#include "eulerian_open_boundary.h"

#include "base_particles.hpp"
namespace SPH
{
//=================================================================================================//
EulerianWeaklyCompressibleBoundaryState
makeEulerianWeaklyCompressibleBoundaryState(BaseParticles &particles)
{
    EulerianWeaklyCompressibleBoundaryState state;
    state.rho_ = particles.getVariableDataByName<Real>("Density");
    state.p_ = particles.getVariableDataByName<Real>("Pressure");
    state.mass_ = particles.getVariableDataByName<Real>("Mass");
    state.Vol_ = particles.getVariableDataByName<Real>("VolumetricMeasure");
    state.vel_ = particles.getVariableDataByName<Vecd>("Velocity");
    state.mom_ = particles.getVariableDataByName<Vecd>("Momentum");
    // The owner argument (&particles) is only used for diagnostic type names
    // inside DynamicCast; the cast target is the matter material of the body.
    WeaklyCompressibleFluid &fluid =
        DynamicCast<WeaklyCompressibleFluid>(&particles, particles.getSPHBody().getMatterMaterial());
    state.fluid_ = &fluid;
    return state;
}
//=================================================================================================//
void syncEulerianWeaklyCompressibleState(EulerianWeaklyCompressibleBoundaryState &state, size_t index_i)
{
    // Mass = rho * Vol. The Eulerian integrator does not move particles, so Vol
    // is treated as fixed here. No division is involved, so a zero Vol yields
    // zero Mass rather than a NaN; the boundary layer is expected to carry real
    // volumes, and TinyReal is therefore not introduced to avoid biasing Mass.
    state.mass_[index_i] = state.rho_[index_i] * state.Vol_[index_i];
    // Momentum = Mass * vel, matching the integrator relation vel = mom / mass.
    state.mom_[index_i] = state.mass_[index_i] * state.vel_[index_i];
    // Pressure from the linear EOS p0_ * (rho / rho0_ - 1).
    state.p_[index_i] = state.fluid_->getPressure(state.rho_[index_i]);
}
//=================================================================================================//
void setEulerianWeaklyCompressibleVelocity(EulerianWeaklyCompressibleBoundaryState &state,
                                           size_t index_i, const Vecd &target_velocity)
{
    state.vel_[index_i] = target_velocity;
    syncEulerianWeaklyCompressibleState(state, index_i);
}
//=================================================================================================//
} // namespace SPH
