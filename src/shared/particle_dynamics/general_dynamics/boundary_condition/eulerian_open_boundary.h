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
 * @file 	eulerian_open_boundary.h
 * @brief 	Consistent state synchronization helpers for weakly compressible
 * 			Eulerian SPH open boundaries. The Eulerian weakly compressible
 * 			integrator keeps Velocity = Momentum / Mass internally, so a
 * 			boundary condition that only overwrites Velocity leaves the
 * 			Momentum/Mass pair stale and the next integration step splits.
 * 			These helpers keep the rho/p/mass/vel/mom quintuple coherent.
 * @author 	KIYOYOZU, Xiangyu Hu
 */
#ifndef EULERIAN_OPEN_BOUNDARY_H
#define EULERIAN_OPEN_BOUNDARY_H

#include "base_general_dynamics.h"
#include "weakly_compressible_fluid.h"

namespace SPH
{
/**
 * @struct EulerianWeaklyCompressibleBoundaryState
 * @brief Aggregate of the particle field pointers and EOS handle needed by a
 *        weakly compressible Eulerian boundary update.
 *        Extracted once from BaseParticles to avoid repeated name lookups
 *        during boundary application.
 */
struct EulerianWeaklyCompressibleBoundaryState
{
    Real *rho_;                      /**< Density */
    Real *p_;                        /**< Pressure */
    Real *mass_;                     /**< Mass */
    Real *Vol_;                      /**< VolumetricMeasure */
    Vecd *vel_;                      /**< Velocity */
    Vecd *mom_;                      /**< Momentum */
    WeaklyCompressibleFluid *fluid_; /**< EOS handle for pressure recovery */
};

/**
 * @brief Build the state aggregate from BaseParticles and its base material.
 *        The base material is dynamically cast to WeaklyCompressibleFluid so
 *        that the linear EOS is available for pressure recovery.
 * @param particles The fluid particles carrying rho/p/mass/Vol/vel/mom.
 * @return Aggregated state; fields point directly into the particle storage.
 */
EulerianWeaklyCompressibleBoundaryState
makeEulerianWeaklyCompressibleBoundaryState(BaseParticles &particles);

/**
 * @brief Consistent synchronization of the boundary quintuple.
 *        Recomputes Mass = rho * Vol, Momentum = Mass * vel, and
 *        p = EOS(rho) from the current rho_ and vel_ values.
 *        Idempotent: the function only reads rho_ and vel_, so applying it
 *        twice on an unchanged state yields identical Mass/Momentum/Pressure.
 * @param state   Aggregated field pointers and EOS handle.
 * @param index_i Target particle index.
 */
void syncEulerianWeaklyCompressibleState(EulerianWeaklyCompressibleBoundaryState &state, size_t index_i);

/**
 * @brief Impose a target velocity at a boundary particle, then run sync so
 *        that Mass/Momentum/Pressure stay consistent with the new velocity.
 * @param state           Aggregated field pointers and EOS handle.
 * @param index_i         Target particle index.
 * @param target_velocity Velocity to write into vel_[index_i].
 */
void setEulerianWeaklyCompressibleVelocity(EulerianWeaklyCompressibleBoundaryState &state,
                                           size_t index_i, const Vecd &target_velocity);
} // namespace SPH
#endif // EULERIAN_OPEN_BOUNDARY_H
