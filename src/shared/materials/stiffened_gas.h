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
 * @file    stiffened_gas.h
 * @brief   Stiffened gas equation of state for compressible multiphase flows.
 *          EOS: p = (gamma - 1) * rho * e - gamma * p_inf
 *          Sound speed: c = sqrt(gamma * (p + p_inf) / rho)
 * @author  KIYOYOZU
 */

#ifndef STIFFENED_GAS_H
#define STIFFENED_GAS_H

#include "base_material.h"

namespace SPH
{
/**
 * @class StiffenedGas
 * @brief Stiffened gas equation of state.
 *        Widely used for modeling water and other liquids in
 *        compressible multiphase flow simulations.
 */
class StiffenedGas : public Fluid
{
  protected:
    Real gamma_;  /**< heat capacity ratio */
    Real p_inf_;  /**< reference (stiffness) pressure */

  public:
    explicit StiffenedGas(Real gamma, Real p_inf);
    virtual ~StiffenedGas();

    virtual Real ReferenceDensity() const override { return 1.0; };
    virtual Real ReferenceSoundSpeed() const override { return 1.0; };

    Real HeatCapacityRatio() const { return gamma_; };
    Real ReferencePressure() const { return p_inf_; };

    /** Pressure from density and internal energy per unit volume:
     *  p = (gamma - 1) * rho_e - gamma * p_inf */
    virtual Real getPressure(Real rho, Real rho_e) override;
    virtual Real getPressure(Real rho) override { return 0.0; };
    virtual Real DensityFromPressure(Real p) override { return 0.0; };

    /** Sound speed: c = sqrt(gamma * (p + p_inf) / rho) */
    virtual Real getSoundSpeed(Real p, Real rho) override;

    /** Internal energy per unit volume from density and pressure:
     *  rho_e = (p + gamma * p_inf) / (gamma - 1) */
    Real InternalEnergyPerVolume(Real rho, Real p) const;
};
} // namespace SPH

#endif // STIFFENED_GAS_H
