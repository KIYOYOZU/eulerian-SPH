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
 * @file    multiphase_mixture.h
 * @brief   Mixture model for the Kapila five-equation two-phase flow model.
 *          Both phases obey the stiffened gas EOS; mixture pressure closure
 *          assumes instantaneous pressure equilibrium (p1 = p2 = p).
 *          Mixture gamma and p_inf are volume-fraction weighted.
 * @author  KIYOYOZU
 */

#ifndef MULTIPHASE_MIXTURE_H
#define MULTIPHASE_MIXTURE_H

#include "stiffened_gas.h"

namespace SPH
{
/**
 * @class MultiphaseMixture
 * @brief Two-phase mixture for the Kapila five-equation model.
 *
 *  Conservation variables: alpha1*rho1, alpha2*rho2, rho*u, rho*E, alpha1.
 *  Mixture density:  rho = alpha1*rho1 + alpha2*rho2.
 *
 *  Mixture EOS parameters (isobaric closure):
 *    1/(gamma_mix - 1) = alpha1/(gamma1 - 1) + alpha2/(gamma2 - 1)
 *    gamma_mix*p_inf_mix/(gamma_mix - 1)
 *        = alpha1*gamma1*p_inf1/(gamma1 - 1) + alpha2*gamma2*p_inf2/(gamma2 - 1)
 *
 *  Mixture pressure from internal energy per unit volume (rho_e):
 *    p = (rho_e - B) / A
 *    where A = alpha1/(gamma1-1) + alpha2/(gamma2-1)
 *          B = alpha1*gamma1*p_inf1/(gamma1-1) + alpha2*gamma2*p_inf2/(gamma2-1)
 *
 *  Wood sound speed:
 *    1/(rho*c^2) = alpha1/(gamma1*(p+p_inf1)) + alpha2/(gamma2*(p+p_inf2))
 */
class MultiphaseMixture
{
  protected:
    StiffenedGas &phase_1_; /**< phase 1 material (e.g. gas) */
    StiffenedGas &phase_2_; /**< phase 2 material (e.g. liquid) */

  public:
    MultiphaseMixture(StiffenedGas &phase_1, StiffenedGas &phase_2);
    virtual ~MultiphaseMixture() = default;

    StiffenedGas &GetPhase1() { return phase_1_; };
    StiffenedGas &GetPhase2() { return phase_2_; };

    /**
     * Helper: A(alpha1) = alpha1/(gamma1-1) + (1-alpha1)/(gamma2-1).
     * This equals 1/(gamma_mix - 1).
     */
    Real MixtureA(Real alpha1) const;

    /**
     * Helper: B(alpha1) = alpha1*gamma1*p_inf1/(gamma1-1)
     *                    + (1-alpha1)*gamma2*p_inf2/(gamma2-1).
     * This equals gamma_mix*p_inf_mix/(gamma_mix - 1).
     */
    Real MixtureB(Real alpha1) const;

    /** Mixture heat capacity ratio from volume fraction alpha1. */
    Real MixtureGamma(Real alpha1) const;

    /** Mixture stiffness pressure from volume fraction alpha1. */
    Real MixturePInf(Real alpha1) const;

    /**
     * Mixture pressure from internal energy per unit volume.
     * p = (rho_e - B) / A
     */
    Real MixturePressure(Real alpha1, Real rho_e) const;

    /**
     * Mixture internal energy per unit volume from pressure.
     * rho_e = A * p + B
     */
    Real MixtureInternalEnergyPerVolume(Real alpha1, Real p) const;

    /**
     * Wood mixture sound speed.
     * 1/(rho*c^2) = alpha1/(gamma1*(p+p_inf1)) + alpha2/(gamma2*(p+p_inf2))
     * @param alpha1 volume fraction of phase 1
     * @param rho mixture density
     * @param p mixture pressure
     */
    Real MixtureSoundSpeed(Real alpha1, Real rho, Real p) const;
};
} // namespace SPH

#endif // MULTIPHASE_MIXTURE_H
