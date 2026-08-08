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
 * @file    eulerian_multiphase_integration.h
 * @brief   Eulerian integration for the Kapila five-equation two-phase model:
 *          first-order (piecewise constant) and second-order (MUSCL-HLLC)
 *          variants. Interface states are computed in
 *          eulerian_multiphase_riemann_solver.h (first-order HLLC solver /
 *          MUSCL bridge with Wood sound speed).
 *
 *          Conservation variables per particle:
 *            - Mass (mixture), Momentum, TotalEnergy, VolumeFraction (alpha1)
 *          The mixture pressure is recovered from the mixture internal energy
 *          using MultiphaseMixture::MixturePressure.
 *
 *          Both Inner<> and Contact<Wall> specializations are provided so that
 *          reflective solid walls can be used at the tube ends, matching the
 *          single-phase shock-tube paradigm.
 * @author  KIYOYOZU
 */

#ifndef EULERIAN_MULTIPHASE_INTEGRATION_H
#define EULERIAN_MULTIPHASE_INTEGRATION_H

#include "base_general_dynamics.h"
#include "eulerian_multiphase_riemann_solver.h"
#include "multiphase_mixture.h"
#include "fluid_integration.hpp"
#include "fluid_time_step.h"
#include "muscl_reconstruction.hpp"

namespace SPH
{
namespace fluid_dynamics
{

/**
 * @class BaseIntegrationInMultiphase
 * @brief Base class (inner delegation) registering the extra state variables
 *        required by the five-equation model and holding a mixture copy.
 */
class BaseIntegrationInMultiphase : public BaseIntegration<DataDelegateInner>
{
  public:
    explicit BaseIntegrationInMultiphase(BaseInnerRelation &inner_relation,
                                           MultiphaseMixture &mixture);
    virtual ~BaseIntegrationInMultiphase() = default;

  protected:
    // Stored by value: the mixture reaches the integrator through DynamicsArgs,
    // which copies it into a temporary tuple. A reference to that tuple element
    // would dangle once construction ends. MultiphaseMixture is a small
    // stateless object (two phase references), so a copy is safe.
    MultiphaseMixture mixture_;
    Real *Vol_, *E_, *dE_dt_, *dmass_dt_;
    Real *alpha_, *dalpha_dt_;
    Vecd *mom_, *force_, *force_prior_;
};

/**
 * @class BaseIntegrationInMultiphaseForWall
 * @brief Base class (contact delegation) for the wall specializations.
 *        Registers the same extra state variables. The mixture reference is
 *        carried by the derived integration classes (mirroring how the
 *        single-phase compressible integrators carry their material).
 *
 *        MSVC-friendly adapter for InteractionWithWall: it must be a
 *        variadic class template (`template <typename...> class`) so that it
 *        matches InteractionWithWall's template-template parameter. A plain
 *        class (or a single-type-parameter template) triggers MSVC C3200.
 */
template <typename...>
class BaseIntegrationInMultiphaseForWall : public BaseIntegration<DataDelegateContact>
{
  public:
    explicit BaseIntegrationInMultiphaseForWall(BaseContactRelation &wall_contact_relation)
        : BaseIntegration<DataDelegateContact>(wall_contact_relation),
          Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
          E_(this->particles_->template registerStateVariableData<Real>("TotalEnergy")),
          dE_dt_(this->particles_->template registerStateVariableData<Real>("TotalEnergyChangeRate")),
          dmass_dt_(this->particles_->template registerStateVariableData<Real>("MassChangeRate")),
          alpha_(this->particles_->template registerStateVariableData<Real>("VolumeFraction")),
          dalpha_dt_(this->particles_->template registerStateVariableData<Real>("VolumeFractionChangeRate")),
          mom_(this->particles_->template registerStateVariableData<Vecd>("Momentum")),
          force_(this->particles_->template registerStateVariableData<Vecd>("Force")),
          force_prior_(this->particles_->template registerStateVariableData<Vecd>("ForcePrior"))
    {
    }
    virtual ~BaseIntegrationInMultiphaseForWall() = default;

  protected:
    Real *Vol_, *E_, *dE_dt_, *dmass_dt_;
    Real *alpha_, *dalpha_dt_;
    Vecd *mom_, *force_, *force_prior_;
};

//----------------------------------------------------------------------
//	First half step: momentum equation
//	d(rho*u)/dt = -div(rho*u (x) u + p*I)
//----------------------------------------------------------------------
template <typename... InteractionTypes>
class EulerianMultiphaseIntegration1stHalf;

template <>
class EulerianMultiphaseIntegration1stHalf<Inner<>> : public BaseIntegrationInMultiphase
{
  public:
    explicit EulerianMultiphaseIntegration1stHalf(BaseInnerRelation &inner_relation,
                                            MultiphaseMixture &mixture);

    template <typename BodyRelationType, typename FirstArg>
    explicit EulerianMultiphaseIntegration1stHalf(
        DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : EulerianMultiphaseIntegration1stHalf(
              parameters.identifier_, std::get<0>(parameters.others_))
    {
    }

    virtual ~EulerianMultiphaseIntegration1stHalf() = default;

    MultiphaseHLLCRiemannSolver riemann_solver_;

    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};

template <>
class EulerianMultiphaseIntegration1stHalf<Contact<Wall>>
    : public InteractionWithWall<BaseIntegrationInMultiphaseForWall>
{
  public:
    explicit EulerianMultiphaseIntegration1stHalf(BaseContactRelation &contact_relation,
                                            MultiphaseMixture &mixture);

    template <typename BodyRelationType, typename FirstArg>
    explicit EulerianMultiphaseIntegration1stHalf(
        DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : EulerianMultiphaseIntegration1stHalf(
              parameters.identifier_, std::get<0>(parameters.others_))
    {
    }

    virtual ~EulerianMultiphaseIntegration1stHalf() = default;

    MultiphaseHLLCRiemannSolver riemann_solver_;

    void interaction(size_t index_i, Real dt = 0.0);
};

//----------------------------------------------------------------------
//	Second half step: mass, energy, and volume fraction equations
//	d(rho)/dt   = -div(rho*u)
//	d(rho*E)/dt = -div((rho*E + p)*u)
//	d(alpha)/dt = -u . grad(alpha)   (advected with the mixture velocity)
//----------------------------------------------------------------------
template <typename... InteractionTypes>
class EulerianMultiphaseIntegration2ndHalf;

template <>
class EulerianMultiphaseIntegration2ndHalf<Inner<>> : public BaseIntegrationInMultiphase
{
  public:
    explicit EulerianMultiphaseIntegration2ndHalf(BaseInnerRelation &inner_relation,
                                            MultiphaseMixture &mixture);

    template <typename BodyRelationType, typename FirstArg>
    explicit EulerianMultiphaseIntegration2ndHalf(
        DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : EulerianMultiphaseIntegration2ndHalf(
              parameters.identifier_, std::get<0>(parameters.others_))
    {
    }

    virtual ~EulerianMultiphaseIntegration2ndHalf() = default;

    MultiphaseHLLCRiemannSolver riemann_solver_;

    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};

template <>
class EulerianMultiphaseIntegration2ndHalf<Contact<Wall>>
    : public InteractionWithWall<BaseIntegrationInMultiphaseForWall>
{
  public:
    explicit EulerianMultiphaseIntegration2ndHalf(BaseContactRelation &contact_relation,
                                            MultiphaseMixture &mixture);

    template <typename BodyRelationType, typename FirstArg>
    explicit EulerianMultiphaseIntegration2ndHalf(
        DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : EulerianMultiphaseIntegration2ndHalf(
              parameters.identifier_, std::get<0>(parameters.others_))
    {
    }

    virtual ~EulerianMultiphaseIntegration2ndHalf() = default;

    MultiphaseHLLCRiemannSolver riemann_solver_;

    void interaction(size_t index_i, Real dt = 0.0);
};

//----------------------------------------------------------------------
//	With-wall aliases (inner + wall contact).
//----------------------------------------------------------------------
using EulerianMultiphaseIntegration1stHalfWithWall =
    ComplexInteraction<EulerianMultiphaseIntegration1stHalf<Inner<>, Contact<Wall>>>;

using EulerianMultiphaseIntegration2ndHalfWithWall =
    ComplexInteraction<EulerianMultiphaseIntegration2ndHalf<Inner<>, Contact<Wall>>>;

//----------------------------------------------------------------------
//	Volume fraction gradient for the MUSCL reconstruction. Plain SPH
//	gradient (raw stencil), stored as "VolumeFractionGradient".
//----------------------------------------------------------------------
class VolumeFractionGradient : public LocalDynamics, public DataDelegateInner
{
  public:
    explicit VolumeFractionGradient(BaseInnerRelation &inner_relation);
    virtual ~VolumeFractionGradient() = default;

    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0) {}

  protected:
    Real *Vol_;
    Real *alpha_;
    Vecd *alpha_grad_;
};

//----------------------------------------------------------------------
//	Second-order (MUSCL-HLLC) variants. Interface states come from
//	MultiphaseMUSCLBridge (see eulerian_multiphase_riemann_solver.h).
//----------------------------------------------------------------------
template <typename... InteractionTypes>
class EulerianMultiphaseIntegration1stHalfMUSCL;

template <>
class EulerianMultiphaseIntegration1stHalfMUSCL<Inner<>> : public BaseIntegrationInMultiphase
{
  public:
    EulerianMultiphaseIntegration1stHalfMUSCL(BaseInnerRelation &inner_relation,
                                        MultiphaseMixture &mixture,
                                        const SecondOrderConfig &cfg);

    template <typename BodyRelationType, typename FirstArg, typename SecondArg>
    explicit EulerianMultiphaseIntegration1stHalfMUSCL(
        DynamicsArgs<BodyRelationType, FirstArg, SecondArg> parameters)
        : EulerianMultiphaseIntegration1stHalfMUSCL(
              parameters.identifier_, std::get<0>(parameters.others_),
              std::get<1>(parameters.others_))
    {
    }

    virtual ~EulerianMultiphaseIntegration1stHalfMUSCL() = default;

    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    MultiphaseMUSCLBridge bridge_;
    Vecd *rho_grad_, *p_grad_, *alpha_grad_;
    Matd *vel_grad_;
};

template <>
class EulerianMultiphaseIntegration1stHalfMUSCL<Contact<Wall>>
    : public InteractionWithWall<BaseIntegrationInMultiphaseForWall>
{
  public:
    EulerianMultiphaseIntegration1stHalfMUSCL(BaseContactRelation &contact_relation,
                                        MultiphaseMixture &mixture,
                                        const SecondOrderConfig &cfg);

    template <typename BodyRelationType, typename FirstArg, typename SecondArg>
    explicit EulerianMultiphaseIntegration1stHalfMUSCL(
        DynamicsArgs<BodyRelationType, FirstArg, SecondArg> parameters)
        : EulerianMultiphaseIntegration1stHalfMUSCL(
              parameters.identifier_, std::get<0>(parameters.others_),
              std::get<1>(parameters.others_))
    {
    }

    virtual ~EulerianMultiphaseIntegration1stHalfMUSCL() = default;

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    MultiphaseMUSCLBridge bridge_;
    Vecd *rho_grad_, *p_grad_, *alpha_grad_;
    Matd *vel_grad_;
};

template <typename... InteractionTypes>
class EulerianMultiphaseIntegration2ndHalfMUSCL;

template <>
class EulerianMultiphaseIntegration2ndHalfMUSCL<Inner<>> : public BaseIntegrationInMultiphase
{
  public:
    EulerianMultiphaseIntegration2ndHalfMUSCL(BaseInnerRelation &inner_relation,
                                        MultiphaseMixture &mixture,
                                        const SecondOrderConfig &cfg);

    template <typename BodyRelationType, typename FirstArg, typename SecondArg>
    explicit EulerianMultiphaseIntegration2ndHalfMUSCL(
        DynamicsArgs<BodyRelationType, FirstArg, SecondArg> parameters)
        : EulerianMultiphaseIntegration2ndHalfMUSCL(
              parameters.identifier_, std::get<0>(parameters.others_),
              std::get<1>(parameters.others_))
    {
    }

    virtual ~EulerianMultiphaseIntegration2ndHalfMUSCL() = default;

    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    MultiphaseMUSCLBridge bridge_;
    Vecd *rho_grad_, *p_grad_, *alpha_grad_;
    Matd *vel_grad_;
};

template <>
class EulerianMultiphaseIntegration2ndHalfMUSCL<Contact<Wall>>
    : public InteractionWithWall<BaseIntegrationInMultiphaseForWall>
{
  public:
    EulerianMultiphaseIntegration2ndHalfMUSCL(BaseContactRelation &contact_relation,
                                        MultiphaseMixture &mixture,
                                        const SecondOrderConfig &cfg);

    template <typename BodyRelationType, typename FirstArg, typename SecondArg>
    explicit EulerianMultiphaseIntegration2ndHalfMUSCL(
        DynamicsArgs<BodyRelationType, FirstArg, SecondArg> parameters)
        : EulerianMultiphaseIntegration2ndHalfMUSCL(
              parameters.identifier_, std::get<0>(parameters.others_),
              std::get<1>(parameters.others_))
    {
    }

    virtual ~EulerianMultiphaseIntegration2ndHalfMUSCL() = default;

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    MultiphaseMUSCLBridge bridge_;
    Vecd *rho_grad_, *p_grad_, *alpha_grad_;
    Matd *vel_grad_;
};

//----------------------------------------------------------------------
//	MUSCL with-wall aliases (inner + wall contact).
//----------------------------------------------------------------------
using EulerianMultiphaseIntegration1stHalfMUSCLWithWall =
    ComplexInteraction<EulerianMultiphaseIntegration1stHalfMUSCL<Inner<>, Contact<Wall>>>;

using EulerianMultiphaseIntegration2ndHalfMUSCLWithWall =
    ComplexInteraction<EulerianMultiphaseIntegration2ndHalfMUSCL<Inner<>, Contact<Wall>>>;

/**
 * @class EulerianMultiphaseAcousticTimeStepSize
 * @brief Time step based on the Wood mixture sound speed.
 */
class EulerianMultiphaseAcousticTimeStepSize : public AcousticTimeStep
{
  protected:
    Real *rho_, *p_, *alpha_;
    Vecd *vel_;
    Real smoothing_length_;
    // By value for the same lifetime reason as the integrators.
    MultiphaseMixture mixture_;

  public:
    explicit EulerianMultiphaseAcousticTimeStepSize(SPHBody &sph_body,
                                              MultiphaseMixture &mixture,
                                              Real acousticCFL = 0.6);
    virtual ~EulerianMultiphaseAcousticTimeStepSize() = default;

    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;
};

} // namespace fluid_dynamics
} // namespace SPH

#endif // EULERIAN_MULTIPHASE_INTEGRATION_H
