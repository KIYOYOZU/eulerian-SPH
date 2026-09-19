/**
 * @file 	time_step_local_h.cpp
 * @author 	KIYOYOZU
 */

#include "time_step_local_h.h"

#include "adaptation.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
EulerianMultiphaseAcousticTimeStepSizeLocalH::EulerianMultiphaseAcousticTimeStepSizeLocalH(
    SPHBody &sph_body, MultiphaseMixture &mixture, Real acousticCFL)
    : AcousticTimeStep(sph_body),
      mixture_(mixture),
      rho_local_(particles_->getVariableDataByName<Real>("Density")),
      p_local_(particles_->getVariableDataByName<Real>("Pressure")),
      alpha_(particles_->getVariableDataByName<Real>("VolumeFraction")),
      h_ratio_(particles_->getVariableDataByName<Real>("SmoothingLengthRatio")),
      vel_local_(particles_->getVariableDataByName<Vecd>("Velocity")),
      h_ref_(sph_body.getSPHAdaptation().ReferenceSmoothingLength())
{
    acousticCFL_ = acousticCFL;
}
//=================================================================================================//
Real EulerianMultiphaseAcousticTimeStepSizeLocalH::reduce(size_t index_i, Real dt)
{
    Real h_local = h_ref_ / h_ratio_[index_i];
    return (mixture_.MixtureSoundSpeed(alpha_[index_i], rho_local_[index_i], p_local_[index_i])
            + vel_local_[index_i].norm())
           / h_local;
}
//=================================================================================================//
Real EulerianMultiphaseAcousticTimeStepSizeLocalH::outputResult(Real reduced_value)
{
    return acousticCFL_ / Dimensions / (reduced_value + TinyReal);
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
