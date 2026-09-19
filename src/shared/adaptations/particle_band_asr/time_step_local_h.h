/**
 * @file 	time_step_local_h.h
 * @brief 	Acoustic time step size using the local smoothing length,
 * 			required by the particle-band ASR method (the paper leaves the
 * 			time stepping unspecified; the framework's Eulerian multiphase
 * 			version uses the global reference smoothing length, which is
 * 			unstable for refined particles with h_local << h_ref).
 * @author 	KIYOYOZU
 */

#ifndef TIME_STEP_LOCAL_H_H
#define TIME_STEP_LOCAL_H_H

#include "fluid_time_step.h"
#include "multiphase_mixture.h"

namespace SPH
{
namespace fluid_dynamics
{
/**
 * @class EulerianMultiphaseAcousticTimeStepSizeLocalH
 * @brief dt = (CFL/d) / max_i [(c_i + |u_i|) / h_local_i],
 * 			with h_local_i = h_ref / h_ratio_i and the mixture (Wood) sound speed.
 */
class EulerianMultiphaseAcousticTimeStepSizeLocalH : public AcousticTimeStep
{
  public:
    EulerianMultiphaseAcousticTimeStepSizeLocalH(
        SPHBody &sph_body, MultiphaseMixture &mixture, Real acousticCFL);
    virtual ~EulerianMultiphaseAcousticTimeStepSizeLocalH() {};

    // Non-virtual by the same name-hiding pattern as the framework's
    // EulerianMultiphaseAcousticTimeStepSize (resolved via ReduceDynamics<T>).
    Real reduce(size_t index_i, Real dt = 0.0);
    Real outputResult(Real reduced_value) override;

  protected:
    MultiphaseMixture &mixture_;
    Real *rho_local_, *p_local_, *alpha_, *h_ratio_;
    Vecd *vel_local_;
    Real h_ref_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // TIME_STEP_LOCAL_H_H
