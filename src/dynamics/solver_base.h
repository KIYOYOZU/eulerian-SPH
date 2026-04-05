#pragma once

#include "../core/simulation_context.h"

namespace eulerian_sph
{
namespace dynamics
{
class EulerianSolver
{
  public:
    virtual ~EulerianSolver() = default;

    virtual std::string name() const = 0;

    virtual void initialize(SimulationContext &context)
    {
#if defined(EULERIAN_SPH_USE_OPENMP)
#pragma omp parallel for schedule(static)
#endif
        for (std::int64_t index_i = 0; index_i < static_cast<std::int64_t>(context.particles.size()); ++index_i)
        {
            ParticleState &particle = context.particles[static_cast<std::size_t>(index_i)];
            synchronizeConservativeVariables(particle);
        }
    }

    virtual bool supportsSplitAdvance() const { return false; }

    virtual void advanceFirstHalf(SimulationContext &context, Scalar dt)
    {
        advance(context, dt);
    }

    virtual void advanceSecondHalf(SimulationContext &context, Scalar dt)
    {
        (void)context;
        (void)dt;
    }

    virtual Scalar acceptedTimeStep(const SimulationContext &context, Scalar dt) const
    {
        return context.accepted_dt > 0.0 ? context.accepted_dt : dt;
    }

    virtual Scalar computeTimeStep(const SimulationContext &context) const = 0;
    virtual void advance(SimulationContext &context, Scalar dt)
    {
        (void)context;
        (void)dt;
    }
};
} // namespace dynamics
} // namespace eulerian_sph
