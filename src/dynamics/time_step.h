#pragma once

#include "../core/parallel.h"
#include "../core/simulation_context.h"

#include <algorithm>

namespace eulerian_sph
{
namespace dynamics
{
class AcousticTimeStepEstimator
{
  public:
    Scalar estimate(const SimulationContext &context, const Scalar sound_speed, const Scalar acoustic_cfl) const
    {
        Scalar max_speed = 0.0;
#if defined(EULERIAN_SPH_USE_OPENMP)
#pragma omp parallel if(shouldParallelize(context.particles.size()))
        {
            Scalar local_max_speed = 0.0;
#pragma omp for nowait
            for (LoopIndex index_i = 0; index_i < static_cast<LoopIndex>(context.particles.size()); ++index_i)
            {
                local_max_speed = std::max(
                    local_max_speed,
                    context.particles[static_cast<std::size_t>(index_i)].vel.norm());
            }
#pragma omp critical
            { max_speed = std::max(max_speed, local_max_speed); }
        }
#else
        for (const ParticleState &particle : context.particles.data())
        {
            max_speed = std::max(max_speed, particle.vel.norm());
        }
#endif

        const Scalar h = context.kernel ? context.kernel->smoothingLength()
                                        : context.config.spatial.particle_spacing * context.config.spatial.smoothing_length_ratio;
        return acoustic_cfl * h * safeInverse(sound_speed + max_speed + kEpsilon);
    }
};
} // namespace dynamics
} // namespace eulerian_sph
