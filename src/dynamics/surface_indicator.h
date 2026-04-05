#pragma once

#include "../core/simulation_context.h"

#include <cmath>
#include <cstdint>

namespace eulerian_sph
{
namespace dynamics
{
inline void updateFreeSurfaceIndicator(SimulationContext &context, const bool include_wall_contact)
{
    const std::size_t real_particle_count = context.particles.size();
    context.position_divergence.assign(real_particle_count, 0.0);

    const Scalar smoothing_length =
        context.kernel ? context.kernel->smoothingLength()
                       : context.config.spatial.particle_spacing * context.config.spatial.smoothing_length_ratio;
    constexpr Scalar threshold = 0.75 * 2.0;

#if defined(EULERIAN_SPH_USE_OPENMP)
#pragma omp parallel for schedule(static)
#endif
    for (std::int64_t index_i = 0; index_i < static_cast<std::int64_t>(real_particle_count); ++index_i)
    {
        Scalar position_divergence = 0.0;

        for (const Neighbor &neighbor : context.neighbors[static_cast<std::size_t>(index_i)])
        {
            const ParticleState &particle_j = context.interaction_particles[neighbor.index];
            position_divergence -= neighbor.dweight_dr * particle_j.volume * neighbor.distance;
        }

        if (include_wall_contact)
        {
            for (const Neighbor &neighbor : context.wall_neighbors[static_cast<std::size_t>(index_i)])
            {
                const WallParticle &wall_particle = context.wall_particles[neighbor.index];
                position_divergence -= neighbor.dweight_dr * wall_particle.volume * neighbor.distance;
            }
        }

        context.position_divergence[static_cast<std::size_t>(index_i)] = position_divergence;
    }

#if defined(EULERIAN_SPH_USE_OPENMP)
#pragma omp parallel for schedule(static)
#endif
    for (std::int64_t index_i = 0; index_i < static_cast<std::int64_t>(real_particle_count); ++index_i)
    {
        ParticleState &particle = context.particles[static_cast<std::size_t>(index_i)];
        particle.position_divergence = context.position_divergence[static_cast<std::size_t>(index_i)];
        particle.indicator = 1;

        if (context.position_divergence[static_cast<std::size_t>(index_i)] > threshold)
        {
            bool is_near_surface = false;
            for (const Neighbor &neighbor : context.neighbors[static_cast<std::size_t>(index_i)])
            {
                if (context.position_divergence[neighbor.index] < threshold &&
                    neighbor.distance < smoothing_length)
                {
                    is_near_surface = true;
                    break;
                }
            }

            if (!is_near_surface)
            {
                particle.indicator = 0;
            }
        }
    }
}

inline void updateSmearedSurfaceIndicator(SimulationContext &context)
{
    const std::size_t real_particle_count = context.particles.size();

#if defined(EULERIAN_SPH_USE_OPENMP)
#pragma omp parallel for schedule(static)
#endif
    for (std::int64_t index_i = 0; index_i < static_cast<std::int64_t>(real_particle_count); ++index_i)
    {
        ParticleState &particle = context.particles[static_cast<std::size_t>(index_i)];
        particle.smeared_surface = 0;

        for (const Neighbor &neighbor : context.neighbors[static_cast<std::size_t>(index_i)])
        {
            if (neighbor.index < real_particle_count &&
                context.particles[neighbor.index].indicator == 1)
            {
                particle.smeared_surface = 1;
                break;
            }
        }
    }
}
} // namespace dynamics
} // namespace eulerian_sph
