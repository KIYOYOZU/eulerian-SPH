#pragma once

#include "../core/parallel.h"
#include "../core/simulation_context.h"

namespace eulerian_sph
{
namespace dynamics
{
using PositionNormalProvider = std::function<Vec2(const Vec2 &)>;
using SignedDistanceProvider = std::function<Scalar(const Vec2 &)>;

inline Vec2 normalizedOr(const Vec2 &candidate, const Vec2 &fallback = Vec2::UnitX())
{
    return candidate.norm() <= kEpsilon ? fallback : candidate.normalized();
}

inline void updateParticleShapeState(ParticleSet &particles, const PositionNormalProvider &normal_provider,
                                     const SignedDistanceProvider &signed_distance_provider)
{
    parallel::forRange(particles.size(), [&](const std::size_t index)
    {
        ParticleState &particle = particles[index];
        particle.normal = normalizedOr(normal_provider(particle.pos));
        particle.signed_distance = signed_distance_provider(particle.pos);
    });
}

inline void updateWallShapeState(WallParticleSet &wall_particles, const PositionNormalProvider &normal_provider,
                                 const SignedDistanceProvider &signed_distance_provider)
{
    parallel::forRange(wall_particles.size(), [&](const std::size_t index)
    {
        WallParticle &particle = wall_particles[index];
        particle.normal = normalizedOr(normal_provider(particle.pos));
        particle.signed_distance = signed_distance_provider(particle.pos);
    });
}

inline void freeSurfaceIndicationInner(SimulationContext &context)
{
    const std::size_t real_particle_count = context.particles.size();
    const Scalar threshold = 1.5;
    const Scalar smoothing_length =
        context.kernel ? context.kernel->smoothingLength()
                       : context.config.spatial.particle_spacing * context.config.spatial.smoothing_length_ratio;
    std::vector<Scalar> position_divergence(real_particle_count, 0.0);

    parallel::forRange(real_particle_count, [&](const std::size_t index_i)
    {
        Scalar pos_divergence = 0.0;
        for (const Neighbor &neighbor : context.neighbors[index_i])
        {
            if (neighbor.index >= real_particle_count)
            {
                continue;
            }
            const ParticleState &particle_j = context.interaction_particles[neighbor.index];
            pos_divergence -= neighbor.dweight_dr * particle_j.volume * neighbor.distance;
        }
        position_divergence[index_i] = pos_divergence;
    });

    parallel::forRange(real_particle_count, [&](const std::size_t index_i)
    {
        ParticleState &particle = context.particles[index_i];
        particle.position_divergence = position_divergence[index_i];
        particle.indicator = 1;
        if (position_divergence[index_i] <= threshold)
        {
            return;
        }

        bool is_near_surface = false;
        for (const Neighbor &neighbor : context.neighbors[index_i])
        {
            if (neighbor.index >= real_particle_count)
            {
                continue;
            }
            if (position_divergence[neighbor.index] < threshold && neighbor.distance < smoothing_length)
            {
                is_near_surface = true;
                break;
            }
        }
        if (!is_near_surface)
        {
            particle.indicator = 0;
        }
    });
}

inline void freeSurfaceIndicationComplex(SimulationContext &context)
{
    const std::size_t real_particle_count = context.particles.size();
    const Scalar threshold = 1.5;
    const Scalar smoothing_length =
        context.kernel ? context.kernel->smoothingLength()
                       : context.config.spatial.particle_spacing * context.config.spatial.smoothing_length_ratio;
    std::vector<Scalar> position_divergence(real_particle_count, 0.0);

    parallel::forRange(real_particle_count, [&](const std::size_t index_i)
    {
        Scalar pos_divergence = 0.0;
        for (const Neighbor &neighbor : context.neighbors[index_i])
        {
            if (neighbor.index >= real_particle_count)
            {
                continue;
            }
            const ParticleState &particle_j = context.interaction_particles[neighbor.index];
            pos_divergence -= neighbor.dweight_dr * particle_j.volume * neighbor.distance;
        }
        for (const Neighbor &neighbor : context.wall_neighbors[index_i])
        {
            const WallParticle &wall_particle = context.wall_particles[neighbor.index];
            pos_divergence -= neighbor.dweight_dr * wall_particle.volume * neighbor.distance;
        }
        position_divergence[index_i] = pos_divergence;
    });

    parallel::forRange(real_particle_count, [&](const std::size_t index_i)
    {
        ParticleState &particle = context.particles[index_i];
        particle.position_divergence = position_divergence[index_i];
        particle.indicator = 1;
        if (position_divergence[index_i] <= threshold)
        {
            return;
        }

        bool is_near_surface = false;
        for (const Neighbor &neighbor : context.neighbors[index_i])
        {
            if (neighbor.index >= real_particle_count)
            {
                continue;
            }
            if (position_divergence[neighbor.index] < threshold && neighbor.distance < smoothing_length)
            {
                is_near_surface = true;
                break;
            }
        }
        if (!is_near_surface)
        {
            particle.indicator = 0;
        }
    });
}

inline void smearedSurfaceIndication(SimulationContext &context)
{
    parallel::forRange(context.particles.size(), [&](const std::size_t index_i)
    {
        ParticleState &particle = context.particles[index_i];
        particle.smeared_surface = 0;
        for (const Neighbor &neighbor : context.neighbors[index_i])
        {
            if (neighbor.index >= context.particles.size())
            {
                continue;
            }
            if (context.particles[neighbor.index].indicator == 1)
            {
                particle.smeared_surface = 1;
                break;
            }
        }
    });
}
} // namespace dynamics
} // namespace eulerian_sph
