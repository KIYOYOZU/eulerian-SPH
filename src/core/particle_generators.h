#pragma once

#include "particle_state.h"

namespace eulerian_sph
{
template <class Predicate>
ParticleSet makeCartesianLattice(const Vec2 &lower, const Vec2 &upper, const Scalar spacing, Predicate predicate)
{
    ParticleSet particles;
    const Scalar volume = spacing * spacing;
    const Vec2 extent = upper - lower;
    const int cells_x = std::max(1, static_cast<int>(std::ceil(extent.x() / spacing)));
    const int cells_y = std::max(1, static_cast<int>(std::ceil(extent.y() / spacing)));

    for (int i = 0; i != cells_x; ++i)
    {
        for (int j = 0; j != cells_y; ++j)
        {
            ParticleState particle;
            particle.pos = lower + Vec2((static_cast<Scalar>(i) + 0.5) * spacing,
                                        (static_cast<Scalar>(j) + 0.5) * spacing);
            particle.volume = volume;
            if (predicate(particle.pos))
            {
                particles.addParticle(particle);
            }
        }
    }

    return particles;
}

inline ParticleSet makeCartesianLattice(const Vec2 &lower, const Vec2 &upper, const Scalar spacing)
{
    return makeCartesianLattice(lower, upper, spacing, [](const Vec2 &) { return true; });
}
} // namespace eulerian_sph
