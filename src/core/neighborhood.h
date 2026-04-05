#pragma once

#include "types.h"

namespace eulerian_sph
{
struct Neighbor
{
    std::size_t index = 0;
    Vec2 direction = Vec2::Zero();
    Scalar distance = 0.0;
    Scalar weight = 0.0;
    Scalar dweight_dr = 0.0;
};

using NeighborList = std::vector<std::vector<Neighbor>>;
} // namespace eulerian_sph
