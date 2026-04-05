#pragma once

#include "neighborhood.h"
#include "parallel.h"
#include "particle_state.h"

#include <cmath>
#include <unordered_map>

namespace eulerian_sph
{
class UniformGrid
{
  public:
    explicit UniformGrid(const Scalar cell_size = 0.1) : cell_size_(cell_size) {}

    void setCellSize(const Scalar cell_size)
    {
        cell_size_ = std::max(cell_size, kEpsilon);
    }

    Scalar cellSize() const { return cell_size_; }

    void rebuild(const ParticleSet &particles)
    {
        buckets_.clear();
        for (std::size_t i = 0; i != particles.size(); ++i)
        {
            buckets_[cellOf(particles[i].pos)].push_back(i);
        }
    }

    void queryNeighbors(const ParticleSet &particles, const Scalar search_radius, NeighborList &neighbors,
                        const std::size_t target_count = std::numeric_limits<std::size_t>::max()) const
    {
        const std::size_t active_count = std::min(target_count, particles.size());
        neighbors.assign(active_count, {});
        const Scalar radius_squared = search_radius * search_radius;

        const LoopIndex loop_active_count = static_cast<LoopIndex>(active_count);
#if defined(EULERIAN_SPH_USE_OPENMP)
#pragma omp parallel for schedule(static) if(shouldParallelize(active_count))
#endif
        for (LoopIndex loop_i = 0; loop_i < loop_active_count; ++loop_i)
        {
            const std::size_t i = static_cast<std::size_t>(loop_i);
            const CellCoord cell = cellOf(particles[i].pos);
            auto &particle_neighbors = neighbors[i];

            for (int dx = -1; dx <= 1; ++dx)
            {
                for (int dy = -1; dy <= 1; ++dy)
                {
                    const CellCoord candidate{cell.x + dx, cell.y + dy};
                    const auto found = buckets_.find(candidate);
                    if (found == buckets_.end())
                    {
                        continue;
                    }

                    for (const std::size_t j : found->second)
                    {
                        if (i == j)
                        {
                            continue;
                        }

                        const Vec2 displacement = particles[j].pos - particles[i].pos;
                        const Scalar distance_squared = displacement.squaredNorm();
                        if (distance_squared > radius_squared || distance_squared <= kEpsilon)
                        {
                            continue;
                        }

                        const Scalar distance = std::sqrt(distance_squared);
                        Neighbor neighbor;
                        neighbor.index = j;
                        neighbor.distance = distance;
                        neighbor.direction = displacement / distance;
                        particle_neighbors.push_back(neighbor);
                    }
                }
            }
        }
    }

  private:
    struct CellCoord
    {
        int x = 0;
        int y = 0;

        bool operator==(const CellCoord &other) const
        {
            return x == other.x && y == other.y;
        }
    };

    struct CellCoordHash
    {
        std::size_t operator()(const CellCoord &coord) const
        {
            const auto hx = static_cast<std::size_t>(std::hash<int>{}(coord.x));
            const auto hy = static_cast<std::size_t>(std::hash<int>{}(coord.y));
            return hx ^ (hy << 1);
        }
    };

    CellCoord cellOf(const Vec2 &position) const
    {
        return CellCoord{
            static_cast<int>(std::floor(position.x() / cell_size_)),
            static_cast<int>(std::floor(position.y() / cell_size_))};
    }

    Scalar cell_size_;
    std::unordered_map<CellCoord, std::vector<std::size_t>, CellCoordHash> buckets_;
};
} // namespace eulerian_sph
