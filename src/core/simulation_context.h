#pragma once

#include "config.h"
#include "grid.h"
#include "wall_particles.h"

#include "../kernel/kernel_base.h"

namespace eulerian_sph
{
struct GhostBoundaryPair
{
    std::size_t real_index = 0;
    std::size_t ghost_local_index = 0;
    Vec2 direction = Vec2::Zero();
};

struct SimulationContext
{
    struct GhostParticleRelation
    {
        std::size_t real_index = 0;
        std::size_t ghost_local_index = 0;
        Vec2 direction = Vec2::UnitX();
        Scalar signed_distance = 0.0;
    };

    SimulationConfig config;
    ParticleSet particles;
    ParticleSet ghost_particles;
    ParticleSet interaction_particles;
    WallParticleSet wall_particles;
    NeighborList neighbors;
    NeighborList wall_neighbors;
    NeighborList wall_inner_neighbors;
    NeighborList wall_contact_neighbors;
    std::vector<Mat2> correction_matrices;
    std::vector<Mat2> wall_correction_matrices;
    std::vector<GhostParticleRelation> ghost_relations;
    UniformGrid grid;
    std::shared_ptr<kernel::KernelBase> kernel;
    std::function<void(SimulationContext &)> mid_step_callback;
    std::array<bool, 2> periodic_axes{false, false};
    Vec2 periodic_lower = Vec2::Zero();
    Vec2 periodic_upper = Vec2::Zero();
    std::vector<GhostBoundaryPair> ghost_boundary_pairs;
    std::function<ParticleState(const ParticleState &)> ghost_projector;
    std::vector<Scalar> position_divergence;
    Scalar current_time = 0.0;
    Scalar accepted_dt = 0.0;
    std::size_t current_step = 0;
};
} // namespace eulerian_sph
