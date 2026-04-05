#pragma once

#include "boundary_base.h"

#include "../dynamics/surface_indicator.h"

namespace eulerian_sph
{
namespace boundary
{
using BoundaryNormalProvider = std::function<Vec2(const ParticleState &)>;
using SignedDistanceProvider = std::function<Scalar(const ParticleState &)>;
using GhostStateUpdater =
    std::function<void(ParticleState &, const SimulationContext &, std::size_t, const Vec2 &)>;

class GhostBoundary : public BoundaryCondition
{
  public:
    GhostBoundary(std::string boundary_name, ParticleSelector selector, BoundaryNormalProvider normal_provider,
                  SignedDistanceProvider signed_distance_provider, GhostStateUpdater state_updater,
                  BoundaryTypeProvider boundary_type_provider = {},
                  const bool use_surface_indicator = true)
        : boundary_name_(std::move(boundary_name)),
          selector_(std::move(selector)),
          normal_provider_(std::move(normal_provider)),
          signed_distance_provider_(std::move(signed_distance_provider)),
          state_updater_(std::move(state_updater)),
          boundary_type_provider_(std::move(boundary_type_provider)),
          use_surface_indicator_(use_surface_indicator) {}

    std::string name() const override { return boundary_name_; }

    void applyBeforeStep(SimulationContext &context) override
    {
        dynamics::updateFreeSurfaceIndicator(context, false);

        real_indices_.clear();
        ghost_normals_.clear();
        context.ghost_relations.clear();
        context.ghost_particles.reserve(context.particles.size());
        for (ParticleState &particle : context.particles.data())
        {
            particle.boundary_type = 0;
        }

        for (std::size_t index_i = 0; index_i != context.particles.size(); ++index_i)
        {
            ParticleState &particle = context.particles[index_i];
            if (!selector_(particle))
            {
                continue;
            }
            if (use_surface_indicator_ && particle.indicator != 1)
            {
                continue;
            }

            particle.boundary_type =
                boundary_type_provider_ ? boundary_type_provider_(particle) : particle.boundary_type;

            Vec2 shape_normal = normal_provider_(particle);
            if (shape_normal.squaredNorm() <= kEpsilon)
            {
                continue;
            }
            shape_normal.normalize();

            Vec2 gradient_summation = Vec2::Zero();
            for (const Neighbor &neighbor : context.neighbors[index_i])
            {
                if (neighbor.index >= context.particles.size())
                {
                    continue;
                }

                const ParticleState &neighbor_particle = context.particles[neighbor.index];
                gradient_summation += neighbor.dweight_dr * neighbor_particle.volume * neighbor.direction;
            }

            Vec2 ghost_normal = shape_normal;
            if (gradient_summation.squaredNorm() > kEpsilon)
            {
                ghost_normal = gradient_summation.normalized();
                if (ghost_normal.dot(shape_normal) < 0.0)
                {
                    ghost_normal = -ghost_normal;
                }
            }

            particle.normal = ghost_normal;
            particle.signed_distance = signed_distance_provider_(particle);

            ParticleState ghost_particle = particle;
            const Scalar distance = std::max(std::abs(particle.signed_distance), context.config.spatial.particle_spacing);
            ghost_particle.pos = particle.pos + distance * ghost_normal;
            ghost_particle.normal = ghost_normal;
            ghost_particle.signed_distance = distance;
            ghost_particle.boundary_type = particle.boundary_type;
            state_updater_(ghost_particle, context, index_i, ghost_normal);
            synchronizeConservativeVariables(ghost_particle);

            SimulationContext::GhostParticleRelation relation;
            relation.real_index = index_i;
            relation.ghost_local_index = context.ghost_particles.size();
            relation.direction = ghost_normal;
            relation.signed_distance = distance;
            context.ghost_relations.push_back(relation);
            context.ghost_particles.addParticle(ghost_particle);
            real_indices_.push_back(index_i);
            ghost_normals_.push_back(ghost_normal);
        }
    }

    void applyMidStep(SimulationContext &context) override
    {
        const std::size_t real_particle_count = context.particles.size();
        context.ghost_relations.clear();
        for (std::size_t local_index = 0; local_index != real_indices_.size(); ++local_index)
        {
            ParticleState &ghost_particle = context.ghost_particles[local_index];
            state_updater_(ghost_particle, context, real_indices_[local_index], ghost_normals_[local_index]);
            synchronizeConservativeVariables(ghost_particle);

            SimulationContext::GhostParticleRelation relation;
            relation.real_index = real_indices_[local_index];
            relation.ghost_local_index = local_index;
            relation.direction = ghost_normals_[local_index];
            relation.signed_distance = ghost_particle.signed_distance;
            context.ghost_relations.push_back(relation);
            if (real_particle_count + local_index < context.interaction_particles.size())
            {
                context.interaction_particles[real_particle_count + local_index] = ghost_particle;
            }
        }
    }

  private:
    std::string boundary_name_;
    ParticleSelector selector_;
    BoundaryNormalProvider normal_provider_;
    SignedDistanceProvider signed_distance_provider_;
    GhostStateUpdater state_updater_;
    BoundaryTypeProvider boundary_type_provider_;
    bool use_surface_indicator_;
    std::vector<std::size_t> real_indices_;
    std::vector<Vec2> ghost_normals_;
};
} // namespace boundary
} // namespace eulerian_sph
