#pragma once

#include "boundary_base.h"

#include "../dynamics/surface_indicator.h"

namespace eulerian_sph
{
namespace boundary
{
using NormalProvider = std::function<Vec2(const ParticleState &)>;

class FarFieldBoundary : public BoundaryCondition
{
  public:
    FarFieldBoundary(std::string boundary_name, ParticleSelector selector, NormalProvider normal_provider,
                     const ParticleState &reference_state, Scalar sound_speed,
                     const bool use_surface_indicator = false, const bool include_wall_contact = false)
        : boundary_name_(std::move(boundary_name)),
          selector_(std::move(selector)),
          normal_provider_(std::move(normal_provider)),
          reference_state_(reference_state),
          sound_speed_(sound_speed),
          use_surface_indicator_(use_surface_indicator),
          include_wall_contact_(include_wall_contact) {}

    std::string name() const override { return boundary_name_; }

    void finalizeInitialization(SimulationContext &context) override
    {
        collectBoundaryStatistics(context);
        applyBoundaryReset(context);
    }

    void applyAfterStep(SimulationContext &context) override
    {
        collectBoundaryStatistics(context);
        applyBoundaryReset(context);
    }

  private:
    void collectBoundaryStatistics(SimulationContext &context)
    {
        if (use_surface_indicator_)
        {
            dynamics::updateFreeSurfaceIndicator(context, include_wall_contact_);
            dynamics::updateSmearedSurfaceIndicator(context);
        }

        const std::size_t real_particle_count = context.particles.size();
        active_.assign(real_particle_count, 0);
        inner_weight_summation_.assign(real_particle_count, 0.0);
        rho_average_.assign(real_particle_count, 0.0);
        vel_normal_average_.assign(real_particle_count, 0.0);
        vel_tangential_average_.assign(real_particle_count, Vec2::Zero());
        vel_average_.assign(real_particle_count, Vec2::Zero());
        const Scalar w0 = context.kernel ? context.kernel->weight(0.0) : 1.0;

#if defined(EULERIAN_SPH_USE_OPENMP)
#pragma omp parallel for schedule(static)
#endif
        for (std::int64_t index_i = 0; index_i < static_cast<std::int64_t>(real_particle_count); ++index_i)
        {
            ParticleState &particle = context.particles[static_cast<std::size_t>(index_i)];
            if (!selector_(particle))
            {
                continue;
            }

            Vec2 normal = normal_provider_(particle);
            if (normal.squaredNorm() <= kEpsilon)
            {
                normal = Vec2::UnitX();
            }
            particle.normal = normal.normalized();

            if (use_surface_indicator_ &&
                particle.indicator != 1 &&
                particle.smeared_surface != 1)
            {
                continue;
            }

            active_[static_cast<std::size_t>(index_i)] = 1;
            Scalar inner_weight_summation = w0 * particle.volume;
            Scalar rho_summation = 0.0;
            Scalar vel_normal_summation = 0.0;
            Vec2 vel_tangential_summation = Vec2::Zero();
            Vec2 vel_summation = Vec2::Zero();
            std::size_t total_inner_neighbors = 0;

            for (const Neighbor &neighbor : context.neighbors[static_cast<std::size_t>(index_i)])
            {
                if (neighbor.index >= real_particle_count)
                {
                    continue;
                }

                const ParticleState &neighbor_particle = context.interaction_particles[neighbor.index];
                if (neighbor_particle.indicator == 1)
                {
                    continue;
                }

                inner_weight_summation += neighbor.weight * neighbor_particle.volume;
                rho_summation += neighbor_particle.rho;
                vel_normal_summation += neighbor_particle.vel.dot(particle.normal);
                vel_tangential_summation +=
                    neighbor_particle.vel - neighbor_particle.vel.dot(particle.normal) * particle.normal;
                vel_summation += neighbor_particle.vel;
                total_inner_neighbors += 1;
            }

            inner_weight_summation_[static_cast<std::size_t>(index_i)] = inner_weight_summation;
            if (total_inner_neighbors == 0)
            {
                continue;
            }

            const Scalar inverse_count = 1.0 / static_cast<Scalar>(total_inner_neighbors);
            rho_average_[static_cast<std::size_t>(index_i)] = rho_summation * inverse_count;
            vel_normal_average_[static_cast<std::size_t>(index_i)] = vel_normal_summation * inverse_count;
            vel_tangential_average_[static_cast<std::size_t>(index_i)] = vel_tangential_summation * inverse_count;
            vel_average_[static_cast<std::size_t>(index_i)] = vel_summation * inverse_count;
        }
    }

    void applyBoundaryReset(SimulationContext &context) const
    {
        const std::size_t real_particle_count = context.particles.size();

#if defined(EULERIAN_SPH_USE_OPENMP)
#pragma omp parallel for schedule(static)
#endif
        for (std::int64_t index_i = 0; index_i < static_cast<std::int64_t>(real_particle_count); ++index_i)
        {
            if (active_[static_cast<std::size_t>(index_i)] == 0)
            {
                continue;
            }

            ParticleState &particle = context.particles[static_cast<std::size_t>(index_i)];
            const Vec2 normal = particle.normal;
            const Scalar velocity_farfield_normal = reference_state_.vel.dot(normal);
            const Scalar velocity_boundary_normal = particle.vel.dot(normal);
            const bool treat_as_inflow = normal.x() <= 0.0 || std::abs(normal.y()) > std::abs(normal.x());

            if (treat_as_inflow)
            {
                if (std::abs(velocity_boundary_normal) >= sound_speed_)
                {
                    assignReferenceState(particle);
                }
                else
                {
                    particle.rho = rho_average_[static_cast<std::size_t>(index_i)] *
                                       inner_weight_summation_[static_cast<std::size_t>(index_i)] +
                                   reference_state_.rho *
                                       (1.0 - inner_weight_summation_[static_cast<std::size_t>(index_i)]);
                    particle.p = reference_state_.p +
                                 sound_speed_ * sound_speed_ * (particle.rho - reference_state_.rho);
                    const Scalar vel_normal =
                        vel_normal_average_[static_cast<std::size_t>(index_i)] *
                            inner_weight_summation_[static_cast<std::size_t>(index_i)] +
                        velocity_farfield_normal *
                            (1.0 - inner_weight_summation_[static_cast<std::size_t>(index_i)]);
                    particle.vel = vel_normal * normal +
                                   (reference_state_.vel - velocity_farfield_normal * normal);
                    synchronizeConservativeVariables(particle);
                }
            }
            else
            {
                if (std::abs(velocity_boundary_normal) >= sound_speed_)
                {
                    particle.rho = std::max(rho_average_[static_cast<std::size_t>(index_i)], reference_state_.rho * 0.25);
                    particle.vel = vel_average_[static_cast<std::size_t>(index_i)];
                    synchronizeConservativeVariables(particle);
                    particle.p = reference_state_.p +
                                 sound_speed_ * sound_speed_ * (particle.rho - reference_state_.rho);
                }
                else
                {
                    particle.rho = rho_average_[static_cast<std::size_t>(index_i)] *
                                       inner_weight_summation_[static_cast<std::size_t>(index_i)] +
                                   reference_state_.rho *
                                       (1.0 - inner_weight_summation_[static_cast<std::size_t>(index_i)]);
                    particle.p = reference_state_.p +
                                 sound_speed_ * sound_speed_ * (particle.rho - reference_state_.rho);
                    const Scalar vel_normal =
                        vel_normal_average_[static_cast<std::size_t>(index_i)] *
                            inner_weight_summation_[static_cast<std::size_t>(index_i)] +
                        velocity_farfield_normal *
                            (1.0 - inner_weight_summation_[static_cast<std::size_t>(index_i)]);
                    particle.vel = vel_normal * normal +
                                   vel_tangential_average_[static_cast<std::size_t>(index_i)];
                    synchronizeConservativeVariables(particle);
                }
            }

            particle.energy = reference_state_.energy * particle.volume;
        }
    }

    void assignReferenceState(ParticleState &particle) const
    {
        const Vec2 original_position = particle.pos;
        const Scalar original_volume = particle.volume;
        particle = reference_state_;
        particle.pos = original_position;
        particle.volume = original_volume;
        synchronizeConservativeVariables(particle);
        particle.energy = reference_state_.energy * original_volume;
    }

    std::string boundary_name_;
    ParticleSelector selector_;
    NormalProvider normal_provider_;
    ParticleState reference_state_;
    Scalar sound_speed_;
    bool use_surface_indicator_;
    bool include_wall_contact_;
    std::vector<int> active_;
    std::vector<Scalar> inner_weight_summation_;
    std::vector<Scalar> rho_average_;
    std::vector<Scalar> vel_normal_average_;
    std::vector<Vec2> vel_tangential_average_;
    std::vector<Vec2> vel_average_;
};
} // namespace boundary
} // namespace eulerian_sph
