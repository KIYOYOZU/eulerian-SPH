#pragma once

#include "../core/parallel.h"
#include "../core/simulation_context.h"

namespace eulerian_sph
{
namespace dynamics
{
class LinearGradientCorrection
{
  public:
    void rebuild(SimulationContext &context) const
    {
        const std::size_t real_particle_count = context.particles.size();
        context.correction_matrices.assign(real_particle_count, Mat2::Identity());
        context.wall_correction_matrices.assign(context.wall_particles.size(), Mat2::Identity());

        const auto compute_correction_matrix =
            [](const Mat2 &local_configuration)
            {
                const Scalar determinant = local_configuration.determinant();
                const Mat2 regularized = local_configuration + 1.0e-8 * Mat2::Identity();
                const Mat2 inverse = regularized.inverse();
                return (std::abs(determinant) > 1.0e-10 && inverse.allFinite() && inverse.norm() < 50.0)
                           ? inverse
                           : Mat2::Identity();
            };

        const LoopIndex loop_particle_count = static_cast<LoopIndex>(real_particle_count);
#if defined(EULERIAN_SPH_USE_OPENMP)
#pragma omp parallel for schedule(static) if(shouldParallelize(real_particle_count))
#endif
        for (LoopIndex loop_i = 0; loop_i < loop_particle_count; ++loop_i)
        {
            const std::size_t index_i = static_cast<std::size_t>(loop_i);
            Mat2 local_configuration = Mat2::Zero();
            for (const Neighbor &neighbor : context.neighbors[index_i])
            {
                const ParticleState &particle_j = context.interaction_particles[neighbor.index];
                const Vec2 grad_w_ij = neighbor.dweight_dr * particle_j.volume * neighbor.direction;
                const Vec2 r_ji = neighbor.distance * neighbor.direction;
                local_configuration -= r_ji * grad_w_ij.transpose();
            }
            for (const Neighbor &neighbor : context.wall_neighbors[index_i])
            {
                const WallParticle &wall_particle = context.wall_particles[neighbor.index];
                const Vec2 grad_w_ij = neighbor.dweight_dr * wall_particle.volume * neighbor.direction;
                const Vec2 r_ji = neighbor.distance * neighbor.direction;
                local_configuration -= r_ji * grad_w_ij.transpose();
            }

            context.correction_matrices[index_i] = compute_correction_matrix(local_configuration);
        }

        const LoopIndex wall_particle_count = static_cast<LoopIndex>(context.wall_particles.size());
#if defined(EULERIAN_SPH_USE_OPENMP)
#pragma omp parallel for schedule(static) if(shouldParallelize(context.wall_particles.size()))
#endif
        for (LoopIndex loop_i = 0; loop_i < wall_particle_count; ++loop_i)
        {
            const std::size_t index_i = static_cast<std::size_t>(loop_i);
            Mat2 local_configuration = Mat2::Zero();

            for (const Neighbor &neighbor : context.wall_inner_neighbors[index_i])
            {
                const WallParticle &wall_particle = context.wall_particles[neighbor.index];
                const Vec2 grad_w_ij = neighbor.dweight_dr * wall_particle.volume * neighbor.direction;
                const Vec2 r_ji = neighbor.distance * neighbor.direction;
                local_configuration -= r_ji * grad_w_ij.transpose();
            }

            for (const Neighbor &neighbor : context.wall_contact_neighbors[index_i])
            {
                const ParticleState &particle = context.particles[neighbor.index];
                const Vec2 grad_w_ij = neighbor.dweight_dr * particle.volume * neighbor.direction;
                const Vec2 r_ji = neighbor.distance * neighbor.direction;
                local_configuration -= r_ji * grad_w_ij.transpose();
            }

            context.wall_correction_matrices[index_i] = compute_correction_matrix(local_configuration);
        }

        const LoopIndex correction_count = static_cast<LoopIndex>(real_particle_count);
#if defined(EULERIAN_SPH_USE_OPENMP)
#pragma omp parallel for schedule(static) if(shouldParallelize(real_particle_count))
#endif
        for (LoopIndex loop_i = 0; loop_i < correction_count; ++loop_i)
        {
            const std::size_t index_i = static_cast<std::size_t>(loop_i);
            for (Neighbor &neighbor : context.neighbors[index_i])
            {
                const Vec2 displacement = neighbor.distance * neighbor.direction;
                const Mat2 &b_i = context.correction_matrices[index_i];
                const Mat2 b_j = neighbor.index < real_particle_count
                                     ? context.correction_matrices[neighbor.index]
                                     : Mat2::Identity();
                const Vec2 corrected_direction = 0.5 * (b_i + b_j) * neighbor.direction;
                const Scalar corrected_norm = corrected_direction.norm();
                if (corrected_norm <= kEpsilon)
                {
                    continue;
                }
                const Scalar scaling = std::clamp(corrected_norm, Scalar(0.5), Scalar(2.0));

                neighbor.dweight_dr *= scaling;
                neighbor.direction = corrected_direction / corrected_norm;
                neighbor.distance = std::max(std::abs(displacement.dot(neighbor.direction)), kEpsilon);
            }

            if (context.particles[index_i].indicator == 1)
            {
                Vec2 gradient_summation = Vec2::Zero();
                for (const Neighbor &neighbor : context.neighbors[index_i])
                {
                    if (neighbor.index < real_particle_count)
                    {
                        gradient_summation +=
                            neighbor.dweight_dr *
                            context.interaction_particles[neighbor.index].volume *
                            neighbor.direction;
                    }
                }

                const Vec2 ghost_gradient = -gradient_summation;
                const Scalar ghost_gradient_norm = ghost_gradient.norm();
                if (ghost_gradient_norm > kEpsilon)
                {
                    for (Neighbor &neighbor : context.neighbors[index_i])
                    {
                        if (neighbor.index < real_particle_count)
                        {
                            continue;
                        }

                        const Scalar ghost_volume = context.interaction_particles[neighbor.index].volume;
                        const Scalar corrected_dweight = -ghost_gradient_norm * safeInverse(ghost_volume);
                        neighbor.dweight_dr = corrected_dweight;
                        neighbor.direction = ghost_gradient / (ghost_gradient_norm + kEpsilon);
                    }
                }
            }

            for (Neighbor &neighbor : context.wall_neighbors[index_i])
            {
                const Vec2 displacement = neighbor.distance * neighbor.direction;
                const Mat2 &b_i = context.correction_matrices[index_i];
                const Mat2 b_j = neighbor.index < context.wall_correction_matrices.size()
                                     ? context.wall_correction_matrices[neighbor.index]
                                     : Mat2::Identity();
                const Vec2 corrected_direction = 0.5 * (b_i + b_j) * neighbor.direction;
                const Scalar corrected_norm = corrected_direction.norm();
                if (corrected_norm <= kEpsilon)
                {
                    continue;
                }
                const Scalar scaling = std::clamp(corrected_norm, Scalar(0.5), Scalar(2.0));
                neighbor.dweight_dr *= scaling;
                neighbor.direction = corrected_direction / corrected_norm;
                neighbor.distance = std::max(std::abs(displacement.dot(neighbor.direction)), kEpsilon);
            }
        }

        for (const SimulationContext::GhostParticleRelation &relation : context.ghost_relations)
        {
            if (relation.real_index >= context.neighbors.size())
            {
                continue;
            }

            Vec2 gradient_summation = Vec2::Zero();
            for (const Neighbor &neighbor : context.neighbors[relation.real_index])
            {
                if (neighbor.index >= real_particle_count)
                {
                    continue;
                }
                const ParticleState &particle = context.interaction_particles[neighbor.index];
                gradient_summation += neighbor.dweight_dr * particle.volume * neighbor.direction;
            }

            const Scalar gradient_norm = gradient_summation.norm();
            if (gradient_norm <= kEpsilon)
            {
                continue;
            }

            const std::size_t ghost_interaction_index = real_particle_count + relation.ghost_local_index;
            if (ghost_interaction_index >= context.interaction_particles.size())
            {
                continue;
            }

            for (Neighbor &neighbor : context.neighbors[relation.real_index])
            {
                if (neighbor.index != ghost_interaction_index)
                {
                    continue;
                }

                const Scalar ghost_dweight_times_volume = -gradient_norm;
                const Scalar ghost_volume = std::max(context.interaction_particles[ghost_interaction_index].volume, kEpsilon);
                neighbor.dweight_dr = ghost_dweight_times_volume * safeInverse(ghost_volume);
                neighbor.direction = gradient_summation / gradient_norm;
                neighbor.distance = std::max(2.0 * relation.signed_distance, kEpsilon);
                break;
            }
        }
    }
};
} // namespace dynamics
} // namespace eulerian_sph
