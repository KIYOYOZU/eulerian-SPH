#pragma once

#include "../core/parallel.h"

#include "kernel_correction.h"
#include "solver_base.h"
#include "time_step.h"

#include "../math/eos.h"
#include "../math/riemann.h"

#include <vector>

namespace eulerian_sph
{
namespace dynamics
{
class WeaklyCompressibleEulerianSolver : public EulerianSolver
{
  public:
    explicit WeaklyCompressibleEulerianSolver(const MaterialConfig &material)
        : eos_(material.reference_density, material.sound_speed, material.reference_pressure),
          reference_density_(material.reference_density),
          sound_speed_(material.sound_speed),
          dynamic_viscosity_(material.dynamic_viscosity) {}

    std::string name() const override { return "WeaklyCompressibleEulerianSolver"; }
    bool supportsSplitAdvance() const override { return true; }

    void initialize(SimulationContext &context) override
    {
        EulerianSolver::initialize(context);
        correction_.rebuild(context);
        parallel::forRange(context.particles.size(),
                           [&](const std::size_t index_i)
                           {
                               ParticleState &particle = context.particles[index_i];
                               particle.rho = std::max(particle.rho, 0.25 * reference_density_);
                               particle.mass = particle.rho * particle.volume;
                               particle.mom = particle.mass * particle.vel;
                               particle.p = eos_.pressureFromDensity(particle.rho);
                           });
    }

    Scalar computeTimeStep(const SimulationContext &context) const override
    {
        return estimator_.estimate(context, eos_.soundSpeed(), context.config.time.acoustic_cfl);
    }

    void advanceFirstHalf(SimulationContext &context, Scalar dt) override
    {
        correction_.rebuild(context);
        math::AcousticRiemannSolver riemann(reference_density_, sound_speed_, 15.0);
        const Scalar smoothing_length =
            context.kernel ? context.kernel->smoothingLength()
                           : context.config.spatial.particle_spacing * context.config.spatial.smoothing_length_ratio;
        const std::size_t real_particle_count = context.particles.size();
        momentum_change_rates_.assign(real_particle_count, Vec2::Zero());
        force_prior_.assign(real_particle_count, Vec2::Zero());

        const LoopIndex loop_particle_count = static_cast<LoopIndex>(real_particle_count);
#if defined(EULERIAN_SPH_USE_OPENMP)
#pragma omp parallel for schedule(static) if(shouldParallelize(real_particle_count))
#endif
        for (LoopIndex loop_i = 0; loop_i < loop_particle_count; ++loop_i)
        {
            const std::size_t index_i = static_cast<std::size_t>(loop_i);
            const ParticleState &particle_i = context.particles[index_i];
            const math::AcousticState state_i{particle_i.rho, particle_i.vel, particle_i.p};
            Vec2 momentum_change_rate = Vec2::Zero();
            Vec2 viscous_force_inner = Vec2::Zero();
            Vec2 viscous_force_wall = Vec2::Zero();
            const Mat2 &b_i = context.correction_matrices[index_i];

            for (const Neighbor &neighbor : context.neighbors[index_i])
            {
                const ParticleState &particle_j = context.interaction_particles[neighbor.index];
                const math::AcousticState state_j{particle_j.rho, particle_j.vel, particle_j.p};
                const math::AcousticStarState interface_state =
                    riemann.interfaceState(state_i, state_j, neighbor.direction);
                const Mat2 convect_flux =
                    interface_state.rho * interface_state.vel * interface_state.vel.transpose();
                const Scalar d_w_ij_v_j = neighbor.dweight_dr * particle_j.volume;

                momentum_change_rate -=
                    2.0 * particle_i.volume *
                    (convect_flux + interface_state.p * Mat2::Identity()) *
                    neighbor.direction * d_w_ij_v_j;

                const Mat2 b_j = neighbor.index < real_particle_count
                                     ? context.correction_matrices[neighbor.index]
                                     : Mat2::Identity();
                const Scalar correction_factor =
                    neighbor.direction.dot((b_i + b_j) * neighbor.direction);
                const Vec2 velocity_derivative =
                    (particle_i.vel - particle_j.vel) /
                    (neighbor.distance + 0.01 * smoothing_length);
                viscous_force_inner +=
                    correction_factor * dynamic_viscosity_ * velocity_derivative *
                    neighbor.dweight_dr * particle_j.volume;
            }

            for (const Neighbor &neighbor : context.wall_neighbors[index_i])
            {
                const WallParticle &wall_particle = context.wall_particles[neighbor.index];
                const math::AcousticState state_j{
                    particle_i.rho,
                    2.0 * wall_particle.velocity_average - particle_i.vel,
                    particle_i.p};
                const math::AcousticStarState interface_state =
                    riemann.interfaceState(state_i, state_j, wall_particle.normal);
                const Mat2 convect_flux =
                    interface_state.rho * interface_state.vel * interface_state.vel.transpose();
                const Scalar d_w_ij_v_j = neighbor.dweight_dr * wall_particle.volume;

                momentum_change_rate -=
                    2.0 * particle_i.volume *
                    (convect_flux + interface_state.p * Mat2::Identity()) *
                    neighbor.direction * d_w_ij_v_j;

                const Scalar correction_factor =
                    2.0 * neighbor.direction.dot(b_i * neighbor.direction);
                const Vec2 velocity_derivative =
                    2.0 * (particle_i.vel - wall_particle.velocity_average) /
                    (neighbor.distance + 0.01 * smoothing_length);
                viscous_force_wall +=
                    correction_factor * dynamic_viscosity_ * velocity_derivative *
                    neighbor.dweight_dr * wall_particle.volume;
            }

            force_prior_[index_i] = (viscous_force_inner + viscous_force_wall) * particle_i.volume;
            momentum_change_rates_[index_i] = momentum_change_rate;
        }

        context.accepted_dt = dt;
        parallel::forRange(real_particle_count,
                           [&](const std::size_t index_i)
                           {
                               ParticleState &particle = context.particles[index_i];
                               particle.mom += (momentum_change_rates_[index_i] + force_prior_[index_i]) * dt;
                               particle.vel = particle.mom * safeInverse(std::max(particle.mass, kEpsilon));
                               if (!particle.vel.allFinite())
                               {
                                   particle.vel.setZero();
                                   particle.mom = particle.mass * particle.vel;
                               }
                           });
    }

    void advanceSecondHalf(SimulationContext &context, Scalar dt) override
    {
        math::AcousticRiemannSolver riemann(reference_density_, sound_speed_, 15.0);
        const std::size_t real_particle_count = context.particles.size();
        mass_change_rates_.assign(real_particle_count, 0.0);

        const LoopIndex loop_particle_count = static_cast<LoopIndex>(real_particle_count);
#if defined(EULERIAN_SPH_USE_OPENMP)
#pragma omp parallel for schedule(static) if(shouldParallelize(real_particle_count))
#endif
        for (LoopIndex loop_i = 0; loop_i < loop_particle_count; ++loop_i)
        {
            const std::size_t index_i = static_cast<std::size_t>(loop_i);
            const ParticleState &particle_i = context.particles[index_i];
            const math::AcousticState state_i{particle_i.rho, particle_i.vel, particle_i.p};
            Scalar mass_change_rate = 0.0;

            for (const Neighbor &neighbor : context.neighbors[index_i])
            {
                const ParticleState &particle_j = context.interaction_particles[neighbor.index];
                const math::AcousticState state_j{particle_j.rho, particle_j.vel, particle_j.p};
                const math::AcousticStarState interface_state =
                    riemann.interfaceState(state_i, state_j, neighbor.direction);
                const Scalar d_w_ij_v_j = neighbor.dweight_dr * particle_j.volume;
                mass_change_rate -=
                    2.0 * particle_i.volume *
                    (interface_state.rho * interface_state.vel).dot(neighbor.direction) *
                    d_w_ij_v_j;
            }

            for (const Neighbor &neighbor : context.wall_neighbors[index_i])
            {
                const WallParticle &wall_particle = context.wall_particles[neighbor.index];
                const math::AcousticState state_j{
                    particle_i.rho,
                    2.0 * wall_particle.velocity_average - particle_i.vel,
                    particle_i.p};
                const math::AcousticStarState interface_state =
                    riemann.interfaceState(state_i, state_j, wall_particle.normal);
                const Scalar d_w_ij_v_j = neighbor.dweight_dr * wall_particle.volume;
                mass_change_rate -=
                    2.0 * particle_i.volume *
                    (interface_state.rho * interface_state.vel).dot(neighbor.direction) *
                    d_w_ij_v_j;
            }

            mass_change_rates_[index_i] = mass_change_rate;
        }

        context.accepted_dt = dt;
        parallel::forRange(real_particle_count,
                           [&](const std::size_t index_i)
                           {
                               ParticleState &particle = context.particles[index_i];
                               const Scalar mass_floor = 0.25 * reference_density_ * particle.volume;
                               particle.mass = std::max(mass_floor, particle.mass + mass_change_rates_[index_i] * dt);
                               particle.rho = particle.mass * safeInverse(particle.volume);
                               particle.vel = particle.mom * safeInverse(std::max(particle.mass, kEpsilon));
                               if (!particle.vel.allFinite())
                               {
                                   particle.vel.setZero();
                                   particle.mom = particle.mass * particle.vel;
                               }
                               particle.p = eos_.pressureFromDensity(particle.rho);
                           });
    }

  private:
    math::WeaklyCompressibleEquationOfState eos_;
    Scalar reference_density_;
    Scalar sound_speed_;
    Scalar dynamic_viscosity_;
    AcousticTimeStepEstimator estimator_;
    LinearGradientCorrection correction_;
    std::vector<Vec2> momentum_change_rates_;
    std::vector<Vec2> force_prior_;
    std::vector<Scalar> mass_change_rates_;
};
} // namespace dynamics
} // namespace eulerian_sph
