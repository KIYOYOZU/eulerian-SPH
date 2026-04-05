#pragma once

#include "../core/parallel.h"

#include "kernel_correction.h"
#include "solver_base.h"
#include "time_step.h"

#include "../math/eos.h"
#include "../math/riemann.h"

#include <algorithm>
#include <vector>

namespace eulerian_sph
{
namespace dynamics
{
class CompressibleEulerianSolver : public EulerianSolver
{
  public:
    explicit CompressibleEulerianSolver(const MaterialConfig &material)
        : eos_(material.heat_capacity_ratio),
          riemann_(eos_, 5.0),
          reference_density_(material.reference_density),
          pressure_floor_(std::max(material.reference_pressure, 1.0e-6)),
          dynamic_viscosity_(material.dynamic_viscosity) {}

    std::string name() const override { return "CompressibleEulerianSolver"; }

    void initialize(SimulationContext &context) override
    {
        EulerianSolver::initialize(context);
        correction_.rebuild(context);
        for (ParticleState &particle : context.particles.data())
        {
            particle.rho = std::max(particle.rho, 0.25 * reference_density_);
            particle.mass = particle.rho * particle.volume;
            particle.mom = particle.mass * particle.vel;
            updatePrimitiveState(particle);
        }
    }

    bool supportsSplitAdvance() const override { return true; }

    Scalar computeTimeStep(const SimulationContext &context) const override
    {
        const Scalar max_sound_speed = parallel::reduceMax(
            context.particles.size(),
            context.config.material.sound_speed,
            [&](const std::size_t index_i)
            {
                const ParticleState &particle = context.particles[index_i];
                return eos_.soundSpeed(std::max(0.0, particle.p), std::max(particle.rho, kEpsilon));
            });
        return 0.5 * estimator_.estimate(context, max_sound_speed, context.config.time.acoustic_cfl);
    }

    void advance(SimulationContext &context, Scalar dt) override
    {
        advanceFirstHalf(context, dt);
        advanceSecondHalf(context, dt);
    }

    void advanceFirstHalf(SimulationContext &context, Scalar dt) override
    {
        correction_.rebuild(context);
        updateGhostKernelGradient(context);

        const std::size_t real_particle_count = context.particles.size();
        force_rates_cache_.assign(real_particle_count, Vec2::Zero());
        force_prior_cache_.assign(real_particle_count, Vec2::Zero());
        std::vector<Scalar> mass_change_rates(real_particle_count, 0.0);
        std::vector<Scalar> energy_change_rates(real_particle_count, 0.0);

        computeMomentumChangeRates(context, force_rates_cache_, force_prior_cache_);
        computeMassAndEnergyChangeRates(context, force_prior_cache_, mass_change_rates, energy_change_rates);
        context.accepted_dt = dt;
        const Scalar accepted_dt = dt;

        parallel::forRange(real_particle_count, [&](const std::size_t index_i)
        {
            ParticleState &particle = context.particles[index_i];
            particle.mom += force_rates_cache_[index_i] * accepted_dt;
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
        const std::size_t real_particle_count = context.particles.size();
        if (force_rates_cache_.size() != real_particle_count)
        {
            force_rates_cache_.assign(real_particle_count, Vec2::Zero());
        }

        std::vector<Scalar> mass_change_rates(real_particle_count, 0.0);
        std::vector<Scalar> energy_change_rates(real_particle_count, 0.0);
        computeMassAndEnergyChangeRates(context, force_prior_cache_, mass_change_rates, energy_change_rates);

        const Scalar accepted_dt = acceptedTimeStep(context, dt);
        parallel::forRange(real_particle_count, [&](const std::size_t index_i)
        {
            ParticleState &particle = context.particles[index_i];
            const Scalar mass_floor = 0.25 * reference_density_ * particle.volume;
            const Scalar internal_energy_floor = pressure_floor_ * particle.volume / (eos_.gamma() - 1.0);
            const Scalar kinetic_total = 0.5 * particle.mass * particle.vel.squaredNorm();
            particle.energy = std::max(
                kinetic_total + internal_energy_floor,
                particle.energy + energy_change_rates[index_i] * accepted_dt);
            particle.mass = std::max(mass_floor, particle.mass + mass_change_rates[index_i] * accepted_dt);
            particle.rho = particle.mass * safeInverse(particle.volume);
            particle.vel = particle.mom * safeInverse(std::max(particle.mass, kEpsilon));
            if (!particle.vel.allFinite())
            {
                particle.vel.setZero();
                particle.mom = particle.mass * particle.vel;
            }
            updatePrimitiveState(particle);
        });
    }

  private:
    void updateGhostKernelGradient(SimulationContext &context) const
    {
        const std::size_t real_particle_count = context.particles.size();
        if (context.ghost_particles.empty())
        {
            return;
        }

        parallel::forRange(real_particle_count, [&](const std::size_t index_i)
        {
            if (context.particles[index_i].indicator != 1)
            {
                return;
            }

            Vec2 gradient_summation = Vec2::Zero();
            for (const Neighbor &neighbor : context.neighbors[index_i])
            {
                if (neighbor.index >= real_particle_count)
                {
                    continue;
                }

                gradient_summation +=
                    neighbor.dweight_dr *
                    context.interaction_particles[neighbor.index].volume *
                    neighbor.direction;
            }

            if (gradient_summation.squaredNorm() <= kEpsilon)
            {
                return;
            }

            for (Neighbor &neighbor : context.neighbors[index_i])
            {
                if (neighbor.index < real_particle_count)
                {
                    continue;
                }

                const ParticleState &ghost_particle = context.interaction_particles[neighbor.index];
                const Vec2 ghost_gradient = -gradient_summation;
                const Scalar ghost_dweight_volume = -ghost_gradient.norm();
                neighbor.dweight_dr = ghost_dweight_volume * safeInverse(ghost_particle.volume);
                neighbor.direction = ghost_gradient / (ghost_dweight_volume + kEpsilon);
                const Vec2 displacement = ghost_particle.pos - context.particles[index_i].pos;
                neighbor.distance = std::max(std::abs(displacement.dot(neighbor.direction)), kEpsilon);
                if (context.kernel)
                {
                    neighbor.weight = context.kernel->weight(neighbor.distance);
                }
            }
        });
    }

    void computeMomentumChangeRates(const SimulationContext &context, std::vector<Vec2> &force_rates,
                                    std::vector<Vec2> &force_prior) const
    {
        const std::size_t real_particle_count = context.particles.size();
        const Scalar smoothing_length =
            context.kernel ? context.kernel->smoothingLength()
                           : context.config.spatial.particle_spacing * context.config.spatial.smoothing_length_ratio;

        parallel::forRange(real_particle_count, [&](const std::size_t index_i)
        {
            const ParticleState &particle_i = context.particles[index_i];
            const math::CompressibleState state_i{
                particle_i.rho,
                particle_i.vel,
                particle_i.p,
                particle_i.energy * safeInverse(std::max(particle_i.volume, kEpsilon))};

            Vec2 viscous_force = Vec2::Zero();
            for (const Neighbor &neighbor : context.neighbors[index_i])
            {
                const ParticleState &particle_j = context.interaction_particles[neighbor.index];
                const Mat2 &b_i = context.correction_matrices[index_i];
                const Mat2 b_j = neighbor.index < real_particle_count
                                     ? context.correction_matrices[neighbor.index]
                                     : Mat2::Identity();
                const Scalar correction_factor =
                    neighbor.direction.dot((b_i + b_j) * neighbor.direction);
                const Vec2 velocity_derivative =
                    (particle_i.vel - particle_j.vel) /
                    (neighbor.distance + 0.01 * smoothing_length);
                viscous_force +=
                    correction_factor * dynamic_viscosity_ * velocity_derivative *
                    neighbor.dweight_dr * particle_j.volume;
            }

            const Vec2 force_prior_rate = viscous_force * particle_i.volume;
            Vec2 momentum_change_rate = force_prior_rate;
            for (const Neighbor &neighbor : context.neighbors[index_i])
            {
                const ParticleState &particle_j = context.interaction_particles[neighbor.index];
                const math::CompressibleState state_j{
                    particle_j.rho,
                    particle_j.vel,
                    particle_j.p,
                    particle_j.energy * safeInverse(std::max(particle_j.volume, kEpsilon))};
                const math::CompressibleStarState interface_state =
                    riemann_.interfaceState(state_i, state_j, neighbor.direction);
                const Mat2 convect_flux =
                    interface_state.rho * interface_state.vel * interface_state.vel.transpose();
                const Scalar d_w_ij_v_j = neighbor.dweight_dr * particle_j.volume;

                momentum_change_rate -=
                    2.0 * particle_i.volume * d_w_ij_v_j *
                    (convect_flux + interface_state.p * Mat2::Identity()) *
                    neighbor.direction;
            }

            force_rates[index_i] = momentum_change_rate;
            force_prior[index_i] = force_prior_rate;
        });
    }

    void computeMassAndEnergyChangeRates(const SimulationContext &context, const std::vector<Vec2> &force_prior,
                                         std::vector<Scalar> &mass_change_rates,
                                         std::vector<Scalar> &energy_change_rates) const
    {
        const std::size_t real_particle_count = context.particles.size();
        parallel::forRange(real_particle_count, [&](const std::size_t index_i)
        {
            const ParticleState &particle_i = context.particles[index_i];
            const math::CompressibleState state_i{
                particle_i.rho,
                particle_i.vel,
                particle_i.p,
                particle_i.energy * safeInverse(std::max(particle_i.volume, kEpsilon))};

            Scalar mass_change_rate = 0.0;
            Scalar energy_change_rate = force_prior[index_i].dot(particle_i.vel);

            for (const Neighbor &neighbor : context.neighbors[index_i])
            {
                const ParticleState &particle_j = context.interaction_particles[neighbor.index];
                const math::CompressibleState state_j{
                    particle_j.rho,
                    particle_j.vel,
                    particle_j.p,
                    particle_j.energy * safeInverse(std::max(particle_j.volume, kEpsilon))};
                const math::CompressibleStarState interface_state =
                    riemann_.interfaceState(state_i, state_j, neighbor.direction);
                const Scalar d_w_ij_v_j = neighbor.dweight_dr * particle_j.volume;

                mass_change_rate -=
                    2.0 * particle_i.volume * d_w_ij_v_j *
                    (interface_state.rho * interface_state.vel).dot(neighbor.direction);
                energy_change_rate -=
                    2.0 * particle_i.volume * d_w_ij_v_j *
                    ((interface_state.energy + interface_state.p) * interface_state.vel).dot(neighbor.direction);
            }

            mass_change_rates[index_i] = mass_change_rate;
            energy_change_rates[index_i] = energy_change_rate;
        });
    }

    void updatePrimitiveState(ParticleState &particle) const
    {
        synchronizeConservativeVariables(particle);
        const Scalar kinetic_energy_density = 0.5 * particle.rho * particle.vel.squaredNorm();
        const Scalar total_energy_density = particle.energy * safeInverse(std::max(particle.volume, kEpsilon));
        const Scalar internal_energy_density = std::max(
            total_energy_density - kinetic_energy_density,
            pressure_floor_ / (eos_.gamma() - 1.0));
        particle.p = std::max(
            pressure_floor_,
            eos_.pressureFromDensityInternalEnergy(particle.rho, internal_energy_density));
        particle.energy = internal_energy_density * particle.volume + 0.5 * particle.mass * particle.vel.squaredNorm();
    }

    math::IdealGasEquationOfState eos_;
    math::HllcLimiterRiemannSolver riemann_;
    Scalar reference_density_;
    Scalar pressure_floor_;
    Scalar dynamic_viscosity_;
    AcousticTimeStepEstimator estimator_;
    LinearGradientCorrection correction_;
    std::vector<Vec2> force_rates_cache_;
    std::vector<Vec2> force_prior_cache_;
};
} // namespace dynamics
} // namespace eulerian_sph
