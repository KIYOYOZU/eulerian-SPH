#include "eulerian_sph.h"

#include <cmath>
#include <iostream>
#include <memory>

using namespace eulerian_sph;

namespace
{
ParticleState makeSupersonicReference(const SimulationConfig &config, const Scalar mach_number)
{
    ParticleState state;
    state.rho = config.material.reference_density;
    state.p = config.material.reference_pressure;
    const Scalar sound_speed = std::sqrt(config.material.heat_capacity_ratio * state.p / state.rho);
    state.vel = Vec2(mach_number * sound_speed, 0.0);
    synchronizeConservativeVariables(state);
    const Scalar rho_e = state.p / (config.material.heat_capacity_ratio - 1.0);
    state.energy = rho_e + 0.5 * state.rho * state.vel.squaredNorm();
    return state;
}

void initializeSupersonicField(ParticleSet &particles, const SimulationConfig &config, const Scalar mach_number)
{
    ParticleState reference = makeSupersonicReference(config, mach_number);
    for (ParticleState &particle : particles.data())
    {
        particle.rho = reference.rho;
        particle.p = reference.p;
        particle.vel = reference.vel;
        synchronizeConservativeVariables(particle);
        particle.energy = reference.p / (config.material.heat_capacity_ratio - 1.0) * particle.volume +
                          0.5 * particle.mass * particle.vel.squaredNorm();
    }
}

Scalar distanceToFluidBlockBoundary(const Vec2 &position, const Vec2 &calculation_center, const Scalar calculation_radius,
                                    const Vec2 &insert_center, const Scalar insert_radius, Vec2 &normal)
{
    const Vec2 outer_radial = position - calculation_center;
    const Scalar outer_distance = std::abs(outer_radial.norm() - calculation_radius);
    Vec2 boundary_normal =
        outer_radial.norm() <= kEpsilon ? Vec2(-1.0, 0.0) : Vec2(outer_radial.normalized());
    Scalar minimum_distance = outer_distance;

    const Vec2 insert_radial = position - insert_center;
    const Scalar insert_distance = std::abs(insert_radial.norm() - insert_radius);
    if (insert_distance < minimum_distance)
    {
        minimum_distance = insert_distance;
        boundary_normal =
            insert_radial.norm() <= kEpsilon ? Vec2(-1.0, 0.0) : Vec2(-insert_radial.normalized());
    }

    const Scalar diameter_distance = std::abs(calculation_center.x() - position.x());
    if (diameter_distance < minimum_distance)
    {
        minimum_distance = diameter_distance;
        boundary_normal = Vec2::UnitX();
    }

    normal = boundary_normal;
    return minimum_distance;
}
} // namespace

int main()
{
    const Scalar spacing = 1.0 / 7.0;
    const Scalar calculation_radius = 11.0;
    const Scalar radius_with_band = calculation_radius + 4.0 * spacing;
    const Vec2 calculation_center(calculation_radius, 0.0);
    const Vec2 insert_center(7.0, 0.0);
    const Scalar insert_radius = 1.0;
    const Scalar mach_number = 2.0;

    SimulationConfig config;
    config.case_name = "supersonic_flow_new_bc";
    config.spatial.domain_lower = Vec2(0.0, -radius_with_band);
    config.spatial.domain_upper = Vec2(radius_with_band, radius_with_band);
    config.spatial.particle_spacing = spacing;
    config.spatial.ghost_band_width = 2.0 * spacing;
    config.material.reference_density = 1.0;
    config.material.reference_pressure = 1.0 / 1.4;
    config.material.sound_speed = 1.0;
    config.material.heat_capacity_ratio = 1.4;
    config.time.acoustic_cfl = 0.1;
    config.time.end_time = 40.0;
    config.time.output_interval = 0.0;
    config.time.max_steps = 40000;
    config.time.screen_output_interval = 500;
    config.output.directory = "cases/supersonic_flow_new_bc/output";

    ParticleSet particles = makeCartesianLattice(
        config.spatial.domain_lower, config.spatial.domain_upper, config.spatial.particle_spacing,
        [&](const Vec2 &position)
        {
            const bool inside_outer_circle = (position - calculation_center).norm() <= calculation_radius;
            const bool outside_insert = (position - insert_center).norm() >= insert_radius;
            const bool left_half_disk = position.x() <= calculation_center.x();
            return inside_outer_circle && outside_insert && left_half_disk;
        });
    initializeSupersonicField(particles, config, mach_number);

    auto kernel = std::make_shared<kernel::TabulatedKernel>(
        std::make_shared<kernel::LaguerreGaussKernel>(
            config.spatial.particle_spacing * config.spatial.smoothing_length_ratio),
        20);
    auto solver = std::make_unique<dynamics::CompressibleEulerianSolver>(config.material);

    Simulation simulation(config, std::move(particles), kernel, std::move(solver));
    const auto fluid_block_normal = [&](const Vec2 &position)
    {
        Vec2 normal = Vec2::UnitX();
        distanceToFluidBlockBoundary(position, calculation_center, calculation_radius, insert_center, insert_radius, normal);
        return normal;
    };
    const auto fluid_block_signed_distance = [&](const Vec2 &position)
    {
        Vec2 normal = Vec2::UnitX();
        return distanceToFluidBlockBoundary(position, calculation_center, calculation_radius, insert_center, insert_radius, normal);
    };
    simulation.addPostInitializeHook(
        [&](SimulationContext &context)
        {
            dynamics::updateParticleShapeState(context.particles, fluid_block_normal, fluid_block_signed_distance);
            dynamics::updateFreeSurfaceIndicator(context, false);
        });

    ParticleState farfield_state = makeSupersonicReference(config, mach_number);
    simulation.addBoundary(std::make_unique<boundary::GhostBoundary>(
        "supersonic_ghost_boundary",
        [](const ParticleState &) { return true; },
        [&](const ParticleState &particle) { return particle.normal; },
        [&](const ParticleState &particle) { return particle.signed_distance; },
        [&](ParticleState &ghost, const SimulationContext &context, const std::size_t index_i, const Vec2 &ghost_normal)
        {
            const ParticleState &particle = context.particles[index_i];
            const bool reflective_wall = particle.boundary_type == 3;
            const Scalar sound_speed =
                std::sqrt(config.material.heat_capacity_ratio * farfield_state.p / farfield_state.rho);

            auto updateEnergy = [&](ParticleState &state)
            {
                synchronizeConservativeVariables(state);
                state.energy = state.p / (config.material.heat_capacity_ratio - 1.0) * state.volume +
                               0.5 * state.mass * state.vel.squaredNorm();
            };

            if (reflective_wall)
            {
                ghost.rho = particle.rho;
                ghost.p = particle.p;
                ghost.vel = particle.vel - 2.0 * particle.vel.dot(ghost_normal) * ghost_normal;
                updateEnergy(ghost);
                return;
            }

            const Scalar velocity_farfield_normal = farfield_state.vel.dot(ghost_normal);
            const Scalar velocity_boundary_normal = particle.vel.dot(ghost_normal);
            const bool treat_as_inflow =
                ghost_normal.x() <= 0.0 || std::abs(ghost_normal.y()) > std::abs(ghost_normal.x());

            Scalar inner_weight_summation = context.kernel ? context.kernel->weight(0.0) * particle.volume : particle.volume;
            Scalar rho_summation = 0.0;
            Scalar p_summation = 0.0;
            Scalar vel_normal_summation = 0.0;
            Vec2 vel_tangential_summation = Vec2::Zero();
            std::size_t total_inner_neighbors = 0;

            for (const Neighbor &neighbor : context.neighbors[index_i])
            {
                if (neighbor.index >= context.particles.size())
                {
                    continue;
                }

                const ParticleState &neighbor_particle = context.particles[neighbor.index];
                inner_weight_summation += neighbor.weight * neighbor_particle.volume;
                rho_summation += neighbor_particle.rho;
                p_summation += neighbor_particle.p;
                vel_normal_summation += neighbor_particle.vel.dot(ghost_normal);
                vel_tangential_summation +=
                    neighbor_particle.vel - neighbor_particle.vel.dot(ghost_normal) * ghost_normal;
                total_inner_neighbors += 1;
            }

            if (total_inner_neighbors == 0)
            {
                ghost.rho = farfield_state.rho;
                ghost.p = farfield_state.p;
                ghost.vel = farfield_state.vel;
                updateEnergy(ghost);
                return;
            }

            const Scalar inverse_count = 1.0 / static_cast<Scalar>(total_inner_neighbors);
            const Scalar rho_average = rho_summation * inverse_count;
            const Scalar p_average = p_summation * inverse_count;
            const Scalar vel_normal_average = vel_normal_summation * inverse_count;
            const Vec2 vel_tangential_average = vel_tangential_summation * inverse_count;

            if (treat_as_inflow)
            {
                if (std::abs(velocity_boundary_normal) >= sound_speed)
                {
                    ghost.rho = farfield_state.rho;
                    ghost.p = farfield_state.p;
                    ghost.vel = farfield_state.vel;
                }
                else
                {
                    ghost.p = p_average * inner_weight_summation + farfield_state.p * (1.0 - inner_weight_summation);
                    ghost.rho = rho_average * inner_weight_summation + farfield_state.rho * (1.0 - inner_weight_summation);
                    const Scalar vel_normal =
                        vel_normal_average * inner_weight_summation +
                        velocity_farfield_normal * (1.0 - inner_weight_summation);
                    ghost.vel = vel_normal * ghost_normal +
                                (farfield_state.vel - velocity_farfield_normal * ghost_normal);
                }
            }
            else
            {
                if (std::abs(velocity_boundary_normal) >= sound_speed)
                {
                    ghost.rho = particle.rho;
                    ghost.p = particle.p;
                    ghost.vel = particle.vel;
                }
                else
                {
                    ghost.p = p_average * inner_weight_summation + farfield_state.p * (1.0 - inner_weight_summation);
                    ghost.rho = rho_average * inner_weight_summation + farfield_state.rho * (1.0 - inner_weight_summation);
                    const Scalar vel_normal =
                        vel_normal_average * inner_weight_summation +
                        velocity_farfield_normal * (1.0 - inner_weight_summation);
                    ghost.vel = vel_normal * ghost_normal + vel_tangential_average;
                }
            }

            updateEnergy(ghost);
        },
        [&](const ParticleState &particle)
        {
            return (particle.pos - insert_center).norm() <= insert_radius + 5.0 * spacing ? 3 : 9;
        }));
    simulation.addDatRecorder(std::make_unique<io::DatRecorder>(
        config.output.directory / (config.case_name + "_MaximumSpeed.dat"),
        "MaximumSpeed",
        [](const SimulationContext &context)
        {
            Scalar maximum_speed = 0.0;
            for (const ParticleState &particle : context.particles.data())
            {
                maximum_speed = std::max(maximum_speed, particle.vel.norm());
            }
            return std::vector<Scalar>{maximum_speed};
        }));

    simulation.run();

    std::cout << "Finished supersonic_flow_new_bc case. final_time="
              << simulation.context().current_time
              << " final_step="
              << simulation.context().current_step
              << std::endl;
    return 0;
}
