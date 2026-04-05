#include "eulerian_sph.h"

#include <cmath>
#include <iostream>
#include <memory>

using namespace eulerian_sph;

namespace
{
void initializeTaylorGreen(ParticleSet &particles, const SimulationConfig &config)
{
    const Scalar gamma = config.material.heat_capacity_ratio;
    const Scalar rho0 = config.material.reference_density;
    const Scalar c0 = config.material.sound_speed;

    for (ParticleState &particle : particles.data())
    {
        particle.rho = rho0;
        particle.p = c0 * c0 * particle.rho / gamma;
        particle.vel.x() = -std::cos(2.0 * kPi * particle.pos.x()) * std::sin(2.0 * kPi * particle.pos.y());
        particle.vel.y() = std::sin(2.0 * kPi * particle.pos.x()) * std::cos(2.0 * kPi * particle.pos.y());
        synchronizeConservativeVariables(particle);

        const Scalar rho_e = particle.p / (gamma - 1.0);
        particle.energy = rho_e * particle.volume + 0.5 * particle.mass * particle.vel.squaredNorm();
    }
}
} // namespace

int main()
{
    SimulationConfig config;
    config.case_name = "taylor_green_lg";
    config.spatial.domain_lower = Vec2(0.0, 0.0);
    config.spatial.domain_upper = Vec2(1.0, 1.0);
    config.spatial.particle_spacing = 1.0 / 50.0;
    config.material.reference_density = 1.0;
    config.material.sound_speed = 10.0;
    config.material.dynamic_viscosity = 0.01;
    config.material.heat_capacity_ratio = 1.4;
    config.time.acoustic_cfl = 0.6;
    config.time.end_time = 0.5;
    config.time.output_interval = 0.0;
    config.time.max_steps = 4000;
    config.output.directory = "cases/taylor_green/output";

    ParticleSet particles = makeCartesianLattice(config.spatial.domain_lower, config.spatial.domain_upper,
                                                 config.spatial.particle_spacing);
    initializeTaylorGreen(particles, config);

    auto kernel = std::make_shared<kernel::TabulatedKernel>(
        std::make_shared<kernel::LaguerreGaussKernel>(
            config.spatial.particle_spacing * config.spatial.smoothing_length_ratio),
        20);
    auto solver = std::make_unique<dynamics::CompressibleEulerianSolver>(config.material);

    Simulation simulation(config, std::move(particles), kernel, std::move(solver));
    simulation.addBoundary(std::make_unique<boundary::PeriodicBoundary>(
        boundary::Axis::x, config.spatial.domain_lower.x(), config.spatial.domain_upper.x()));
    simulation.addBoundary(std::make_unique<boundary::PeriodicBoundary>(
        boundary::Axis::y, config.spatial.domain_lower.y(), config.spatial.domain_upper.y()));
    simulation.addDatRecorder(std::make_unique<io::DatRecorder>(
        config.output.directory / (config.case_name + "_TotalKineticEnergy.dat"),
        "TotalKineticEnergy",
        [](const SimulationContext &context)
        {
            Scalar total_kinetic_energy = 0.0;
            for (const ParticleState &particle : context.particles.data())
            {
                total_kinetic_energy += 0.5 * particle.mass * particle.vel.squaredNorm();
            }
            return std::vector<Scalar>{total_kinetic_energy};
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

    std::cout << "Finished taylor_green_lg case. final_time="
              << simulation.context().current_time
              << " final_step="
              << simulation.context().current_step
              << std::endl;
    return 0;
}
