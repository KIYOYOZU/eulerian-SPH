#include "eulerian_sph.h"

#include <cmath>
#include <iostream>
#include <memory>

using namespace eulerian_sph;

namespace
{
void initializeFreestream(ParticleSet &particles, const SimulationConfig &config, const Vec2 &velocity)
{
    for (ParticleState &particle : particles.data())
    {
        particle.rho = config.material.reference_density;
        particle.vel = velocity;
        synchronizeConservativeVariables(particle);
        particle.p = config.material.reference_pressure;
    }
}

WallParticleSet createCylinderWallParticles(const Vec2 &center, const Scalar radius, const Scalar spacing)
{
    WallParticleSet wall_particles;
    const Vec2 lower = center - Vec2(radius, radius);
    const Vec2 upper = center + Vec2(radius, radius);
    const Vec2 extent = upper - lower;
    const int cells_x = std::max(1, static_cast<int>(std::ceil(extent.x() / spacing)));
    const int cells_y = std::max(1, static_cast<int>(std::ceil(extent.y() / spacing)));

    for (int i = 0; i != cells_x; ++i)
    {
        for (int j = 0; j != cells_y; ++j)
        {
            const Vec2 position = lower + Vec2((static_cast<Scalar>(i) + 0.5) * spacing,
                                               (static_cast<Scalar>(j) + 0.5) * spacing);
            const Vec2 radial = position - center;
            if (radial.squaredNorm() > radius * radius)
            {
                continue;
            }

            WallParticle particle;
            particle.pos = position;
            particle.normal = radial.norm() <= kEpsilon ? Vec2::UnitX() : radial.normalized();
            particle.volume = spacing * spacing;
            wall_particles.addParticle(particle);
        }
    }

    return wall_particles;
}

Scalar distanceToOuterBox(const Vec2 &position, const Vec2 &lower, const Vec2 &upper, Vec2 &normal)
{
    const Scalar distance_left = std::abs(position.x() - lower.x());
    const Scalar distance_right = std::abs(upper.x() - position.x());
    const Scalar distance_bottom = std::abs(position.y() - lower.y());
    const Scalar distance_top = std::abs(upper.y() - position.y());

    Scalar minimum_distance = distance_left;
    normal = -Vec2::UnitX();
    if (distance_right < minimum_distance)
    {
        minimum_distance = distance_right;
        normal = Vec2::UnitX();
    }
    if (distance_bottom < minimum_distance)
    {
        minimum_distance = distance_bottom;
        normal = -Vec2::UnitY();
    }
    if (distance_top < minimum_distance)
    {
        minimum_distance = distance_top;
        normal = Vec2::UnitY();
    }
    return minimum_distance;
}
} // namespace

int main()
{
    const Scalar dl = 15.0;
    const Scalar dh = 10.0;
    const Scalar resolution = 1.0 / 4.0;
    const Scalar dl_sponge = resolution * 2.0;
    const Scalar dh_sponge = resolution * 2.0;
    const Vec2 cylinder_center(4.0, dh / 2.0);
    const Scalar cylinder_radius = 1.0;
    const Scalar cylinder_wall_spacing = resolution / 2.0;
    const Vec2 freestream_velocity(1.0, 0.0);

    SimulationConfig config;
    config.case_name = "flow_around_cylinder_lg";
    config.spatial.domain_lower = Vec2(-dl_sponge, -dh_sponge);
    config.spatial.domain_upper = Vec2(dl, dh + dh_sponge);
    config.spatial.particle_spacing = resolution;
    config.material.reference_density = 1.0;
    config.material.sound_speed = 10.0;
    config.material.dynamic_viscosity = 0.02;
    config.material.reference_pressure = 0.0;
    config.time.end_time = 80.0;
    config.time.output_interval = 0.0;
    config.time.max_steps = 20000;
    config.time.screen_output_interval = 1000;
    config.output.directory = "cases/flow_around_cylinder_lg/output";

    ParticleSet particles = makeCartesianLattice(
        config.spatial.domain_lower, config.spatial.domain_upper, config.spatial.particle_spacing,
        [&](const Vec2 &position)
        {
            return (position - cylinder_center).norm() > cylinder_radius;
        });
    initializeFreestream(particles, config, freestream_velocity);

    auto kernel = std::make_shared<kernel::TabulatedKernel>(
        std::make_shared<kernel::LaguerreGaussKernel>(
            config.spatial.particle_spacing * config.spatial.smoothing_length_ratio),
        20);
    auto solver = std::make_unique<dynamics::WeaklyCompressibleEulerianSolver>(config.material);

    Simulation simulation(config, std::move(particles), kernel, std::move(solver));
    simulation.setWallParticles(createCylinderWallParticles(
        cylinder_center, cylinder_radius, cylinder_wall_spacing));

    const auto water_block_normal = [&](const Vec2 &position) -> Vec2
    {
        Vec2 outer_normal = Vec2::UnitX();
        const Scalar outer_distance =
            distanceToOuterBox(position, config.spatial.domain_lower, config.spatial.domain_upper, outer_normal);
        const Vec2 radial = position - cylinder_center;
        const Scalar cylinder_distance = std::abs(radial.norm() - cylinder_radius);
        if (cylinder_distance < outer_distance && radial.norm() > kEpsilon)
        {
            return Vec2(-radial.normalized());
        }
        return outer_normal;
    };
    const auto water_block_signed_distance = [&](const Vec2 &position)
    {
        Vec2 outer_normal = Vec2::UnitX();
        const Scalar outer_distance =
            distanceToOuterBox(position, config.spatial.domain_lower, config.spatial.domain_upper, outer_normal);
        const Scalar cylinder_distance = std::abs((position - cylinder_center).norm() - cylinder_radius);
        return std::min(outer_distance, cylinder_distance);
    };
    const auto cylinder_normal = [&](const Vec2 &position) -> Vec2
    {
        const Vec2 radial = position - cylinder_center;
        return radial.norm() <= kEpsilon ? Vec2::UnitX() : Vec2(radial.normalized());
    };
    const auto cylinder_signed_distance = [&](const Vec2 &position)
    {
        return std::abs((position - cylinder_center).norm() - cylinder_radius);
    };
    simulation.addPostInitializeHook(
        [&](SimulationContext &context)
        {
            dynamics::updateParticleShapeState(context.particles, water_block_normal, water_block_signed_distance);
            dynamics::updateWallShapeState(context.wall_particles, cylinder_normal, cylinder_signed_distance);
            dynamics::updateFreeSurfaceIndicator(context, true);
            dynamics::updateSmearedSurfaceIndicator(context);
        });

    ParticleState farfield_state;
    farfield_state.rho = config.material.reference_density;
    farfield_state.vel = freestream_velocity;
    farfield_state.p = config.material.reference_pressure;
    const Scalar farfield_band = 2.0 * config.spatial.particle_spacing;

    simulation.addBoundary(std::make_unique<boundary::FarFieldBoundary>(
        "far_field",
        [&](const ParticleState &particle)
        {
            Vec2 outer_normal = Vec2::UnitX();
            const Scalar outer_distance =
                distanceToOuterBox(particle.pos, config.spatial.domain_lower, config.spatial.domain_upper, outer_normal);
            return outer_distance <= farfield_band;
        },
        [](const ParticleState &particle) { return particle.normal; },
        farfield_state,
        config.material.sound_speed,
        true,
        true));
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

    std::cout << "Finished flow_around_cylinder_lg case. final_time="
              << simulation.context().current_time
              << " final_step="
              << simulation.context().current_step
              << std::endl;
    return 0;
}
