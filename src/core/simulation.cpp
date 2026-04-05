#include "simulation.h"

#include "parallel.h"

#include "../io/csv_writer.h"
#include "../io/vtk_writer.h"

#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace eulerian_sph
{
Simulation::Simulation(SimulationConfig config, ParticleSet particles, std::shared_ptr<kernel::KernelBase> kernel,
                       std::unique_ptr<dynamics::EulerianSolver> solver)
    : solver_(std::move(solver))
{
    context_.config = std::move(config);
    context_.particles = std::move(particles);
    context_.kernel = std::move(kernel);
    const Scalar cell_size = context_.kernel ? context_.kernel->supportRadius()
                                             : context_.config.spatial.particle_spacing * context_.config.spatial.smoothing_length_ratio;
    context_.grid.setCellSize(cell_size);
}

void Simulation::addBoundary(std::unique_ptr<boundary::BoundaryCondition> boundary)
{
    boundaries_.push_back(std::move(boundary));
}

void Simulation::addDatRecorder(std::unique_ptr<io::DatRecorder> recorder)
{
    dat_recorders_.push_back(std::move(recorder));
}

void Simulation::addPostInitializeHook(ContextHook hook)
{
    post_initialize_hooks_.push_back(std::move(hook));
}

void Simulation::setWallParticles(WallParticleSet wall_particles)
{
    context_.wall_particles = std::move(wall_particles);
}

void Simulation::initialize()
{
    if (!solver_)
    {
        throw std::runtime_error("Simulation requires a valid solver.");
    }

    std::filesystem::create_directories(context_.config.output.directory);
    context_.ghost_particles.clear();
    context_.ghost_relations.clear();
    context_.ghost_boundary_pairs.clear();
    context_.ghost_projector = {};
    context_.mid_step_callback = {};
    context_.periodic_axes = {false, false};
    for (auto &boundary : boundaries_)
    {
        boundary->initialize(context_);
        boundary->prepareStep(context_);
    }
    rebuildNeighbors();
    for (const auto &hook : post_initialize_hooks_)
    {
        hook(context_);
    }
    for (auto &boundary : boundaries_)
    {
        boundary->finalizeInitialization(context_);
    }
    for (auto &boundary : boundaries_)
    {
        boundary->applyBeforeStep(context_);
    }
    if (!context_.ghost_particles.empty())
    {
        rebuildNeighbors();
    }
    solver_->initialize(context_);
    next_output_time_ = resolveOutputInterval();
    writeFrame();
}

void Simulation::step()
{
    context_.ghost_particles.clear();
    context_.ghost_relations.clear();
    context_.ghost_boundary_pairs.clear();
    context_.ghost_projector = {};
    context_.mid_step_callback = {};
    context_.periodic_axes = {false, false};
    for (auto &boundary : boundaries_)
    {
        boundary->prepareStep(context_);
    }
    rebuildNeighbors();
    for (auto &boundary : boundaries_)
    {
        boundary->applyBeforeStep(context_);
    }
    if (!context_.ghost_particles.empty())
    {
        rebuildNeighbors();
    }

    const Scalar dt = solver_->computeTimeStep(context_);
    if (solver_->supportsSplitAdvance())
    {
        solver_->advanceFirstHalf(context_, dt);
        for (auto &boundary : boundaries_)
        {
            boundary->applyMidStep(context_);
        }
        solver_->advanceSecondHalf(context_, dt);
    }
    else
    {
        solver_->advance(context_, dt);
    }
    const Scalar accepted_dt = solver_->acceptedTimeStep(context_, dt);
    context_.current_time += accepted_dt;
    context_.current_step += 1;

    for (auto &boundary : boundaries_)
    {
        boundary->applyAfterStep(context_);
    }
}

void Simulation::run()
{
    initialize();

    while (context_.current_time < context_.config.time.end_time &&
           context_.current_step < context_.config.time.max_steps)
    {
        step();
        if (context_.config.time.screen_output_interval > 0 &&
            context_.current_step % context_.config.time.screen_output_interval == 0)
        {
            const Scalar maximum_speed = parallel::reduceMax(
                context_.particles.size(),
                0.0,
                [&](const std::size_t index_i)
                {
                    return context_.particles[index_i].vel.norm();
                });
            std::cout << "step=" << context_.current_step
                      << " time=" << context_.current_time
                      << " dt=" << context_.accepted_dt
                      << " max_speed=" << maximum_speed
                      << std::endl;
        }
        while (next_output_time_ > 0.0 &&
               context_.current_time + kEpsilon >= next_output_time_)
        {
            writeFrame();
            next_output_time_ += resolveOutputInterval();
        }
    }

    if (context_.current_step != last_output_step_)
    {
        writeFrame();
    }
}

Scalar Simulation::resolveOutputInterval() const
{
    if (context_.config.time.output_interval > 0.0)
    {
        return context_.config.time.output_interval;
    }

    if (!context_.config.output.write_vtp || context_.config.output.vtp_output_count <= 1 ||
        context_.config.time.end_time <= 0.0)
    {
        return 0.0;
    }

    return context_.config.time.end_time /
           static_cast<Scalar>(context_.config.output.vtp_output_count - 1);
}

void Simulation::rebuildNeighbors()
{
    const Scalar search_radius = context_.kernel ? context_.kernel->supportRadius()
                                                 : context_.config.spatial.particle_spacing * context_.config.spatial.smoothing_length_ratio;
    context_.interaction_particles = context_.particles;
    for (const ParticleState &ghost_particle : context_.ghost_particles.data())
    {
        context_.interaction_particles.addParticle(ghost_particle);
    }

    if (context_.periodic_axes[0] || context_.periodic_axes[1])
    {
        const std::size_t real_particle_count = context_.particles.size();
        context_.neighbors.assign(real_particle_count, {});
        const Scalar radius_squared = search_radius * search_radius;
        const Vec2 periodic_span = context_.periodic_upper - context_.periodic_lower;

        const LoopIndex particle_count = static_cast<LoopIndex>(real_particle_count);
#if defined(EULERIAN_SPH_USE_OPENMP)
#pragma omp parallel for schedule(static) if(shouldParallelize(real_particle_count))
#endif
        for (LoopIndex index_i = 0; index_i < particle_count; ++index_i)
        {
            auto &particle_neighbors = context_.neighbors[static_cast<std::size_t>(index_i)];
            for (std::size_t index_j = 0; index_j != context_.interaction_particles.size(); ++index_j)
            {
                if (static_cast<std::size_t>(index_i) == index_j)
                {
                    continue;
                }

                Vec2 displacement =
                    context_.interaction_particles[index_j].pos - context_.particles[static_cast<std::size_t>(index_i)].pos;
                for (int axis = 0; axis != 2; ++axis)
                {
                    if (!context_.periodic_axes[axis])
                    {
                        continue;
                    }

                    const Scalar span = periodic_span[axis];
                    if (displacement[axis] > 0.5 * span)
                    {
                        displacement[axis] -= span;
                    }
                    else if (displacement[axis] < -0.5 * span)
                    {
                        displacement[axis] += span;
                    }
                }

                const Scalar distance_squared = displacement.squaredNorm();
                if (distance_squared > radius_squared || distance_squared <= kEpsilon)
                {
                    continue;
                }

                Neighbor neighbor;
                neighbor.index = index_j;
                neighbor.distance = std::sqrt(distance_squared);
                neighbor.direction = displacement / neighbor.distance;
                particle_neighbors.push_back(neighbor);
            }
        }
    }
    else
    {
        context_.grid.rebuild(context_.interaction_particles);
        context_.grid.queryNeighbors(context_.interaction_particles, search_radius, context_.neighbors, context_.particles.size());
    }

    context_.wall_neighbors.assign(context_.particles.size(), {});
    if (!context_.wall_particles.empty())
    {
        const Scalar radius_squared = search_radius * search_radius;
        const std::size_t particle_count = context_.particles.size();
        const LoopIndex loop_particle_count = static_cast<LoopIndex>(particle_count);
#if defined(EULERIAN_SPH_USE_OPENMP)
#pragma omp parallel for schedule(static) if(shouldParallelize(particle_count))
#endif
        for (LoopIndex index_i = 0; index_i < loop_particle_count; ++index_i)
        {
            auto &particle_neighbors = context_.wall_neighbors[static_cast<std::size_t>(index_i)];
            for (std::size_t index_j = 0; index_j != context_.wall_particles.size(); ++index_j)
            {
                const Vec2 displacement =
                    context_.wall_particles[index_j].pos - context_.particles[static_cast<std::size_t>(index_i)].pos;
                const Scalar distance_squared = displacement.squaredNorm();
                if (distance_squared > radius_squared || distance_squared <= kEpsilon)
                {
                    continue;
                }

                Neighbor neighbor;
                neighbor.index = index_j;
                neighbor.distance = std::sqrt(distance_squared);
                neighbor.direction = displacement / neighbor.distance;
                particle_neighbors.push_back(neighbor);
            }
        }
    }
    else
    {
        context_.wall_inner_neighbors.clear();
        context_.wall_contact_neighbors.clear();
        context_.wall_correction_matrices.clear();
    }

    if (!context_.wall_particles.empty())
    {
        const Scalar radius_squared = search_radius * search_radius;
        const std::size_t wall_count = context_.wall_particles.size();
        context_.wall_inner_neighbors.assign(wall_count, {});
        context_.wall_contact_neighbors.assign(wall_count, {});
        const LoopIndex loop_wall_count = static_cast<LoopIndex>(wall_count);
#if defined(EULERIAN_SPH_USE_OPENMP)
#pragma omp parallel for schedule(static) if(shouldParallelize(wall_count))
#endif
        for (LoopIndex index_i = 0; index_i < loop_wall_count; ++index_i)
        {
            auto &wall_inner_neighbors = context_.wall_inner_neighbors[static_cast<std::size_t>(index_i)];
            auto &wall_contact_neighbors = context_.wall_contact_neighbors[static_cast<std::size_t>(index_i)];

            for (std::size_t index_j = 0; index_j != wall_count; ++index_j)
            {
                if (static_cast<std::size_t>(index_i) == index_j)
                {
                    continue;
                }

                const Vec2 displacement =
                    context_.wall_particles[index_j].pos - context_.wall_particles[static_cast<std::size_t>(index_i)].pos;
                const Scalar distance_squared = displacement.squaredNorm();
                if (distance_squared > radius_squared || distance_squared <= kEpsilon)
                {
                    continue;
                }

                Neighbor neighbor;
                neighbor.index = index_j;
                neighbor.distance = std::sqrt(distance_squared);
                neighbor.direction = displacement / neighbor.distance;
                wall_inner_neighbors.push_back(neighbor);
            }

            for (std::size_t index_j = 0; index_j != context_.particles.size(); ++index_j)
            {
                const Vec2 displacement =
                    context_.particles[index_j].pos - context_.wall_particles[static_cast<std::size_t>(index_i)].pos;
                const Scalar distance_squared = displacement.squaredNorm();
                if (distance_squared > radius_squared || distance_squared <= kEpsilon)
                {
                    continue;
                }

                Neighbor neighbor;
                neighbor.index = index_j;
                neighbor.distance = std::sqrt(distance_squared);
                neighbor.direction = displacement / neighbor.distance;
                wall_contact_neighbors.push_back(neighbor);
            }
        }
    }

    if (!context_.kernel)
    {
        return;
    }

    const LoopIndex real_particle_count = static_cast<LoopIndex>(context_.neighbors.size());
#if defined(EULERIAN_SPH_USE_OPENMP)
#pragma omp parallel for schedule(static) if(shouldParallelize(context_.neighbors.size()))
#endif
    for (LoopIndex index_i = 0; index_i < real_particle_count; ++index_i)
    {
        auto &particle_neighbors = context_.neighbors[static_cast<std::size_t>(index_i)];
        for (Neighbor &neighbor : particle_neighbors)
        {
            neighbor.weight = context_.kernel->weight(neighbor.distance);
            neighbor.dweight_dr = context_.kernel->dWeightDr(neighbor.distance);
        }
    }

    const LoopIndex wall_neighbor_count = static_cast<LoopIndex>(context_.wall_neighbors.size());
#if defined(EULERIAN_SPH_USE_OPENMP)
#pragma omp parallel for schedule(static) if(shouldParallelize(context_.wall_neighbors.size()))
#endif
    for (LoopIndex index_i = 0; index_i < wall_neighbor_count; ++index_i)
    {
        auto &particle_neighbors = context_.wall_neighbors[static_cast<std::size_t>(index_i)];
        for (Neighbor &neighbor : particle_neighbors)
        {
            neighbor.weight = context_.kernel->weight(neighbor.distance);
            neighbor.dweight_dr = context_.kernel->dWeightDr(neighbor.distance);
        }
    }

    const LoopIndex wall_inner_count = static_cast<LoopIndex>(context_.wall_inner_neighbors.size());
#if defined(EULERIAN_SPH_USE_OPENMP)
#pragma omp parallel for schedule(static) if(shouldParallelize(context_.wall_inner_neighbors.size()))
#endif
    for (LoopIndex index_i = 0; index_i < wall_inner_count; ++index_i)
    {
        auto &particle_neighbors = context_.wall_inner_neighbors[static_cast<std::size_t>(index_i)];
        for (Neighbor &neighbor : particle_neighbors)
        {
            neighbor.weight = context_.kernel->weight(neighbor.distance);
            neighbor.dweight_dr = context_.kernel->dWeightDr(neighbor.distance);
        }
    }

    const LoopIndex wall_contact_count = static_cast<LoopIndex>(context_.wall_contact_neighbors.size());
#if defined(EULERIAN_SPH_USE_OPENMP)
#pragma omp parallel for schedule(static) if(shouldParallelize(context_.wall_contact_neighbors.size()))
#endif
    for (LoopIndex index_i = 0; index_i < wall_contact_count; ++index_i)
    {
        auto &particle_neighbors = context_.wall_contact_neighbors[static_cast<std::size_t>(index_i)];
        for (Neighbor &neighbor : particle_neighbors)
        {
            neighbor.weight = context_.kernel->weight(neighbor.distance);
            neighbor.dweight_dr = context_.kernel->dWeightDr(neighbor.distance);
        }
    }
}

void Simulation::writeFrame() const
{
    const std::string stem = context_.config.case_name + "_" + std::to_string(frame_index_);
    if (context_.config.output.write_vtp)
    {
        io::writeVtkFrame(context_, context_.config.output.directory / (stem + ".vtp"));
    }
    for (const auto &recorder : dat_recorders_)
    {
        recorder->write(context_);
    }
    last_output_step_ = context_.current_step;
    frame_index_ += 1;
}
} // namespace eulerian_sph
