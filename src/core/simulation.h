#pragma once

#include "simulation_context.h"

#include "../boundary/boundary_base.h"
#include "../dynamics/solver_base.h"
#include "../io/dat_recorder.h"

namespace eulerian_sph
{
class Simulation
{
  public:
    using ContextHook = std::function<void(SimulationContext &)>;

    Simulation(SimulationConfig config, ParticleSet particles, std::shared_ptr<kernel::KernelBase> kernel,
               std::unique_ptr<dynamics::EulerianSolver> solver);

    void addBoundary(std::unique_ptr<boundary::BoundaryCondition> boundary);
    void addDatRecorder(std::unique_ptr<io::DatRecorder> recorder);
    void addPostInitializeHook(ContextHook hook);
    void setWallParticles(WallParticleSet wall_particles);

    void initialize();
    void step();
    void run();

    SimulationContext &context() { return context_; }
    const SimulationContext &context() const { return context_; }

  private:
    Scalar resolveOutputInterval() const;
    void rebuildNeighbors();
    void writeFrame() const;

    SimulationContext context_;
    std::unique_ptr<dynamics::EulerianSolver> solver_;
    std::vector<std::unique_ptr<boundary::BoundaryCondition>> boundaries_;
    std::vector<std::unique_ptr<io::DatRecorder>> dat_recorders_;
    std::vector<ContextHook> post_initialize_hooks_;
    mutable std::size_t frame_index_ = 0;
    mutable std::size_t last_output_step_ = 0;
    Scalar next_output_time_ = 0.0;
};
} // namespace eulerian_sph
