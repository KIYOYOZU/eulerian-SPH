#pragma once

#include "../core/simulation_context.h"

namespace eulerian_sph
{
namespace boundary
{
using ParticleSelector = std::function<bool(const ParticleState &)>;
using ParticleUpdater = std::function<void(ParticleState &)>;
using GhostProjector = std::function<ParticleState(const ParticleState &)>;
using BoundaryTypeProvider = std::function<int(const ParticleState &)>;

class BoundaryCondition
{
  public:
    virtual ~BoundaryCondition() = default;

    virtual std::string name() const = 0;
    virtual void initialize(SimulationContext &context) { (void)context; }
    virtual void finalizeInitialization(SimulationContext &context) { (void)context; }
    virtual void prepareStep(SimulationContext &context) { (void)context; }
    virtual void applyBeforeStep(SimulationContext &context) { (void)context; }
    virtual void applyMidStep(SimulationContext &context) { (void)context; }
    virtual void applyAfterStep(SimulationContext &context) { (void)context; }
};
} // namespace boundary
} // namespace eulerian_sph
