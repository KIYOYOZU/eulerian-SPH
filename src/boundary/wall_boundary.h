#pragma once

#include "boundary_base.h"

namespace eulerian_sph
{
namespace boundary
{
class WallBoundary : public BoundaryCondition
{
  public:
    WallBoundary(std::string boundary_name, ParticleSelector selector, ParticleUpdater updater)
        : boundary_name_(std::move(boundary_name)),
          selector_(std::move(selector)),
          updater_(std::move(updater)) {}

    std::string name() const override { return boundary_name_; }

    void applyAfterStep(SimulationContext &context) override
    {
        for (ParticleState &particle : context.particles.data())
        {
            if (selector_(particle))
            {
                updater_(particle);
                synchronizeConservativeVariables(particle);
            }
        }
    }

  private:
    std::string boundary_name_;
    ParticleSelector selector_;
    ParticleUpdater updater_;
};
} // namespace boundary
} // namespace eulerian_sph
