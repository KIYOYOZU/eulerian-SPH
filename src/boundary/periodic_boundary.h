#pragma once

#include "boundary_base.h"

namespace eulerian_sph
{
namespace boundary
{
enum class Axis
{
    x = 0,
    y = 1
};

class PeriodicBoundary : public BoundaryCondition
{
  public:
    PeriodicBoundary(const Axis axis, const Scalar lower, const Scalar upper)
        : axis_(axis), lower_(lower), upper_(upper), span_(upper - lower) {}

    std::string name() const override { return "PeriodicBoundary"; }

    void prepareStep(SimulationContext &context) override
    {
        const int component = axis_ == Axis::x ? 0 : 1;
        context.periodic_axes[component] = true;
        context.periodic_lower[component] = lower_;
        context.periodic_upper[component] = upper_;
    }

    void applyAfterStep(SimulationContext &context) override
    {
        const int component = axis_ == Axis::x ? 0 : 1;
        for (ParticleState &particle : context.particles.data())
        {
            while (particle.pos[component] < lower_)
            {
                particle.pos[component] += span_;
            }
            while (particle.pos[component] > upper_)
            {
                particle.pos[component] -= span_;
            }
        }
    }

  private:
    Axis axis_;
    Scalar lower_;
    Scalar upper_;
    Scalar span_;
};
} // namespace boundary
} // namespace eulerian_sph
