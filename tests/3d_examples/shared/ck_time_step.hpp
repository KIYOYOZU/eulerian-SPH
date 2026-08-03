#pragma once

#include "sphinxsys.h"

namespace SPH::channel_ck
{
using AdvectionViscousTimeStep =
    ReduceDynamicsCK<execution::ParallelPolicy, fluid_dynamics::AdvectionViscousTimeStepCK>;

template <class FluidType = WeaklyCompressibleFluid>
using AcousticTimeStep =
    ReduceDynamicsCK<execution::ParallelPolicy, fluid_dynamics::AcousticTimeStepCK<FluidType>>;
} // namespace SPH::channel_ck
