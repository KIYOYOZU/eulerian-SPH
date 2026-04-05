#pragma once

#include "../core/types.h"

#include <cmath>

namespace eulerian_sph
{
namespace math
{
class WeaklyCompressibleEquationOfState
{
  public:
    WeaklyCompressibleEquationOfState(const Scalar reference_density, const Scalar sound_speed,
                                      const Scalar reference_pressure = 0.0)
        : rho0_(reference_density), c0_(sound_speed), p0_(reference_pressure) {}

    Scalar pressureFromDensity(const Scalar density) const
    {
        return p0_ + c0_ * c0_ * (density - rho0_);
    }

    Scalar soundSpeed() const { return c0_; }

  private:
    Scalar rho0_;
    Scalar c0_;
    Scalar p0_;
};

class IdealGasEquationOfState
{
  public:
    explicit IdealGasEquationOfState(const Scalar heat_capacity_ratio) : gamma_(heat_capacity_ratio) {}

    Scalar pressureFromDensityInternalEnergy(const Scalar /*density*/, const Scalar internal_energy_density) const
    {
        return (gamma_ - 1.0) * internal_energy_density;
    }

    Scalar soundSpeed(const Scalar pressure, const Scalar density) const
    {
        return std::sqrt(std::max(0.0, gamma_ * pressure * safeInverse(density)));
    }

    Scalar gamma() const { return gamma_; }

  private:
    Scalar gamma_;
};
} // namespace math
} // namespace eulerian_sph
