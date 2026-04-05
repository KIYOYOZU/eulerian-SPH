#pragma once

#include "eos.h"

#include <algorithm>
#include <cmath>

namespace eulerian_sph
{
namespace math
{
struct AcousticState
{
    Scalar rho = 1.0;
    Vec2 vel = Vec2::Zero();
    Scalar p = 0.0;
};

struct AcousticStarState : AcousticState
{
};

class AcousticRiemannSolver
{
  public:
    AcousticRiemannSolver(const Scalar reference_density, const Scalar sound_speed, const Scalar limiter_slope = 3.0)
        : rho0_(reference_density),
          c0_(sound_speed),
          rho0c0_(rho0_ * c0_),
          inv_rho0c0_ave_(1.0 / rho0c0_),
          rho0c0_geo_ave_(rho0c0_),
          limiter_ref_(reference_density * inv_rho0c0_ave_),
          limiter_slope_(limiter_slope) {}

    AcousticStarState interfaceState(const AcousticState &left, const AcousticState &right, const Vec2 &normal) const
    {
        const Scalar rho_star = 0.5 * (left.rho + right.rho);
        const Scalar p_average = 0.5 * (left.p + right.p);
        const Vec2 velocity_average = 0.5 * (left.vel + right.vel);
        const Scalar ul = -normal.dot(left.vel);
        const Scalar ur = -normal.dot(right.vel);
        const Scalar u_jump = ul - ur;
        const Scalar limited_mach_number = limiter(std::max(u_jump, Scalar(0)));

        const Scalar p_star = p_average + 0.5 * rho0c0_geo_ave_ * u_jump * limited_mach_number;
        const Scalar u_dissipative =
            0.5 * (left.p - right.p) * inv_rho0c0_ave_ * limited_mach_number * limited_mach_number;

        AcousticStarState star;
        star.rho = rho_star;
        star.p = p_star;
        star.vel = velocity_average - normal * u_dissipative;
        return star;
    }

  private:
    Scalar limiter(const Scalar measure) const
    {
        return std::min(limiter_slope_ * measure * limiter_ref_, Scalar(1));
    }

    Scalar rho0_;
    Scalar c0_;
    Scalar rho0c0_;
    Scalar inv_rho0c0_ave_;
    Scalar rho0c0_geo_ave_;
    Scalar limiter_ref_;
    Scalar limiter_slope_;
};

struct CompressibleState
{
    Scalar rho = 1.0;
    Vec2 vel = Vec2::Zero();
    Scalar p = 0.0;
    Scalar energy = 0.0;
};

struct CompressibleStarState : CompressibleState
{
};

class HllcLimiterRiemannSolver
{
  public:
    HllcLimiterRiemannSolver(const IdealGasEquationOfState eos, const Scalar limiter_parameter)
        : eos_(eos), limiter_parameter_(limiter_parameter) {}

    CompressibleStarState interfaceState(const CompressibleState &left, const CompressibleState &right,
                                         const Vec2 &normal) const
    {
        const Scalar ul = -normal.dot(left.vel);
        const Scalar ur = -normal.dot(right.vel);

        const Vec2 vl = left.vel - ul * (-normal);
        const Vec2 vr = right.vel - ur * (-normal);
        const Scalar density_ratio = right.rho * safeInverse(left.rho);
        const Scalar u_tilde = (ul + ur * density_ratio) / (1.0 + density_ratio);
        const Scalar v_tilde = (vl.norm() + vr.norm() * density_ratio) / (1.0 + density_ratio);
        const Scalar hl = (left.energy + left.p) * safeInverse(left.rho);
        const Scalar hr = (right.energy + right.p) * safeInverse(right.rho);
        const Scalar h_tilde = (hl + hr * density_ratio) / (1.0 + density_ratio);
        const Scalar sound_tilde = std::sqrt(
            std::max(0.0, (eos_.gamma() - 1.0) * (h_tilde - 0.5 * (u_tilde * u_tilde + v_tilde * v_tilde))));

        const Scalar s_l = std::min(ul - eos_.soundSpeed(left.p, left.rho), u_tilde - sound_tilde);
        const Scalar s_r = std::max(ur + eos_.soundSpeed(right.p, right.rho), u_tilde + sound_tilde);

        const Scalar rhol_cl = eos_.soundSpeed(left.p, left.rho) * left.rho;
        const Scalar rhor_cr = eos_.soundSpeed(right.p, right.rho) * right.rho;
        const Scalar clr = (rhol_cl + rhor_cr) * safeInverse(left.rho + right.rho);
        const Scalar limiter = std::min(limiter_parameter_ * std::max((ul - ur) * safeInverse(clr), 0.0), 1.0);

        const Scalar denominator = left.rho * (s_l - ul) - right.rho * (s_r - ur);
        const Scalar s_star = (right.p - left.p) * limiter * limiter * safeInverse(denominator) +
                              (left.rho * (s_l - ul) * ul - right.rho * (s_r - ur) * ur) * safeInverse(denominator);

        CompressibleStarState state;
        if (0.0 < s_l)
        {
            state.rho = left.rho;
            state.p = left.p;
            state.vel = left.vel;
            state.energy = left.energy;
            return state;
        }

        if (s_l <= 0.0 && 0.0 <= s_star)
        {
            state.p = 0.5 * (left.p + right.p) +
                      0.5 * (left.rho * (s_l - ul) * (s_star - ul) + right.rho * (s_r - ur) * (s_star - ur)) * limiter;
            state.vel = left.vel - normal * (s_star - ul);
            state.rho = left.rho * (s_l - ul) * safeInverse(s_l - s_star);
            state.energy = ((s_l - ul) * left.energy - left.p * ul + state.p * s_star) * safeInverse(s_l - s_star);
            return state;
        }

        if (s_star <= 0.0 && 0.0 <= s_r)
        {
            state.p = 0.5 * (left.p + right.p) +
                      0.5 * (left.rho * (s_l - ul) * (s_star - ul) + right.rho * (s_r - ur) * (s_star - ur)) * limiter;
            state.vel = right.vel - normal * (s_star - ur);
            state.rho = right.rho * (s_r - ur) * safeInverse(s_r - s_star);
            state.energy = ((s_r - ur) * right.energy - right.p * ur + state.p * s_star) * safeInverse(s_r - s_star);
            return state;
        }

        state.rho = right.rho;
        state.p = right.p;
        state.vel = right.vel;
        state.energy = right.energy;
        return state;
    }

  private:
    IdealGasEquationOfState eos_;
    Scalar limiter_parameter_;
};
} // namespace math
} // namespace eulerian_sph
