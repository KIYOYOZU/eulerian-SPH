#pragma once

#include "kernel_base.h"

#include <cmath>

namespace eulerian_sph
{
namespace kernel
{
class LaguerreGaussKernel : public KernelBase
{
  public:
    explicit LaguerreGaussKernel(const Scalar smoothing_length)
        : h_(std::max(smoothing_length, kEpsilon)),
          inv_h_(safeInverse(h_)),
          factor_w_2d_(3.0 * inv_h_ * inv_h_ / kPi),
          factor_dw_2d_(inv_h_ * factor_w_2d_),
          support_radius_(2.0 * h_) {}

    Scalar smoothingLength() const override { return h_; }
    Scalar supportRadius() const override { return support_radius_; }

    Scalar weight(const Scalar distance) const override
    {
        const Scalar q = distance * inv_h_;
        if (q >= 2.0)
        {
            return 0.0;
        }

        const Scalar polynomial = 1.0 - q * q + std::pow(q, 4) / 6.0;
        return factor_w_2d_ * polynomial * std::exp(-(q * q));
    }

    Scalar dWeightDr(const Scalar distance) const override
    {
        const Scalar q = distance * inv_h_;
        if (q >= 2.0)
        {
            return 0.0;
        }

        const Scalar polynomial = -std::pow(q, 5) / 3.0 + 8.0 * std::pow(q, 3) / 3.0 - 4.0 * q;
        return factor_dw_2d_ * polynomial * std::exp(-(q * q));
    }

  private:
    Scalar h_;
    Scalar inv_h_;
    Scalar factor_w_2d_;
    Scalar factor_dw_2d_;
    Scalar support_radius_;
};
} // namespace kernel
} // namespace eulerian_sph
