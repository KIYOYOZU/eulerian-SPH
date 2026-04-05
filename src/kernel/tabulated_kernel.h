#pragma once

#include "kernel_base.h"

#include <algorithm>

namespace eulerian_sph
{
namespace kernel
{
class TabulatedKernel : public KernelBase
{
  public:
    TabulatedKernel(std::shared_ptr<KernelBase> base_kernel, const std::size_t table_resolution = 20)
        : base_kernel_(std::move(base_kernel)),
          support_radius_(base_kernel_ ? base_kernel_->supportRadius() : 0.0),
          smoothing_length_(base_kernel_ ? base_kernel_->smoothingLength() : 0.0),
          table_resolution_(std::max<std::size_t>(2, table_resolution)),
          kernel_size_(smoothing_length_ > kEpsilon ? support_radius_ * safeInverse(smoothing_length_) : 0.0),
          dq_(kernel_size_ * safeInverse(static_cast<Scalar>(table_resolution_))),
          delta_q_0_((-1.0 * dq_) * (-2.0 * dq_) * (-3.0 * dq_)),
          delta_q_1_(dq_ * (-1.0 * dq_) * (-2.0 * dq_)),
          delta_q_2_((2.0 * dq_) * dq_ * (-1.0 * dq_)),
          delta_q_3_((3.0 * dq_) * (2.0 * dq_) * dq_)
    {
        weight_table_.resize(table_resolution_ + 4, 0.0);
        gradient_table_.resize(table_resolution_ + 4, 0.0);
        if (!base_kernel_)
        {
            return;
        }

        for (std::size_t index = 0; index != table_resolution_ + 4; ++index)
        {
            const Scalar q = (static_cast<Scalar>(index) - 1.0) * dq_;
            const Scalar distance = q * smoothing_length_;
            weight_table_[index] = base_kernel_->weight(distance);
            gradient_table_[index] = base_kernel_->dWeightDr(distance);
        }
    }

    Scalar smoothingLength() const override { return smoothing_length_; }
    Scalar supportRadius() const override { return support_radius_; }

    Scalar weight(const Scalar distance) const override
    {
        return interpolate(weight_table_, distance);
    }

    Scalar dWeightDr(const Scalar distance) const override
    {
        return interpolate(gradient_table_, distance);
    }

  private:
    Scalar interpolate(const std::vector<Scalar> &table, const Scalar distance) const
    {
        if (!base_kernel_ || table.empty())
        {
            return table.empty() ? 0.0 : table.front();
        }
        if (distance >= support_radius_)
        {
            return 0.0;
        }

        const Scalar q = std::max(distance, 0.0) * safeInverse(std::max(smoothing_length_, kEpsilon));
        const Scalar location = std::floor(q * safeInverse(dq_));
        const std::size_t i = std::min<std::size_t>(
            static_cast<std::size_t>(std::max(location + 1.0, 1.0)),
            table_resolution_);

        const Scalar fraction_1 = q - location * dq_;
        const Scalar fraction_0 = fraction_1 + dq_;
        const Scalar fraction_2 = fraction_1 - dq_;
        const Scalar fraction_3 = fraction_1 - 2.0 * dq_;

        return (fraction_1 * fraction_2 * fraction_3) * safeInverse(delta_q_0_) * table[i - 1] +
               (fraction_0 * fraction_2 * fraction_3) * safeInverse(delta_q_1_) * table[i] +
               (fraction_0 * fraction_1 * fraction_3) * safeInverse(delta_q_2_) * table[i + 1] +
               (fraction_0 * fraction_1 * fraction_2) * safeInverse(delta_q_3_) * table[i + 2];
    }

    std::shared_ptr<KernelBase> base_kernel_;
    Scalar support_radius_;
    Scalar smoothing_length_;
    std::size_t table_resolution_;
    Scalar kernel_size_;
    Scalar dq_;
    Scalar delta_q_0_;
    Scalar delta_q_1_;
    Scalar delta_q_2_;
    Scalar delta_q_3_;
    std::vector<Scalar> weight_table_;
    std::vector<Scalar> gradient_table_;
};
} // namespace kernel
} // namespace eulerian_sph
