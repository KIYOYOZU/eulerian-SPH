#pragma once

#include "../core/types.h"

namespace eulerian_sph
{
namespace kernel
{
class KernelBase
{
  public:
    virtual ~KernelBase() = default;

    virtual Scalar smoothingLength() const = 0;
    virtual Scalar supportRadius() const = 0;
    virtual Scalar weight(Scalar distance) const = 0;
    virtual Scalar dWeightDr(Scalar distance) const = 0;
};
} // namespace kernel
} // namespace eulerian_sph
