#pragma once

#include "../core/types.h"

#include <cmath>

namespace eulerian_sph
{
namespace math
{
inline Scalar minmod(const Scalar a, const Scalar b)
{
    if (a * b <= 0.0)
    {
        return 0.0;
    }
    return (std::abs(a) < std::abs(b)) ? a : b;
}
} // namespace math
} // namespace eulerian_sph
