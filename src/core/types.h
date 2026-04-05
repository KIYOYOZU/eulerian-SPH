#pragma once

#include <array>
#include <Eigen/Dense>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace eulerian_sph
{
using Scalar = double;
using Vec2 = Eigen::Vector2d;
using Mat2 = Eigen::Matrix2d;

constexpr Scalar kEpsilon = 1.0e-12;
constexpr Scalar kPi = 3.1415926535897932384626433832795;

inline Scalar safeInverse(const Scalar value)
{
    if (std::abs(value) < kEpsilon)
    {
        return Scalar(1) / (value < 0.0 ? -kEpsilon : kEpsilon);
    }
    return Scalar(1) / value;
}
} // namespace eulerian_sph
