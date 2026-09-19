#include "kernel_quadratic.h"

#include <cmath>

namespace SPH
{
//=================================================================================================//
KernelQuadratic::KernelQuadratic(Real h)
    : Kernel(h, 2.0, 2.0, "QuadraticKernel")
{
    factor_W_1D_ = 0.8 * inv_h_;
    factor_W_2D_ = 1.6 * inv_h_ * inv_h_ / Pi;
    factor_W_3D_ = inv_h_ * inv_h_ * inv_h_ / Pi;
    setDerivativeParameters();
}
//=================================================================================================//
Real KernelQuadratic::W_1D(const Real q) const
{
    return 5.0 * (3.0 * q * q - 12.0 * q + 12.0) / 64.0;
}
//=================================================================================================//
Real KernelQuadratic::W_2D(const Real q) const
{
    return W_1D(q);
}
//=================================================================================================//
Real KernelQuadratic::W_3D(const Real q) const
{
    return W_1D(q);
}
//=================================================================================================//
Real KernelQuadratic::dW_1D(const Real q) const
{
    return 15.0 * (q - 2.0) / 32.0;
}
//=================================================================================================//
Real KernelQuadratic::dW_2D(const Real q) const
{
    return dW_1D(q);
}
//=================================================================================================//
Real KernelQuadratic::dW_3D(const Real q) const
{
    return dW_1D(q);
}
//=================================================================================================//
Real KernelQuadratic::d2W_1D(const Real q) const
{
    return 15.0 / 32.0;
}
//=================================================================================================//
Real KernelQuadratic::d2W_2D(const Real q) const
{
    return d2W_1D(q);
}
//=================================================================================================//
Real KernelQuadratic::d2W_3D(const Real q) const
{
    return d2W_1D(q);
}
//=================================================================================================//
} // namespace SPH
