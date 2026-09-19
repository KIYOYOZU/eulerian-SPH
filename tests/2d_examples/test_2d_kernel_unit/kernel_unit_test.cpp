/**
 * @file 	kernel_unit_test.cpp
 * @brief 	Unit tests for kernel shape functions, derivatives, normalization
 * 			and variable-smoothing-length (h_ratio) overloads.
 */
#include "kernel_cubic_B_spline.h"
#include "kernel_hyperbolic.h"
#include "kernel_quadratic.h"
#include "kernel_wendland_c2.h"

#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

using namespace SPH;

static int failures = 0;

static void check(bool ok, const std::string &msg, Real got, Real expect, Real tol)
{
    if (!ok)
    {
        ++failures;
        printf("FAIL: %s (got %.12g, expect %.12g, tol %.2g)\n", msg.c_str(), got, expect, tol);
    }
    else
    {
        printf("PASS: %s\n", msg.c_str());
    }
}

static void checkNear(const std::string &msg, Real got, Real expect, Real tol)
{
    check(std::abs(got - expect) <= tol * std::max(1.0, std::abs(expect)), msg, got, expect, tol);
}

/** Central finite difference of f at q with step eps. */
static Real fd(const std::function<Real(Real)> &f, Real q, Real eps = 1e-6)
{
    return (f(q + eps) - f(q - eps)) / (2.0 * eps);
}

/** Simpson integration of f on [a,b] with n (even) intervals. */
static Real integrate(const std::function<Real(Real)> &f, Real a, Real b, int n = 20000)
{
    Real h = (b - a) / n;
    Real sum = f(a) + f(b);
    for (int i = 1; i < n; ++i)
        sum += (i % 2 == 0 ? 2.0 : 4.0) * f(a + i * h);
    return sum * h / 3.0;
}

static void testKernel(Kernel &kernel, const std::string &name, bool has_branch_at_one)
{
    printf("\n===== %s =====\n", name.c_str());
    const Real h = kernel.SmoothingLength();

    // A. dW_*D matches finite difference of W_*D (shape level).
    const Real q_probes = 1e-6; // dummy to keep formatting simple
    (void)q_probes;
    for (Real q : {0.25, 0.5, 0.75, 0.9, 1.1, 1.25, 1.5, 1.75, 1.95})
    {
        checkNear(name + " dW_1D FD at q=" + std::to_string(q),
                  kernel.dW_1D(q), fd([&](Real x) { return kernel.W_1D(x); }, q), 1e-5);
        checkNear(name + " d2W_1D FD at q=" + std::to_string(q),
                  kernel.d2W_1D(q), fd([&](Real x) { return kernel.dW_1D(x); }, q), 1e-4);
    }

    // B. Branch continuity at q=1 (kernels with piecewise definition).
    if (has_branch_at_one)
    {
        Real eps = 1e-10;
        checkNear(name + " W continuous at q=1", kernel.W_1D(1.0 - eps), kernel.W_1D(1.0 + eps), 1e-6);
        checkNear(name + " dW continuous at q=1", kernel.dW_1D(1.0 - eps), kernel.dW_1D(1.0 + eps), 1e-6);
    }

    // C. Normalization in 1D/2D/3D (shape integrals against the stored factors).
    Real norm_1d = kernel.FactorW1D() * h * 2.0 *
                   integrate([&](Real q) { return kernel.W_1D(q); }, 0.0, 2.0);
    checkNear(name + " 1D normalization", norm_1d, 1.0, 1e-10);

    Real norm_2d = kernel.FactorW2D() * h * h * 2.0 * Pi *
                   integrate([&](Real q) { return kernel.W_2D(q) * q; }, 0.0, 2.0);
    checkNear(name + " 2D normalization", norm_2d, 1.0, 1e-10);

    Real norm_3d = kernel.FactorW3D() * h * h * h * 4.0 * Pi *
                   integrate([&](Real q) { return kernel.W_3D(q) * q * q; }, 0.0, 2.0);
    checkNear(name + " 3D normalization", norm_3d, 1.0, 1e-10);

    // D. h_ratio overloads consistent with shape functions
    //    (q = r * inv_h * h_ratio; W scaled by h_ratio^d, dW by h_ratio^(d+1)).
    for (Real h_ratio : {1.0, 0.5, 2.0, 0.3})
    {
        Real r = 0.7 * h; // probe distance in physical units
        Real q = r / h * h_ratio;
        Vec2d disp(r, 0.0);
        Real d = disp.norm();

        checkNear(name + " W(2D,h_ratio=" + std::to_string(h_ratio) + ")",
                  kernel.W(h_ratio, d, disp),
                  kernel.FactorW2D() * kernel.W_2D(q) * h_ratio * h_ratio, 1e-12);
        checkNear(name + " dW(2D,h_ratio=" + std::to_string(h_ratio) + ")",
                  kernel.dW(h_ratio, d, disp),
                  kernel.FactorW2D() / h * kernel.dW_2D(q) * h_ratio * h_ratio * h_ratio, 1e-12);
        checkNear(name + " W0(2D,h_ratio=" + std::to_string(h_ratio) + ")",
                  kernel.W0(h_ratio, disp),
                  kernel.FactorW2D() * h_ratio * h_ratio, 1e-12);
    }

    // E. dW via dimensional FD of W (end-to-end, catches factor_dW wiring errors).
    {
        Real r = 0.6 * h;
        Real eps = 1e-7 * h;
        Vec2d disp(r, 0.0);
        Real fd_dw = (kernel.W(1.0, r + eps, Vec2d(r + eps, 0.0)) -
                      kernel.W(1.0, r - eps, Vec2d(r - eps, 0.0))) /
                     (2.0 * eps);
        checkNear(name + " dW(2D) matches dimensional FD of W", kernel.dW(1.0, r, disp), fd_dw, 1e-5);
    }
}

int main()
{
    Real h = 0.05;
    KernelHyperbolic hyperbolic(h);
    KernelQuadratic quadratic(h);
    KernelWendlandC2 wendland(h);
    KernelCubicBSpline cubic(h);

    testKernel(hyperbolic, "KernelHyperbolic", true);
    testKernel(quadratic, "KernelQuadratic", false);
    testKernel(wendland, "KernelWendlandC2", false);
    testKernel(cubic, "KernelCubicBSpline", true);

    printf("\n%s (%d failures)\n", failures == 0 ? "ALL KERNEL TESTS PASSED" : "KERNEL TESTS FAILED", failures);
    return failures == 0 ? 0 : 1;
}
