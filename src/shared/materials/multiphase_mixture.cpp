#include "multiphase_mixture.h"

namespace SPH
{
//=============================================================================================//
MultiphaseMixture::MultiphaseMixture(StiffenedGas &phase_1, StiffenedGas &phase_2)
    : phase_1_(phase_1), phase_2_(phase_2) {}
//=============================================================================================//
Real MultiphaseMixture::MixtureA(Real alpha1) const
{
    Real alpha2 = 1.0 - alpha1;
    Real g1m1 = phase_1_.HeatCapacityRatio() - 1.0;
    Real g2m1 = phase_2_.HeatCapacityRatio() - 1.0;
    return alpha1 / g1m1 + alpha2 / g2m1;
}
//=============================================================================================//
Real MultiphaseMixture::MixtureB(Real alpha1) const
{
    Real alpha2 = 1.0 - alpha1;
    Real g1 = phase_1_.HeatCapacityRatio();
    Real g2 = phase_2_.HeatCapacityRatio();
    Real g1m1 = g1 - 1.0;
    Real g2m1 = g2 - 1.0;
    return alpha1 * g1 * phase_1_.ReferencePressure() / g1m1
         + alpha2 * g2 * phase_2_.ReferencePressure() / g2m1;
}
//=============================================================================================//
Real MultiphaseMixture::MixtureGamma(Real alpha1) const
{
    Real A = MixtureA(alpha1);
    return 1.0 + 1.0 / A;
}
//=============================================================================================//
Real MultiphaseMixture::MixturePInf(Real alpha1) const
{
    Real A = MixtureA(alpha1);
    Real B = MixtureB(alpha1);
    Real gamma_mix = 1.0 + 1.0 / A;
    // B = gamma_mix * p_inf_mix / (gamma_mix - 1) = gamma_mix * p_inf_mix * A
    // => p_inf_mix = B / (gamma_mix * A)
    return B / (gamma_mix * A);
}
//=============================================================================================//
Real MultiphaseMixture::MixturePressure(Real alpha1, Real rho_e) const
{
    Real A = MixtureA(alpha1);
    Real B = MixtureB(alpha1);
    return (rho_e - B) / A;
}
//=============================================================================================//
Real MultiphaseMixture::MixtureInternalEnergyPerVolume(Real alpha1, Real p) const
{
    Real A = MixtureA(alpha1);
    Real B = MixtureB(alpha1);
    return A * p + B;
}
//=============================================================================================//
Real MultiphaseMixture::MixtureSoundSpeed(Real alpha1, Real rho, Real p) const
{
    Real alpha2 = 1.0 - alpha1;
    Real g1 = phase_1_.HeatCapacityRatio();
    Real g2 = phase_2_.HeatCapacityRatio();
    Real p_inf1 = phase_1_.ReferencePressure();
    Real p_inf2 = phase_2_.ReferencePressure();

    // Wood sound speed: 1/(rho*c^2) = alpha1/(gamma1*(p+p_inf1)) + alpha2/(gamma2*(p+p_inf2))
    // Floor (p + p_inf) and rho so a transient negative-pressure undershoot keeps
    // the sound speed real instead of producing NaN.
    Real pp1 = SMAX(p + p_inf1, TinyReal);
    Real pp2 = SMAX(p + p_inf2, TinyReal);
    Real inv_rho_c2 = alpha1 / (g1 * pp1) + alpha2 / (g2 * pp2);
    return std::sqrt(1.0 / SMAX(rho * inv_rho_c2, TinyReal));
}
//=============================================================================================//
} // namespace SPH
