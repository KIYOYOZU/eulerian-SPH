#include "stiffened_gas.h"

namespace SPH
{
//=============================================================================================//
StiffenedGas::StiffenedGas(Real gamma, Real p_inf) : Fluid(), gamma_(gamma), p_inf_(p_inf)
{
    material_type_name_ = "StiffenedGas";
}
//=============================================================================================//
StiffenedGas::~StiffenedGas() = default;
//=============================================================================================//
Real StiffenedGas::getPressure(Real rho, Real rho_e)
{
    return (gamma_ - 1.0) * rho_e - gamma_ * p_inf_;
}
//=============================================================================================//
Real StiffenedGas::getSoundSpeed(Real p, Real rho)
{
    // Guard against a transient pressure undershoot that would make (p + p_inf)
    // non-positive: keep the sound speed real so a single bad particle cannot
    // poison the whole field with NaN.
    return std::sqrt(gamma_ * SMAX(p + p_inf_, TinyReal) / SMAX(rho, TinyReal));
}
//=============================================================================================//
Real StiffenedGas::InternalEnergyPerVolume(Real rho, Real p) const
{
    return (p + gamma_ * p_inf_) / (gamma_ - 1.0);
}
//=============================================================================================//
} // namespace SPH
