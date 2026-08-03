#ifndef CYLINDER_3D_COMPRESSIBLE_STATE_HPP
#define CYLINDER_3D_COMPRESSIBLE_STATE_HPP

/**
 * Seven-field conservative state contract for the Ma=0.3 compressible cylinder:
 *
 *     Density, VolumetricMeasure, Mass, Velocity, Momentum, Pressure, TotalEnergy
 *
 * with the invariants
 *     Mass        = Density * VolumetricMeasure
 *     Momentum    = Mass * Velocity
 *     TotalEnergy = [Pressure/(gamma-1) + 0.5*Density*|Velocity|^2] * VolumetricMeasure
 * and positivity of Density, Pressure and internal energy.
 *
 * setCompressiblePrimitiveState() is the single write entry point used by the
 * initial condition, the ghost boundary reset and the restart consistency gate;
 * no code path may write a subset of the seven fields directly. The weakly
 * compressible EulerianWeaklyCompressibleBoundaryState /
 * syncEulerianWeaklyCompressibleState() pair must not be used here: it does not
 * manage TotalEnergy and assumes a different EOS.
 */

#include "cylinder_3d_compressible_data.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <string>

namespace SPH
{
namespace cylinder_3d_compressible
{

inline bool isFiniteRealValue(Real value)
{
    return std::isfinite(static_cast<double>(value));
}

inline bool isFiniteVecValue(const Vecd &value)
{
    for (int axis = 0; axis != Dimensions; ++axis)
    {
        if (!isFiniteRealValue(value[axis]))
        {
            return false;
        }
    }
    return true;
}

/** Raw pointers to the seven conservative fields plus the EOS exponent. */
struct CompressibleConservativeStateView
{
    Real *rho = nullptr;
    Real *Vol = nullptr;
    Real *mass = nullptr;
    Vecd *vel = nullptr;
    Vecd *mom = nullptr;
    Real *p = nullptr;
    Real *E = nullptr;
    Real gamma = 0.0;
};

inline CompressibleConservativeStateView makeCompressibleStateView(BaseParticles &particles, Real gamma)
{
    if (gamma <= 1.0)
    {
        throw std::runtime_error("makeCompressibleStateView: gamma must be greater than 1.");
    }
    CompressibleConservativeStateView state;
    state.rho = particles.registerStateVariableData<Real>("Density");
    state.Vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    state.mass = particles.registerStateVariableData<Real>("Mass");
    state.vel = particles.registerStateVariableData<Vecd>("Velocity");
    state.mom = particles.registerStateVariableData<Vecd>("Momentum");
    state.p = particles.registerStateVariableData<Real>("Pressure");
    state.E = particles.registerStateVariableData<Real>("TotalEnergy");
    state.gamma = gamma;
    if (state.rho == nullptr || state.Vol == nullptr || state.mass == nullptr ||
        state.vel == nullptr || state.mom == nullptr || state.p == nullptr || state.E == nullptr)
    {
        throw std::runtime_error("makeCompressibleStateView: one of the seven conservative fields is missing.");
    }
    return state;
}

/** Register the six evolving flow fields so restart serializes the full state. */
inline void registerCompressibleEvolvingVariables(BaseParticles &particles)
{
    particles.addEvolvingVariable<Real>("Density");
    particles.addEvolvingVariable<Real>("Mass");
    particles.addEvolvingVariable<Real>("Pressure");
    particles.addEvolvingVariable<Real>("TotalEnergy");
    particles.addEvolvingVariable<Vecd>("Velocity");
    particles.addEvolvingVariable<Vecd>("Momentum");
}

/**
 * Velocity implied by the conservative pair (Momentum, Mass).
 *
 * This is the velocity the shared solver itself uses: the second half step
 * updates Mass/Density/TotalEnergy and recovers the pressure from
 * mom_/mass_ (eulerian_compressible_fluid_integration.cpp:176) WITHOUT writing
 * the result back into Velocity. The stored Velocity therefore lags the mass
 * update by one step, while Momentum/Mass is always current. Every invariant
 * that the solver actually maintains must be evaluated against this velocity.
 *
 * On the write path (initial condition, ghost reset, restart) Momentum is set to
 * Mass * Velocity, so the two agree exactly and the distinction is invisible.
 */
inline Vecd conservativeVelocity(const CompressibleConservativeStateView &state, size_t index)
{
    return state.mom[index] / state.mass[index];
}

/**
 * Internal energy per unit volume recovered from the conservative state.
 *
 * Uses the conservative velocity so that the recovered value matches the
 * solver's own energy split bit for bit; using the stored (lagging) Velocity
 * would report a spurious O(dMass/Mass) inconsistency.
 */
inline Real internalEnergyPerVolume(const CompressibleConservativeStateView &state, size_t index)
{
    return state.E[index] / state.Vol[index] -
           0.5 * state.rho[index] * conservativeVelocity(state, index).squaredNorm();
}

/**
 * The only write entry point for the seven-field state.
 *
 * Fail-fast on a non-physical primitive state instead of clamping: a silently
 * repaired ghost or restart state would hide the real defect.
 */
inline void setCompressiblePrimitiveState(CompressibleConservativeStateView &state, size_t index,
                                          Real rho, const Vecd &vel, Real p)
{
    if (!isFiniteRealValue(rho) || !isFiniteRealValue(p) || !isFiniteVecValue(vel))
    {
        throw std::runtime_error("setCompressiblePrimitiveState: non-finite primitive state.");
    }
    if (rho <= 0.0)
    {
        throw std::runtime_error("setCompressiblePrimitiveState: density must be positive.");
    }
    if (p <= 0.0)
    {
        throw std::runtime_error("setCompressiblePrimitiveState: pressure must be positive.");
    }
    const Real volume = state.Vol[index];
    if (!isFiniteRealValue(volume) || volume <= 0.0)
    {
        throw std::runtime_error("setCompressiblePrimitiveState: volumetric measure must be positive.");
    }

    state.rho[index] = rho;
    state.p[index] = p;
    state.vel[index] = vel;
    state.mass[index] = rho * volume;
    state.mom[index] = state.mass[index] * vel;
    state.E[index] = (p / (state.gamma - 1.0) + 0.5 * rho * vel.squaredNorm()) * volume;
}

/** Apply the freestream primitive state to one particle. */
inline void setFreestreamState(CompressibleConservativeStateView &state, size_t index,
                               const CompressibleFreestreamState &freestream)
{
    setCompressiblePrimitiveState(state, index, freestream.rho, freestream.vel, freestream.p);
}

/** Pressure recovered from the conservative state by the ideal-gas EOS. */
inline Real recoveredPressure(const CompressibleConservativeStateView &state, size_t index)
{
    return (state.gamma - 1.0) * internalEnergyPerVolume(state, index);
}

/** Local sound speed sqrt(gamma*p/rho) from the stored pressure. */
inline Real localSoundSpeed(const CompressibleConservativeStateView &state, size_t index)
{
    return std::sqrt(state.gamma * state.p[index] / state.rho[index]);
}

inline Real localMachNumber(const CompressibleConservativeStateView &state, size_t index)
{
    return state.vel[index].norm() / localSoundSpeed(state, index);
}

inline bool isThermodynamicallyAdmissible(const CompressibleConservativeStateView &state, size_t index)
{
    if (!isFiniteRealValue(state.rho[index]) || !isFiniteRealValue(state.Vol[index]) ||
        !isFiniteRealValue(state.mass[index]) || !isFiniteVecValue(state.vel[index]) ||
        !isFiniteVecValue(state.mom[index]) || !isFiniteRealValue(state.p[index]) ||
        !isFiniteRealValue(state.E[index]))
    {
        return false;
    }
    if (state.rho[index] <= 0.0 || state.Vol[index] <= 0.0 || state.mass[index] <= 0.0 ||
        state.p[index] <= 0.0)
    {
        return false;
    }
    return internalEnergyPerVolume(state, index) > 0.0;
}

struct StateConsistencyResult
{
    size_t total = 0;
    size_t first_bad = 0;
    Real max_mass_error = 0.0;
    Real max_momentum_error = 0.0;
    Real max_energy_error = 0.0;
    Real max_pressure_error = 0.0;
    /**
     * Relative gap between the stored Velocity and Momentum/Mass. Reported, not
     * judged: it equals |dMass|*dt/Mass by construction of the shared two-half
     * update and is O(1e-2) on a legitimate run. Only the write path (where the
     * two are set together) may require it to vanish.
     */
    Real max_velocity_lag = 0.0;
    Real rho_min = std::numeric_limits<Real>::max();
    Real p_min = std::numeric_limits<Real>::max();
    Real internal_energy_min = std::numeric_limits<Real>::max();
    Real mach_max = 0.0;
    bool pass = false;
};

/**
 * Verify the seven-field invariants over [begin, end).
 *
 * Errors are relative to the local scale so one tolerance covers both the
 * freestream and the compressed states near the cylinder.
 *
 * TWO DISTINCT CONTRACTS, deliberately not merged:
 *
 *  - The SOLVER contract, always judged. Mass == Density*VolumetricMeasure,
 *    and TotalEnergy / Pressure consistent with the conservative velocity
 *    Momentum/Mass. These the shared scheme maintains exactly at every step.
 *
 *  - The WRITE-PATH contract, judged only when require_velocity_sync is set.
 *    Momentum == Mass*Velocity. This holds after setCompressiblePrimitiveState()
 *    but is structurally broken mid-scheme: the second half step advances Mass
 *    without rewriting Velocity, so the gap is exactly |dMass|*dt/Mass. Judging
 *    it on an evolved state would fail a correct run; ignoring it on the write
 *    path would miss a genuine partial write. Hence the flag rather than a
 *    loosened tolerance -- a loosened tolerance would silently accept a real
 *    partial write on the initial condition too.
 */
inline StateConsistencyResult checkStateConsistency(const CompressibleConservativeStateView &state,
                                                    size_t begin, size_t end,
                                                    const std::string &stage,
                                                    Real velocity_scale,
                                                    Real tolerance = 1.0e-12,
                                                    bool require_velocity_sync = true)
{
    StateConsistencyResult result;
    result.total = end > begin ? end - begin : 0;
    result.first_bad = end;

    for (size_t i = begin; i != end; ++i)
    {
        if (!isThermodynamicallyAdmissible(state, i))
        {
            if (result.first_bad == end)
            {
                result.first_bad = i;
            }
            continue;
        }

        const Real volume = state.Vol[i];
        const Real mass_reference = state.rho[i] * volume;
        const Vecd vel_conservative = conservativeVelocity(state, i);
        // Energy and pressure are checked against the conservative velocity,
        // which is the one the solver's own EOS inversion uses.
        const Real energy_reference =
            (state.p[i] / (state.gamma - 1.0) +
             0.5 * state.rho[i] * vel_conservative.squaredNorm()) *
            volume;
        const Vecd momentum_reference = state.mass[i] * state.vel[i];

        result.max_mass_error = SMAX(
            result.max_mass_error, std::fabs(state.mass[i] - mass_reference) / mass_reference);
        // The scale is floored at velocity_scale (the freestream speed): a
        // stagnation-point particle approaches |u| = 0, and dividing by its own
        // speed would blow the relative error up on a perfectly consistent state.
        const Real velocity_reference_scale =
            mass_reference * SMAX(vel_conservative.norm(), velocity_scale);
        const Real velocity_lag =
            (state.mom[i] - momentum_reference).norm() / velocity_reference_scale;
        result.max_velocity_lag = SMAX(result.max_velocity_lag, velocity_lag);
        if (require_velocity_sync)
        {
            result.max_momentum_error = SMAX(result.max_momentum_error, velocity_lag);
        }
        result.max_energy_error = SMAX(
            result.max_energy_error, std::fabs(state.E[i] - energy_reference) / energy_reference);
        result.max_pressure_error = SMAX(
            result.max_pressure_error, std::fabs(recoveredPressure(state, i) - state.p[i]) / state.p[i]);

        result.rho_min = SMIN(result.rho_min, state.rho[i]);
        result.p_min = SMIN(result.p_min, state.p[i]);
        result.internal_energy_min = SMIN(result.internal_energy_min, internalEnergyPerVolume(state, i));
        result.mach_max = SMAX(result.mach_max, localMachNumber(state, i));
    }

    result.pass = result.first_bad == end &&
                  result.max_mass_error <= tolerance &&
                  result.max_momentum_error <= tolerance &&
                  result.max_energy_error <= tolerance &&
                  result.max_pressure_error <= tolerance;

    std::cout << "[Cylinder3DCompressible][StateGate] " << stage
              << " count=" << result.total
              << " mass_err=" << result.max_mass_error
              << " mom_err=" << result.max_momentum_error
              << " E_err=" << result.max_energy_error
              << " p_err=" << result.max_pressure_error
              << " vel_lag=" << result.max_velocity_lag
              << (require_velocity_sync ? " (judged)" : " (reported)")
              << " rho_min=" << result.rho_min
              << " p_min=" << result.p_min
              << " e_min=" << result.internal_energy_min
              << " Ma_max=" << result.mach_max
              << " pass=" << (result.pass ? "yes" : "NO") << std::endl;
    if (result.first_bad != end)
    {
        std::cout << "[Cylinder3DCompressible][BadState] i=" << result.first_bad
                  << " rho=" << state.rho[result.first_bad]
                  << " Vol=" << state.Vol[result.first_bad]
                  << " mass=" << state.mass[result.first_bad]
                  << " p=" << state.p[result.first_bad]
                  << " E=" << state.E[result.first_bad]
                  << " |u|=" << state.vel[result.first_bad].norm() << std::endl;
    }
    return result;
}

/**
 * Total mass and total energy over [begin, end), plus the time-integrated change
 * rates the solver applied.
 *
 * The domain is NOT closed: four x/y faces are open, so
 *     d/dt integral(rho) = -closed_surface_integral(rho u . n)
 * and the totals are free to change. "sum(Mass) stays constant" is therefore the
 * wrong conservation statement here -- it only holds for a closed domain.
 *
 * accumulated_* are NOT boundary fluxes, despite being the natural place to look
 * for one:
 *   - they re-sum the very arrays the solver used to advance Mass / TotalEnergy,
 *     so their difference from the state change is an identity, not a test
 *     (see reportConservationDrift);
 *   - dE_dt_ is additionally seeded with volumetric work, not flux:
 *         energy_change_rate = force_prior_[i].dot(vel_[i])
 *     (eulerian_compressible_fluid_integration.cpp:123), and with the viscous
 *     force registered as a ForcePrior that term carries the viscous power of
 *     every interior particle.
 * They are kept only so the printed report can state explicitly how much of the
 * change is accounted for, and are never gated.
 */
struct ConservationBudget
{
    Real total_mass = 0.0;
    Real total_energy = 0.0;
    Real accumulated_mass_rate = 0.0;
    Real accumulated_energy_rate = 0.0;
};

inline ConservationBudget accumulateConservationBudget(const CompressibleConservativeStateView &state,
                                                       size_t begin, size_t end)
{
    ConservationBudget budget;
    for (size_t i = begin; i != end; ++i)
    {
        budget.total_mass += state.mass[i];
        budget.total_energy += state.E[i];
    }
    return budget;
}

/**
 * One step's worth of the applied change rates, added into a running total.
 *
 * Call after the second half step has written MassChangeRate /
 * TotalEnergyChangeRate and before the next step overwrites them, with the same
 * dt that advanced the state. Reporting only -- see ConservationBudget for why
 * this is not a conservation measurement.
 */
inline void accumulateAppliedRates(ConservationBudget &budget,
                                   const Real *dmass_dt, const Real *dE_dt,
                                   size_t begin, size_t end, Real dt)
{
    Real mass_rate = 0.0;
    Real energy_rate = 0.0;
    for (size_t i = begin; i != end; ++i)
    {
        mass_rate += dmass_dt[i];
        energy_rate += dE_dt[i];
    }
    budget.accumulated_mass_rate += mass_rate * dt;
    budget.accumulated_energy_rate += energy_rate * dt;
}

} // namespace cylinder_3d_compressible
} // namespace SPH

#endif // CYLINDER_3D_COMPRESSIBLE_STATE_HPP
