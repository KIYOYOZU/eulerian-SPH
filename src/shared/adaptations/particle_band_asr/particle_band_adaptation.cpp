/**
 * @file 	particle_band_adaptation.cpp
 * @brief 	Implementation of the particle-band ASR adaptation policy.
 * @author 	KIYOYOZU
 */

#include "particle_band_adaptation.h"

#include "base_particles.hpp"

#include <cmath>
#include <iostream>

namespace SPH
{
//=================================================================================================//
ParticleBandAdaptation::ParticleBandAdaptation(
    Real global_resolution, Real coarse_to_fine_ratio,
    Real band_coefficient, Real band_width_factor, Real h_spacing_ratio)
    : AdaptiveSmoothingLength(global_resolution * coarse_to_fine_ratio, h_spacing_ratio, 1.0,
                              coarse_to_fine_ratio > 1.0
                                  ? (int)std::ceil(std::log2(coarse_to_fine_ratio) - 1e-10)
                                  : 0),
      ds_min_(global_resolution), band_coefficient_(band_coefficient),
      band_width_factor_(band_width_factor), dv_band_(nullptr), dv_ref_spacing_(nullptr),
      band_(nullptr), ref_spacing_(nullptr)
{
    if (coarse_to_fine_ratio < 1.0)
    {
        std::cout << "\n Error: ParticleBandAdaptation coarse_to_fine_ratio < 1!" << std::endl;
        exit(1);
    }
    // band 0 (finest, ds_min) at the interface, geometric coarsening by C_r
    // (Eq. 34), capped at ds_max; band width dS_k = band_width_factor * ds_k (Eq. 38).
    band_count_ = coarse_to_fine_ratio > 1.0
                      ? (int)std::ceil(std::log(coarse_to_fine_ratio) / std::log(band_coefficient_) - 1e-10)
                      : 0;
    band_spacing_.resize(band_count_ + 1);
    for (int k = 0; k <= band_count_; ++k)
        band_spacing_[k] = SMIN(ds_min_ * pow(band_coefficient_, k), spacing_ref_);
    band_spacing_[0] = ds_min_;
    band_spacing_[band_count_] = spacing_ref_;

    band_boundary_.resize(band_count_ + 2, 0.0);
    for (int k = 0; k <= band_count_; ++k)
        band_boundary_[k + 1] = band_boundary_[k] + band_width_factor_ * band_spacing_[k];

    n_r_ = computeLatticeNeighborNumber();
}
//=================================================================================================//
Real ParticleBandAdaptation::computeLatticeNeighborNumber() const
{
    // neighbor count of a correctly spaced particle on a square lattice of
    // unit spacing within the kernel support of h = h_spacing_ratio_ (r < 2h)
    Real cutoff = kernel_ptr_->KernelSize() * h_spacing_ratio_;
    int depth = (int)std::ceil(cutoff);
    Real count(0);
    for (int j = -depth; j <= depth; ++j)
        for (int i = -depth; i <= depth; ++i)
        {
            Real distance_sq = Real(i * i + j * j);
            if (distance_sq > 0 && distance_sq < cutoff * cutoff)
                count += 1.0;
        }
    return count;
}
//=================================================================================================//
int ParticleBandAdaptation::BandOfDistance(Real distance) const
{
    int k = 0;
    while (k < band_count_ && distance >= band_boundary_[k + 1])
        ++k;
    return k;
}
//=================================================================================================//
int ParticleBandAdaptation::BandOfSpacing(Real spacing) const
{
    // buffer particles have zero volume measure before being spawned
    if (spacing <= TinyReal)
        return band_count_;
    int k = (int)std::round(std::log(spacing / ds_min_) / std::log(band_coefficient_));
    return std::max(0, std::min(band_count_, k));
}
//=================================================================================================//
void ParticleBandAdaptation::initializeAdaptationVariables(BaseParticles &base_particles)
{
    AdaptiveSmoothingLength::initializeAdaptationVariables(base_particles);
    dv_band_ = base_particles.registerStateVariable<int>(
        "ParticleBand", [&](size_t i) -> int
        { return BandOfSpacing(base_particles.ParticleSpacing(i)); });
    band_ = dv_band_->Data();
    dv_ref_spacing_ = base_particles.registerStateVariable<Real>(
        "ReferenceSpacing", [&](size_t i) -> Real
        { return band_spacing_[band_[i]]; });
    ref_spacing_ = dv_ref_spacing_->Data();
    // Spawn/Remove copy the evolving set only: band and spacing target must
    // survive particle split/merge events.
    base_particles.addEvolvingVariable<int>("ParticleBand");
    base_particles.addEvolvingVariable<Real>("ReferenceSpacing");
}
//=================================================================================================//
} // namespace SPH
