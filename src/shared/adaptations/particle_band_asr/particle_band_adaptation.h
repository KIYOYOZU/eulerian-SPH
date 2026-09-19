/**
 * @file 	particle_band_adaptation.h
 * @brief 	Particle-band adaptive spatial resolution (SPH-ASR).
 * @details Implements the band-wise resolution adaptation of
 * 			Yang, Kong & Liu, Phys. Rev. E 104, 055308 (2021), Sec. III.A:
 * 			consecutive particle bands parallel to the material interface with
 * 			geometric spacing refinement ds_{k+1} = ds_k / C_r (C_r = 2^{1/d})
 * 			and band width dS_k = band_width_factor * ds_k (Eq. 38).
 * 			The reference smoothing length is h_ref = h_spacing_ratio * ds_max,
 * 			so that h_ratio = h_ref / h_local >= 1 always holds (required by
 * 			the multi-level cell linked list).
 * @author 	KIYOYOZU
 */

#ifndef PARTICLE_BAND_ADAPTATION_H
#define PARTICLE_BAND_ADAPTATION_H

#include "adaptation.h"

namespace SPH
{
/**
 * @class ParticleBandAdaptation
 * @brief Adaptation policy for the particle-band ASR method.
 * The first constructor argument (global_resolution) is supplied automatically
 * by SPHBody::defineAdaptation and equals the finest particle spacing ds_min.
 */
class ParticleBandAdaptation : public AdaptiveSmoothingLength
{
  public:
    typedef ParticleBandAdaptation CellLinkedListIdentifier;

    /**
     * @param global_resolution finest particle spacing ds_min (injected)
     * @param coarse_to_fine_ratio ds_max / ds_min (1.0 for uniform mode)
     * @param band_coefficient geometric refinement ratio C_r (Eq. 34, 2^{1/d})
     * @param band_width_factor band width in units of local spacing (Eq. 38)
     * @param h_spacing_ratio h_ref / ds_max (paper: 1.5 via h_r = 1.5 V^{1/d})
     */
    ParticleBandAdaptation(Real global_resolution, Real coarse_to_fine_ratio,
                           Real band_coefficient = 1.4142135623730951, Real band_width_factor = 5.0,
                           Real h_spacing_ratio = 1.5);
    virtual ~ParticleBandAdaptation() {};

    virtual void initializeAdaptationVariables(BaseParticles &base_particles) override;
    /** Not used by the banded generator; kept for interface completeness. */
    virtual Real getLocalSpacing(Shape &shape, const Vecd &position) override { return spacing_ref_; };

    int BandCount() const { return band_count_; };
    /** ds_k of band k; band 0 is the finest band at the interface. */
    Real BandSpacing(int k) const { return band_spacing_[k]; };
    /** Cumulative band boundary distance from the interface, S_0 = 0. */
    Real BandBoundary(int k) const { return band_boundary_[k]; };
    Real BandCoefficient() const { return band_coefficient_; };
    Real CoarsestSpacing() const { return band_spacing_[band_count_]; };
    Real FinestSpacing() const { return band_spacing_[0]; };
    /** Band index owning the distance to the interface (clamped). */
    int BandOfDistance(Real distance) const;
    /** Band index whose spacing is geometrically closest to the given one. */
    int BandOfSpacing(Real spacing) const;
    /** Reference neighbor number N_r of Eq. (6), lattice calibrated. */
    Real ReferenceNeighborNumber() const { return n_r_; };
    void SetReferenceNeighborNumber(Real n_r) { n_r_ = n_r; };

    DiscreteVariable<int> *dvParticleBand() { return dv_band_; };
    DiscreteVariable<Real> *dvReferenceSpacing() { return dv_ref_spacing_; };

  protected:
    Real ds_min_, band_coefficient_, band_width_factor_;
    int band_count_;                    /**< highest band index (bands 0..band_count_) */
    StdVec<Real> band_spacing_;         /**< ds_k, k = 0 (finest, interface) .. band_count_ (coarsest) */
    StdVec<Real> band_boundary_;        /**< cumulative S_k from interface, size band_count_+2 */
    Real n_r_;                          /**< reference neighbor number, Eq. (6) */

    DiscreteVariable<int> *dv_band_;
    DiscreteVariable<Real> *dv_ref_spacing_;
    int *band_;
    Real *ref_spacing_;

    Real computeLatticeNeighborNumber() const;
};
} // namespace SPH
#endif // PARTICLE_BAND_ADAPTATION_H
