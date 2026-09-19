/**
 * @file 	update_particle_bands.h
 * @brief 	Interface-following band assignment of the particle-band ASR
 * 			method, Sec. III.A of Yang, Kong & Liu, PRE 104, 055308 (2021).
 * @author 	KIYOYOZU
 */

#ifndef UPDATE_PARTICLE_BANDS_H
#define UPDATE_PARTICLE_BANDS_H

#include "particle_band_adaptation.h"

#include "base_body_relation.h"

namespace SPH
{
/**
 * @class UpdateParticleBands
 * @brief Reassigns the band index of every particle from its distance to the
 * material interface, so that band 0 (finest) follows the interface as it
 * translates. Band boundaries are the cumulative widths S_k = sum_m dS_m with
 * dS_k = band_width_factor * ds_k (Eq. 38); the distance is clamped to the
 * coarsest band.
 *
 * Interface detection has two paths:
 * - planar fast path: column-wise sign change of the equal-volume phase
 *   indicator sum(alpha V) - sum((1 - alpha) V) locates the alpha = 0.5
 *   contour of a single x = const interface to sub-bin accuracy by linear
 *   interpolation (a (2 alpha - 1) V weighting would balance at alpha ~ 0.92
 *   for a 1.4/1000 density pair and stall on smeared profiles);
 * - general fallback: multi-source Dijkstra over the neighbor graph from
 *   interface seeds (alpha within tol of 0.5, or a neighbor alpha jump
 *   above 0.5), edge weights are the cached pair distances.
 *
 * A band change only takes effect when the distance has crossed the relevant
 * band boundary by at least band_hysteresis * ds_target, so that particles
 * jittering around a boundary do not flip their spacing target every update
 * (which would drive a split/merge oscillation through gamma = Vol/ds^2).
 *
 * Optional shock tracking (methodological extension beyond the paper, which
 * refines around the material interface only): with shock_band on, shock seed
 * particles are marked at PARTICLE level by three orthogonal gates, each
 * excluding one class of non-shock feature -- same-pure-phase neighbor pairs
 * (kills the material interface, across which pressure is continuous),
 * negative kernel-gradient velocity divergence (kills rarefactions; the
 * difference form is exactly zero on uniform velocity), and a pair pressure
 * jump above shock_rel_jump of the local level (kills pure-noise fields).
 * Band assignment then runs as ONE multi-source Dijkstra from the union of
 * interface and shock seeds, so every particle is banded by its graph
 * distance to the nearest feature: geometry-free, 3D-ready, multiple shocks
 * admitted, and no argmax competition in which one feature could starve
 * another. Ringing particles can pass the gates; that only refines a noisy
 * region and is benign by the same no-competition property. DetectedShock()
 * reports the seed centroid x for monitoring. With shock_band off the planar
 * fast path / Dijkstra fallback behave exactly as before.
 *
 * Serial execution: the Dijkstra path and the writes to the band field are
 * single-threaded; the planar path scans particles twice.
 */
class UpdateParticleBands
{
  public:
    UpdateParticleBands(RealBody &real_body, BaseInnerRelation &inner_relation,
                        Real interface_alpha_tol = 0.1, Real band_hysteresis = 0.5,
                        bool force_dijkstra = false, bool shock_band = false,
                        Real shock_rel_jump = 0.01);
    /** @return number of particles whose band changed. */
    size_t exec();

    /** x position of the planar interface found by the last exec (NaN if none). */
    Real DetectedInterface() const { return detected_interface_; };
    /** x position of the planar shock found by the last exec (NaN if none / shock_band off). */
    Real DetectedShock() const { return detected_shock_; };

  protected:
    ParticleBandAdaptation &adaptation_;
    BaseParticles &particles_;
    ParticleConfiguration &inner_configuration_;
    Real alpha_tol_, band_hysteresis_;
    bool force_dijkstra_;
    bool shock_band_;
    Real shock_rel_jump_;
    bool warned_no_seeds_ = false;
    Real detected_interface_, detected_shock_;
    Real *alpha_, *Vol_, *ref_spacing_, *p_;
    int *band_;
    Vecd *pos_, *vel_;
    StdVec<char> shock_seed_; /**< per-particle shock seed flags, rebuilt each exec */

    Real detectPlanarInterface(size_t total_real_particles) const;
    /** marks shock seed particles, returns their count, sets detected_shock_ centroid. */
    size_t markShockSeeds();
    size_t assignByDijkstra(size_t total_real_particles);
    /** hysteresis-guarded band update of one particle; returns true on change. */
    bool assignBand(size_t index_i, Real distance);
};

} // namespace SPH

#endif // UPDATE_PARTICLE_BANDS_H
