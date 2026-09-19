/**
 * @file 	particle_split_merge.h
 * @brief 	Band-based particle splitting and merging (SPH-ASR).
 * @details Yang, Kong & Liu, Phys. Rev. E 104, 055308 (2021), Sec. III.B:
 * 			split criterion gamma = Vol_i / ds_band^2 > gamma_s (Eq. 41),
 * 			offset direction e perpendicular to the nearest-neighbor line,
 * 			daughter positions r_i +/- (lambda/2) Vol_i^{1/d} e (Eqs. 45-47);
 * 			merge criterion gamma < gamma_m (Eq. 50) with mutual nearest
 * 			pairing and conservative state combination (Eqs. 51-53).
 * 			Eulerian modification: gamma uses the volume ratio only (no mass
 * 			ratio), because shock compression would false-trigger the paper's
 * 			compressible mass-based criterion.
 * @author 	KIYOYOZU
 */

#ifndef PARTICLE_SPLIT_MERGE_H
#define PARTICLE_SPLIT_MERGE_H

#include "multiphase_mixture.h"
#include "particle_band_adaptation.h"
#include "particle_operation.h"

#include "base_body_relation.h"

namespace SPH
{
/**
 * @class ParticleSplittingByBand
 * @brief Serial split pass over the frozen candidate list; each split spawns
 * one daughter (evolving state copied from the mother) and both particles
 * take half of Vol/mass/momentum/energy, so the event conserves them exactly.
 * Capacity is pre-checked before any spawn (SpawnRealParticle has no safe
 * failure path: it increments first, checks bound after).
 */
class ParticleSplittingByBand
{
  public:
	ParticleSplittingByBand(RealBody &real_body, BaseInnerRelation &inner_relation,
							Real split_threshold, Real offset_factor,
							Real periodic_height = 0.0);
	/** @return number of split events (each adds one particle). */
	size_t exec();

  protected:
	ParticleBandAdaptation &adaptation_;
	BaseParticles &particles_;
	ParticleConfiguration &inner_configuration_;
	Real gamma_split_, lambda_, periodic_height_;
	Vecd *pos_, *vel_, *mom_;
	Real *Vol_, *mass_, *rho_, *p_, *E_, *alpha_;
	int *band_;
	Real *h_ratio_, *ref_spacing_;
	SpawnRealParticle spawn_;

	/** the split offset can push a particle across the y periodic boundary;
	 *  the cell-linked-list periodic machinery requires positions in [0, H) */
	void wrapPeriodicY(Vecd &position) const
	{
		if (periodic_height_ <= 0.0)
			return;
		if (position[1] < 0.0)
			position[1] += periodic_height_;
		else if (position[1] >= periodic_height_)
			position[1] -= periodic_height_;
	}

	size_t findNearestNeighbor(size_t index_i) const;
};

/**
 * @class ParticleMergingByBand
 * @brief Serial merge pass: candidates with gamma < gamma_m pair with their
 * nearest neighbor when the choice is mutual, the neighbor is on the same
 * side of the alpha = 0.5 material split and is not a band-0 (finest) seed.
 * The winner accumulates mass/momentum/energy/Vol; losers are marked in
 * life_status and compressed by RemoveRealParticle in one ascending pass
 * (swap-back invalidates indices, so no neighbor data may be consumed
 * afterwards until the configuration is rebuilt). A full state recovery
 * pass (vel, rho, p from the conserved state) repairs the non-evolving
 * variables of particles relocated by the compression.
 */
class ParticleMergingByBand
{
  public:
	ParticleMergingByBand(RealBody &real_body, BaseInnerRelation &inner_relation,
						  MultiphaseMixture &mixture, Real merge_threshold);
	/** @return number of merge events (each removes one particle). */
	size_t exec();

  protected:
	ParticleBandAdaptation &adaptation_;
	BaseParticles &particles_;
	ParticleConfiguration &inner_configuration_;
	MultiphaseMixture &mixture_;
	Real gamma_merge_;
	Vecd *pos_, *vel_, *mom_;
	Real *Vol_, *mass_, *rho_, *p_, *E_, *alpha_;
	int *band_;
	Real *h_ratio_, *ref_spacing_;
	RemoveRealParticle remover_;
	StdVec<int> life_status_;
	StdVec<int> partner_;

	void recoverState(size_t total_particles);
};

} // namespace SPH

#endif // PARTICLE_SPLIT_MERGE_H
