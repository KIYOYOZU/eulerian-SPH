/**
 * @file 	particle_band_relation.h
 * @brief 	Symmetric inner neighbor relation for the particle-band ASR method.
 * @details Caches the symmetric kernel gradient of Eq. (3) in
 * 			Yang, Kong & Liu, Phys. Rev. E 104, 055308 (2021):
 * 			dW_ij = 0.5 [dW(r, h_i) + dW(r, h_j)], with per-term truncation
 * 			at the respective kernel support (library kernels do not truncate
 * 			themselves; callers must guard q >= kernel_size). The pair is built
 * 			when r < max(rc_i, rc_j) = rc_ref / min(h_ratio_i, h_ratio_j), so
 * 			that both partner lists contain the pair and the cached dW_ij is
 * 			identical from either side, which is the precondition for exact
 * 			global conservation of the pairwise-antisymmetric ESPH fluxes.
 * @author 	KIYOYOZU
 */

#ifndef PARTICLE_BAND_RELATION_H
#define PARTICLE_BAND_RELATION_H

#include "inner_body_relation.h"

namespace SPH
{
/**
 * @class NeighborBuilderInnerSymmetric
 * @brief Inner neighbor builder caching W(r, h_i) and the symmetric dW of Eq. (3).
 */
class NeighborBuilderInnerSymmetric : public NeighborBuilder
{
  public:
    explicit NeighborBuilderInnerSymmetric(SPHBody &body);
    void operator()(Neighborhood &neighborhood,
                    const Vecd &pos_i, size_t index_i, const ListData &list_data_j) override;

  protected:
    Real *h_ratio_;
    Real inv_h_ref_, kernel_size_;
};

/**
 * @class ParticleBandInnerRelation
 * @brief Inner relation building configurations with the symmetric kernel gradient.
 */
class ParticleBandInnerRelation : public AdaptiveInnerRelation
{
  public:
    explicit ParticleBandInnerRelation(RealBody &real_body);
    virtual ~ParticleBandInnerRelation() {};

    virtual void updateConfiguration() override;

  protected:
    NeighborBuilderInnerSymmetric get_symmetric_inner_neighbor_;
};
} // namespace SPH
#endif // PARTICLE_BAND_RELATION_H
