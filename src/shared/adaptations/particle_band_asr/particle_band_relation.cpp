/**
 * @file 	particle_band_relation.cpp
 * @author 	KIYOYOZU
 */

#include "particle_band_relation.h"

#include "adaptation.h"
#include "base_particles.hpp"
#include "cell_linked_list.hpp"

namespace SPH
{
//=================================================================================================//
NeighborBuilderInnerSymmetric::NeighborBuilderInnerSymmetric(SPHBody &body)
    : NeighborBuilder(body.getSPHAdaptation().getKernel()),
      h_ratio_(body.getBaseParticles().getVariableDataByName<Real>("SmoothingLengthRatio")),
      inv_h_ref_(1.0 / kernel_->SmoothingLength()),
      kernel_size_(kernel_->KernelSize()) {}
//=================================================================================================//
void NeighborBuilderInnerSymmetric::
operator()(Neighborhood &neighborhood, const Vecd &pos_i, size_t index_i, const ListData &list_data_j)
{
    size_t index_j = list_data_j.first;
    Vecd displacement = pos_i - list_data_j.second;
    Real distance = displacement.norm();
    Real i_h_ratio = h_ratio_[index_i];
    Real h_ratio_min = SMIN(i_h_ratio, h_ratio_[index_j]);
    // pair support is the larger of both kernels' supports
    if (distance < kernel_->CutOffRadius(h_ratio_min) && index_i != index_j)
    {
        Real j_h_ratio = h_ratio_[index_j];
        Real q_i = distance * inv_h_ref_ * i_h_ratio;
        Real q_j = distance * inv_h_ref_ * j_h_ratio;
        // per-term truncation: library kernels do not zero out beyond support
        Real W_ij = q_i < kernel_size_ ? kernel_->W(i_h_ratio, distance, displacement) : 0.0;
        Real dW_i = q_i < kernel_size_ ? kernel_->dW(i_h_ratio, distance, displacement) : 0.0;
        Real dW_j = q_j < kernel_size_ ? kernel_->dW(j_h_ratio, distance, displacement) : 0.0;
        Real dW_ij = 0.5 * (dW_i + dW_j); // symmetric under i <-> j, Eq. (3)
        Vecd e_ij = displacement / (distance + TinyReal);

        if (neighborhood.current_size_ >= neighborhood.allocated_size_)
        {
            neighborhood.j_.push_back(index_j);
            neighborhood.W_ij_.push_back(W_ij);
            neighborhood.dW_ij_.push_back(dW_ij);
            neighborhood.r_ij_.push_back(distance);
            neighborhood.e_ij_.push_back(e_ij);
            neighborhood.allocated_size_++;
        }
        else
        {
            size_t n = neighborhood.current_size_;
            neighborhood.j_[n] = index_j;
            neighborhood.W_ij_[n] = W_ij;
            neighborhood.dW_ij_[n] = dW_ij;
            neighborhood.r_ij_[n] = distance;
            neighborhood.e_ij_[n] = e_ij;
        }
        neighborhood.current_size_++;
    }
}
//=================================================================================================//
ParticleBandInnerRelation::ParticleBandInnerRelation(RealBody &real_body)
    : AdaptiveInnerRelation(real_body),
      get_symmetric_inner_neighbor_(real_body) {}
//=================================================================================================//
void ParticleBandInnerRelation::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    Mesh *meshes = multi_level_cell_linked_list_.getMeshes();
    for (size_t l = 0; l != multi_level_cell_linked_list_.ResolutionLevels(); ++l)
    {
        multi_level_cell_linked_list_.searchNeighborsByMesh(
            meshes[l], sph_body_, inner_configuration_,
            *get_multi_level_search_depth_[l], get_symmetric_inner_neighbor_);
    }
}
//=================================================================================================//
} // namespace SPH
