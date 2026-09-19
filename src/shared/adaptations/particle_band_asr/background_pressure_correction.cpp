/**
 * @file 	background_pressure_correction.cpp
 * @author 	KIYOYOZU
 */

#include "background_pressure_correction.h"

#include "base_body.h"
#include "base_particles.hpp"
#include "multiphase_mixture.h"
#include "particle_iterators.h"

namespace SPH
{
//=================================================================================================//
MultiphaseBackgroundPressureCorrection::MultiphaseBackgroundPressureCorrection(
	BaseInnerRelation &inner_relation, BaseContactRelation &wall_contact_relation,
	MultiphaseMixture &mixture, bool energy_correlation)
	: particles_(inner_relation.getRelation().real_body_->getBaseParticles()),
	  inner_configuration_(inner_relation.getRelation().inner_configuration_),
	  wall_contact_(wall_contact_relation),
	  mixture_(mixture),
	  energy_correlation_(energy_correlation),
	  Vol_(particles_.getVariableDataByName<Real>("VolumetricMeasure")),
	  mass_(particles_.getVariableDataByName<Real>("Mass")),
	  rho_(particles_.getVariableDataByName<Real>("Density")),
	  p_(particles_.getVariableDataByName<Real>("Pressure")),
	  E_(particles_.getVariableDataByName<Real>("TotalEnergy")),
	  alpha_(particles_.getVariableDataByName<Real>("VolumeFraction")),
	  vel_(particles_.getVariableDataByName<Vecd>("Velocity")),
	  mom_(particles_.getVariableDataByName<Vecd>("Momentum")) {}
//=================================================================================================//
void MultiphaseBackgroundPressureCorrection::snapshotState()
{
	size_t total_real_particles = particles_.TotalRealParticles();
	vel_snap_.assign(vel_, vel_ + total_real_particles);
	rho_snap_.assign(rho_, rho_ + total_real_particles);
	if (energy_correlation_)
	{
		p_snap_.assign(p_, p_ + total_real_particles);
		E_snap_.assign(E_, E_ + total_real_particles);
	}
}
//=================================================================================================//
void MultiphaseBackgroundPressureCorrection::execMomentumBetweenHalves(
	Real dt, bool upwind_select, Real interface_density_ratio)
{
	size_t total_real_particles = particles_.TotalRealParticles();
	StdVec<BaseParticles *> contact_particles = wall_contact_.getContactParticles();
	StdVec<ParticleConfiguration> &contact_configuration = wall_contact_.contact_configuration_;
	StdVec<Real *> wall_Vol;
	for (size_t k = 0; k != contact_particles.size(); ++k)
		wall_Vol.push_back(contact_particles[k]->getVariableDataByName<Real>("VolumetricMeasure"));

	particle_for(execution::ParallelPolicy(), IndexRange(0, total_real_particles),
				 [&](size_t index_i)
				 {
					 Real two_V_i_dt = 2.0 * Vol_[index_i] * dt;
					 // frozen state = step-start snapshot: rho/p/E/alpha are not
					 // touched by the 1st half, and vel must be the pre-kick one
					 const Vecd &vel_i = vel_snap_[index_i];
					 Matd T_i = rho_snap_[index_i] * vel_i * vel_i.transpose() +
								p_[index_i] * Matd::Identity();

					 Vecd mom_corr = Vecd::Zero();
					 Neighborhood &inner_neighborhood = inner_configuration_[index_i];
					 for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
					 {
						 size_t index_j = inner_neighborhood.j_[n];
						 Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
						 const Vecd &e_ij = inner_neighborhood.e_ij_[n];
						 const Vecd &vel_j = vel_snap_[index_j];
						 Matd T_j = rho_snap_[index_j] * vel_j * vel_j.transpose() +
									p_[index_j] * Matd::Identity();
						 Real w_i = 0.5, w_j = 0.5;
						 if (upwind_select &&
							 SMAX(rho_snap_[index_i], rho_snap_[index_j]) >
								 interface_density_ratio *
									 SMIN(rho_snap_[index_i], rho_snap_[index_j]))
						 {
						 // advective upwind endpoint, matching the solver's
						 // branch: its normal axis is -e_ij (i->j), so the
						 // contact s* >= 0 picks the left star state = i.
						 // s* is proxied by the mean normal velocity: exact
						 // for equilibrated contacts with continuous velocity
						 // (the band-boundary target case); a velocity-jump
						 // transient can mis-pick, bounded by O(du) flux error
						 Real u_n = -0.5 * (vel_i + vel_j).dot(e_ij);
						 w_i = u_n >= 0.0 ? 1.0 : 0.0;
						 w_j = 1.0 - w_i;
						 }
						 mom_corr += dW_ijV_j * (w_i * T_i + w_j * T_j) * e_ij;
					 }

					 // wall pairs: the frozen wall flux is p_i e. Using the exact
					 // mirror-star pressure instead (p* = p_i + rho u_n(2 u_n +- c)
					 // under wall-normal flow) is a positive feedback loop: the
					 // correction would push along +-e with gain ~ 2 rho c dt per
					 // step, amplifying any residual wall-normal velocity until
					 // divergence (measured: 1.4 m/s growth at the left wall in a quiescent run).
					 // At rest p* = p_i exactly, and under flow the mismatch is
					 // confined to wall particles, where the mirror-completed
					 // stencil sums (M_inner + M_wall) to ~0.
					 for (size_t k = 0; k != contact_configuration.size(); ++k)
					 {
						 Neighborhood &wall_neighborhood = contact_configuration[k][index_i];
						 for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
						 {
							 Real dW_ijV_j = wall_neighborhood.dW_ij_[n] *
											 wall_Vol[k][wall_neighborhood.j_[n]];
							 mom_corr += dW_ijV_j * p_[index_i] * wall_neighborhood.e_ij_[n];
						 }
					 }

					 mom_[index_i] += two_V_i_dt * mom_corr;
					 vel_[index_i] = mom_[index_i] / mass_[index_i];
				 });
}
//=================================================================================================//
void MultiphaseBackgroundPressureCorrection::execMassEnergyAfterHalves(
	Real dt, bool upwind_select, Real interface_density_ratio)
{
	if (!energy_correlation_)
		return;

	size_t total_real_particles = particles_.TotalRealParticles();
	particle_for(execution::ParallelPolicy(), IndexRange(0, total_real_particles),
				 [&](size_t index_i)
				 {
					 Real two_V_i_dt = 2.0 * Vol_[index_i] * dt;
					 // frozen state = what the 2nd-half interaction consumed:
					 // step-start rho/p/E/vel (the between-halves momentum stage
					 // restores vel to the snapshot up to ~1e-15)
					 const Vecd &vel_i = vel_snap_[index_i];
					 Vecd rho_u_i = rho_snap_[index_i] * vel_i;
					 Vecd xi_i = (E_snap_[index_i] / Vol_[index_i] + p_snap_[index_i]) * vel_i;

					 Real mass_corr_sum = 0.0;
					 Real energy_corr = 0.0;
					 Neighborhood &inner_neighborhood = inner_configuration_[index_i];
					 for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
					 {
						 size_t index_j = inner_neighborhood.j_[n];
						 Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
						 const Vecd &e_ij = inner_neighborhood.e_ij_[n];
						 const Vecd &vel_j = vel_snap_[index_j];
						 Vecd rho_u_j = rho_snap_[index_j] * vel_j;
						 Vecd xi_j = (E_snap_[index_j] / Vol_[index_j] + p_snap_[index_j]) * vel_j;
						 Real w_i = 0.5, w_j = 0.5;
						 if (upwind_select &&
							 SMAX(rho_snap_[index_i], rho_snap_[index_j]) >
								 interface_density_ratio *
									 SMIN(rho_snap_[index_i], rho_snap_[index_j]))
						 {
							 // same selection as the momentum stage: at an
							 // equilibrated contact the star mass flux rho*u and
							 // energy flux (E*+p*)u* equal the upwind endpoint's
							 // exactly (E* = E_up when p* = p_up, s* = u)
							 Real u_n = -0.5 * (vel_i + vel_j).dot(e_ij);
							 w_i = u_n >= 0.0 ? 1.0 : 0.0;
							 w_j = 1.0 - w_i;
						 }
						 mass_corr_sum += dW_ijV_j * (w_i * rho_u_i + w_j * rho_u_j).dot(e_ij);
						 energy_corr += dW_ijV_j * (w_i * xi_i + w_j * xi_j).dot(e_ij);
					 }
					 // wall pairs need no mass/energy correction: the reflective
					 // ghost star state has u* . e = 0, annihilating both fluxes

					 mass_[index_i] += two_V_i_dt * mass_corr_sum;
					 E_[index_i] += two_V_i_dt * energy_corr;

					 // full recovery from the corrected conserved state
					 rho_[index_i] = mass_[index_i] / Vol_[index_i];
					 vel_[index_i] = mom_[index_i] / mass_[index_i];
					 Real rho_e = E_[index_i] / Vol_[index_i] -
								  0.5 * rho_[index_i] * vel_[index_i].squaredNorm();
					 p_[index_i] = mixture_.MixturePressure(alpha_[index_i], rho_e);
				 });
}
//=================================================================================================//
} // namespace SPH
