#include "eulerian_multiphase_integration.h"

#include "adaptation.h"

namespace SPH
{
namespace fluid_dynamics
{

namespace
{
Vecd grad_row(const Matd &g, int r)
{
#if SPH_NDIM == 2
    return Vecd(g(r, 0), g(r, 1));
#else
    return Vecd(g(r, 0), g(r, 1), g(r, 2));
#endif
}
} // namespace
//=================================================================================================//
BaseIntegrationInMultiphase::BaseIntegrationInMultiphase(
    BaseInnerRelation &inner_relation, MultiphaseMixture &mixture)
    : BaseIntegration(inner_relation),
      mixture_(mixture),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      E_(particles_->registerStateVariableData<Real>("TotalEnergy")),
      dE_dt_(particles_->registerStateVariableData<Real>("TotalEnergyChangeRate")),
      dmass_dt_(particles_->registerStateVariableData<Real>("MassChangeRate")),
      alpha_(particles_->registerStateVariableData<Real>("VolumeFraction")),
      dalpha_dt_(particles_->registerStateVariableData<Real>("VolumeFractionChangeRate")),
      mom_(particles_->registerStateVariableData<Vecd>("Momentum")),
      force_(particles_->registerStateVariableData<Vecd>("Force")),
      force_prior_(particles_->registerStateVariableData<Vecd>("ForcePrior")),
      grad_corr_(particles_->registerStateVariableData<Vecd>("GradientCorrection")) {};
//=================================================================================================//
//	First half step: Inner<>
//=================================================================================================//
EulerianMultiphaseIntegration1stHalf<Inner<>>::EulerianMultiphaseIntegration1stHalf(
    BaseInnerRelation &inner_relation, MultiphaseMixture &mixture)
    : BaseIntegrationInMultiphase(inner_relation, mixture),
      riemann_solver_(mixture) {};
//-------------------------------------------------------------------------------------------------//
void EulerianMultiphaseIntegration1stHalf<Inner<>>::interaction(size_t index_i, Real dt)
{
    Real energy_per_volume_i = E_[index_i] / Vol_[index_i];
    MultiphaseFluidState state_i(rho_[index_i], vel_[index_i], p_[index_i],
                                   energy_per_volume_i, alpha_[index_i]);
    Vecd momentum_change_rate = force_prior_[index_i];
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];

        Real energy_per_volume_j = E_[index_j] / Vol_[index_j];
        MultiphaseFluidState state_j(rho_[index_j], vel_[index_j], p_[index_j],
                                       energy_per_volume_j, alpha_[index_j]);
        MultiphaseFluidStarState interface_state =
            riemann_solver_.getInterfaceState(state_i, state_j, e_ij);

        // zeroth-order consistent gradient: dW e V_j - c_i V_j closes the
        // first moment on graded lattices (no-op where grad_corr_ ~ 0)
        Vecd gradW_V_j = dW_ijV_j * e_ij - grad_corr_[index_i] * Vol_[index_j];
        Matd convect_flux = interface_state.rho_ * interface_state.vel_ * interface_state.vel_.transpose();
        momentum_change_rate -= 2.0 * Vol_[index_i] *
                                (convect_flux + interface_state.p_ * Matd::Identity()) * gradW_V_j;
    }
    force_[index_i] = momentum_change_rate;
}
//-------------------------------------------------------------------------------------------------//
void EulerianMultiphaseIntegration1stHalf<Inner<>>::update(size_t index_i, Real dt)
{
    mom_[index_i] += force_[index_i] * dt;
    vel_[index_i] = mom_[index_i] / mass_[index_i];
}
//=================================================================================================//
//	First half step: Contact<Wall>
//=================================================================================================//
// MSVC workaround: naming InteractionWithWall<...> directly in a member
// initializer mis-resolves the template-template argument; use an alias.
using MultiphaseWallBase = InteractionWithWall<BaseIntegrationInMultiphaseForWall>;
//-------------------------------------------------------------------------------------------------//
EulerianMultiphaseIntegration1stHalf<Contact<Wall>>::EulerianMultiphaseIntegration1stHalf(
    BaseContactRelation &contact_relation, MultiphaseMixture &mixture)
    : MultiphaseWallBase(contact_relation),
      riemann_solver_(mixture) {};
//-------------------------------------------------------------------------------------------------//
void EulerianMultiphaseIntegration1stHalf<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    Real energy_per_volume_i = this->E_[index_i] / this->Vol_[index_i];
    MultiphaseFluidState state_i(this->rho_[index_i], this->vel_[index_i], this->p_[index_i],
                                   energy_per_volume_i, this->alpha_[index_i]);
    // Prior ownership contract: Inner<> seeds force_ with ForcePrior and assigns;
    // this wall pass starts from zero and only appends wall flux.
    Vecd momentum_change_rate = Vecd::Zero();
    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        Real *Vol_k = this->wall_Vol_[k];
        Vecd *vel_ave_k = this->wall_vel_ave_[k];
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];
            Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * Vol_k[index_j];

            // Reflective ghost state: mirror velocity, reuse fluid thermodynamics.
            Vecd vel_reflect = 2.0 * vel_ave_k[index_j] - this->vel_[index_i];
            MultiphaseFluidState state_g(this->rho_[index_i], vel_reflect, this->p_[index_i],
                                           energy_per_volume_i, this->alpha_[index_i]);

            MultiphaseFluidStarState interface_state =
                riemann_solver_.getInterfaceState(state_i, state_g, e_ij);

            // fluid-side c_i closes the mirror wall stencil too
            Vecd gradW_V_j = dW_ijV_j * e_ij - grad_corr_[index_i] * Vol_k[index_j];
            Matd convect_flux = interface_state.rho_ * interface_state.vel_ * interface_state.vel_.transpose();
            momentum_change_rate -= 2.0 * this->Vol_[index_i] *
                                    (convect_flux + interface_state.p_ * Matd::Identity()) * gradW_V_j;
        }
    }
    this->force_[index_i] += momentum_change_rate;
}
//=================================================================================================//
//	Second half step: Inner<>
//=================================================================================================//
EulerianMultiphaseIntegration2ndHalf<Inner<>>::EulerianMultiphaseIntegration2ndHalf(
    BaseInnerRelation &inner_relation, MultiphaseMixture &mixture)
    : BaseIntegrationInMultiphase(inner_relation, mixture),
      riemann_solver_(mixture) {};
//-------------------------------------------------------------------------------------------------//
void EulerianMultiphaseIntegration2ndHalf<Inner<>>::interaction(size_t index_i, Real dt)
{
    Real energy_per_volume_i = E_[index_i] / Vol_[index_i];
    MultiphaseFluidState state_i(rho_[index_i], vel_[index_i], p_[index_i],
                                   energy_per_volume_i, alpha_[index_i]);
    Real mass_change_rate = 0.0;
    // Seed with the power of ForcePrior to stay consistent with the momentum
    // equation (matches the single-phase reference; TODO: not conservative
    // formulation). ForcePrior is zero here unless body forces are added.
    Real energy_change_rate = force_prior_[index_i].dot(vel_[index_i]);
    Real alpha_change_rate = 0.0;

    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];

        Real energy_per_volume_j = E_[index_j] / Vol_[index_j];
        MultiphaseFluidState state_j(rho_[index_j], vel_[index_j], p_[index_j],
                                       energy_per_volume_j, alpha_[index_j]);
        MultiphaseFluidStarState interface_state =
            riemann_solver_.getInterfaceState(state_i, state_j, e_ij);

        // zeroth-order consistent gradient for the conservative channels
        Vecd gradW_V_j = dW_ijV_j * e_ij - grad_corr_[index_i] * Vol_[index_j];

        // Mass flux: div(rho * u)
        mass_change_rate -= 2.0 * Vol_[index_i] *
                            (interface_state.rho_ * interface_state.vel_).dot(gradW_V_j);

        // Energy flux: div((rho*E + p) * u)
        energy_change_rate -= 2.0 * Vol_[index_i] *
                              ((interface_state.E_ + interface_state.p_) * interface_state.vel_).dot(gradW_V_j);

        // Volume fraction advection (non-conservative form):
        // d(alpha)/dt = -u . grad(alpha). Difference form (alpha_i - alpha*),
        // already immune to the first moment, so the raw dW_ijV_j e_ij is kept.
        Real u_star_n = interface_state.vel_.dot(e_ij);
        alpha_change_rate += 2.0 * Vol_[index_i] * dW_ijV_j * u_star_n *
                             (alpha_[index_i] - interface_state.alpha_);
    }
    dmass_dt_[index_i] = mass_change_rate;
    dE_dt_[index_i] = energy_change_rate;
    dalpha_dt_[index_i] = alpha_change_rate / Vol_[index_i];
}
//-------------------------------------------------------------------------------------------------//
void EulerianMultiphaseIntegration2ndHalf<Inner<>>::update(size_t index_i, Real dt)
{
    E_[index_i] += dE_dt_[index_i] * dt;
    mass_[index_i] += dmass_dt_[index_i] * dt;
    rho_[index_i] = mass_[index_i] / Vol_[index_i];

    // Update volume fraction and clamp to [0, 1].
    alpha_[index_i] += dalpha_dt_[index_i] * dt;
    alpha_[index_i] = SMAX(0.0, SMIN(1.0, alpha_[index_i]));

    // Recover mixture pressure from internal energy.
    Real rho_e = E_[index_i] / Vol_[index_i]
               - 0.5 * (mom_[index_i] / mass_[index_i]).squaredNorm() * rho_[index_i];
    p_[index_i] = mixture_.MixturePressure(alpha_[index_i], rho_e);
}
//=================================================================================================//
//	Second half step: Contact<Wall>
//=================================================================================================//
EulerianMultiphaseIntegration2ndHalf<Contact<Wall>>::EulerianMultiphaseIntegration2ndHalf(
    BaseContactRelation &contact_relation, MultiphaseMixture &mixture)
    : MultiphaseWallBase(contact_relation),
      riemann_solver_(mixture) {};
//-------------------------------------------------------------------------------------------------//
void EulerianMultiphaseIntegration2ndHalf<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    Real energy_per_volume_i = this->E_[index_i] / this->Vol_[index_i];
    MultiphaseFluidState state_i(this->rho_[index_i], this->vel_[index_i], this->p_[index_i],
                                   energy_per_volume_i, this->alpha_[index_i]);
    Real mass_change_rate = 0.0;
    Real energy_change_rate = 0.0;
    Real alpha_change_rate = 0.0;

    for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
    {
        Real *Vol_k = this->wall_Vol_[k];
        Vecd *vel_ave_k = this->wall_vel_ave_[k];
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];
            Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * Vol_k[index_j];

            // Reflective ghost state.
            Vecd vel_reflect = 2.0 * vel_ave_k[index_j] - this->vel_[index_i];
            MultiphaseFluidState state_g(this->rho_[index_i], vel_reflect, this->p_[index_i],
                                           energy_per_volume_i, this->alpha_[index_i]);

            MultiphaseFluidStarState interface_state =
                riemann_solver_.getInterfaceState(state_i, state_g, e_ij);

            Vecd gradW_V_j = dW_ijV_j * e_ij - grad_corr_[index_i] * Vol_k[index_j];
            mass_change_rate -= 2.0 * this->Vol_[index_i] *
                                (interface_state.rho_ * interface_state.vel_).dot(gradW_V_j);
            energy_change_rate -= 2.0 * this->Vol_[index_i] *
                                  ((interface_state.E_ + interface_state.p_) * interface_state.vel_).dot(gradW_V_j);

            Real u_star_n = interface_state.vel_.dot(e_ij);
            alpha_change_rate += 2.0 * this->Vol_[index_i] * dW_ijV_j * u_star_n *
                                 (this->alpha_[index_i] - interface_state.alpha_);
        }
    }
    this->dmass_dt_[index_i] += mass_change_rate;
    this->dE_dt_[index_i] += energy_change_rate;
    this->dalpha_dt_[index_i] += alpha_change_rate / this->Vol_[index_i];
}
//=================================================================================================//
//	Time step size
//=================================================================================================//
EulerianMultiphaseAcousticTimeStepSize::EulerianMultiphaseAcousticTimeStepSize(
    SPHBody &sph_body, MultiphaseMixture &mixture, Real acousticCFL)
    : AcousticTimeStep(sph_body),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      p_(particles_->getVariableDataByName<Real>("Pressure")),
      alpha_(particles_->getVariableDataByName<Real>("VolumeFraction")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      smoothing_length_(sph_body.getSPHAdaptation().ReferenceSmoothingLength()),
      mixture_(mixture)
{
    acousticCFL_ = acousticCFL;
};
//=================================================================================================//
Real EulerianMultiphaseAcousticTimeStepSize::reduce(size_t index_i, Real dt)
{
    return mixture_.MixtureSoundSpeed(alpha_[index_i], rho_[index_i], p_[index_i])
         + vel_[index_i].norm();
}
//=================================================================================================//
Real EulerianMultiphaseAcousticTimeStepSize::outputResult(Real reduced_value)
{
    return acousticCFL_ / Dimensions * smoothing_length_ / (reduced_value + TinyReal);
}
//=================================================================================================//
//	Volume fraction gradient (raw SPH stencil) for the MUSCL reconstruction
//=================================================================================================//
VolumeFractionGradient::VolumeFractionGradient(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      alpha_(particles_->getVariableDataByName<Real>("VolumeFraction"))
{
    particles_->registerStateVariableData<Vecd>("VolumeFractionGradient");
    alpha_grad_ = particles_->getVariableDataByName<Vecd>("VolumeFractionGradient");
}
//=================================================================================================//
void VolumeFractionGradient::interaction(size_t index_i, Real dt)
{
    Vecd alpha_grad = Vecd::Zero();
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    Real alpha_i = alpha_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd nablaW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
        alpha_grad -= (alpha_i - alpha_[index_j]) * nablaW_ijV_j;
    }
    alpha_grad_[index_i] = alpha_grad;
}
//=================================================================================================//
//	MUSCL second half step: Inner<>
//=================================================================================================//
EulerianMultiphaseIntegration1stHalfMUSCL<Inner<>>::EulerianMultiphaseIntegration1stHalfMUSCL(
    BaseInnerRelation &inner_relation, MultiphaseMixture &mixture,
    const SecondOrderConfig &cfg)
    : BaseIntegrationInMultiphase(inner_relation, mixture),
      bridge_(mixture, cfg),
      rho_grad_(particles_->getVariableDataByName<Vecd>("DensityGradient")),
      p_grad_(particles_->getVariableDataByName<Vecd>("PressureGradient")),
      alpha_grad_(particles_->getVariableDataByName<Vecd>("VolumeFractionGradient")),
      vel_grad_(particles_->getVariableDataByName<Matd>("VelocityGradient")) {}
//=================================================================================================//
void EulerianMultiphaseIntegration1stHalfMUSCL<Inner<>>::interaction(size_t index_i, Real dt)
{
    Real energy_per_volume_i = E_[index_i] / Vol_[index_i];
    MultiphaseFluidState state_i(rho_[index_i], vel_[index_i], p_[index_i],
                                   energy_per_volume_i, alpha_[index_i]);
    Vecd momentum_change_rate = force_prior_[index_i];
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];

    Matd vg_i = vel_grad_[index_i];
    Vecd grad_u_i = grad_row(vg_i, 0);
    Vecd grad_v_i = grad_row(vg_i, 1);
#if SPH_NDIM == 3
    Vecd grad_w_i = grad_row(vg_i, 2);
#endif

    const Vecd &xi = pos_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];

        Real energy_per_volume_j = E_[index_j] / Vol_[index_j];
        MultiphaseFluidState state_j(rho_[index_j], vel_[index_j], p_[index_j],
                                       energy_per_volume_j, alpha_[index_j]);

        Matd vg_j = vel_grad_[index_j];
        const Vecd xj = xi - inner_neighborhood.r_ij_[n] * e_ij;
        Vecd xf = 0.5 * (xi + xj);

        MultiphaseFluidStarState interface_state = bridge_.getInterfaceState(
            state_i, state_j, xi, xj, xf, e_ij,
            rho_grad_[index_i], rho_grad_[index_j],
            grad_u_i, grad_row(vg_j, 0),
            grad_v_i, grad_row(vg_j, 1),
#if SPH_NDIM == 3
            grad_w_i, grad_row(vg_j, 2),
#endif
            p_grad_[index_i], p_grad_[index_j],
            alpha_grad_[index_i], alpha_grad_[index_j]);

        Vecd gradW_V_j = dW_ijV_j * e_ij - grad_corr_[index_i] * Vol_[index_j];
        Matd convect_flux = interface_state.rho_ * interface_state.vel_ * interface_state.vel_.transpose();
        momentum_change_rate -= 2.0 * Vol_[index_i] *
                                (convect_flux + interface_state.p_ * Matd::Identity()) * gradW_V_j;
    }
    force_[index_i] = momentum_change_rate;
}
//=================================================================================================//
void EulerianMultiphaseIntegration1stHalfMUSCL<Inner<>>::update(size_t index_i, Real dt)
{
    mom_[index_i] += force_[index_i] * dt;
    vel_[index_i] = mom_[index_i] / mass_[index_i];
}
//=================================================================================================//
EulerianMultiphaseIntegration2ndHalfMUSCL<Inner<>>::EulerianMultiphaseIntegration2ndHalfMUSCL(
    BaseInnerRelation &inner_relation, MultiphaseMixture &mixture,
    const SecondOrderConfig &cfg)
    : BaseIntegrationInMultiphase(inner_relation, mixture),
      bridge_(mixture, cfg),
      rho_grad_(particles_->getVariableDataByName<Vecd>("DensityGradient")),
      p_grad_(particles_->getVariableDataByName<Vecd>("PressureGradient")),
      alpha_grad_(particles_->getVariableDataByName<Vecd>("VolumeFractionGradient")),
      vel_grad_(particles_->getVariableDataByName<Matd>("VelocityGradient")) {}
//=================================================================================================//
void EulerianMultiphaseIntegration2ndHalfMUSCL<Inner<>>::interaction(size_t index_i, Real dt)
{
    Real energy_per_volume_i = E_[index_i] / Vol_[index_i];
    MultiphaseFluidState state_i(rho_[index_i], vel_[index_i], p_[index_i],
                                   energy_per_volume_i, alpha_[index_i]);
    Real mass_change_rate = 0.0;
    Real energy_change_rate = force_prior_[index_i].dot(vel_[index_i]);
    Real alpha_change_rate = 0.0;
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];

    Matd vg_i = vel_grad_[index_i];
    Vecd grad_u_i = grad_row(vg_i, 0);
    Vecd grad_v_i = grad_row(vg_i, 1);
#if SPH_NDIM == 3
    Vecd grad_w_i = grad_row(vg_i, 2);
#endif

    const Vecd &xi = pos_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];

        Real energy_per_volume_j = E_[index_j] / Vol_[index_j];
        MultiphaseFluidState state_j(rho_[index_j], vel_[index_j], p_[index_j],
                                       energy_per_volume_j, alpha_[index_j]);

        Matd vg_j = vel_grad_[index_j];
        const Vecd xj = xi - inner_neighborhood.r_ij_[n] * e_ij;
        Vecd xf = 0.5 * (xi + xj);

        MultiphaseFluidStarState interface_state = bridge_.getInterfaceState(
            state_i, state_j, xi, xj, xf, e_ij,
            rho_grad_[index_i], rho_grad_[index_j],
            grad_u_i, grad_row(vg_j, 0),
            grad_v_i, grad_row(vg_j, 1),
#if SPH_NDIM == 3
            grad_w_i, grad_row(vg_j, 2),
#endif
            p_grad_[index_i], p_grad_[index_j],
            alpha_grad_[index_i], alpha_grad_[index_j]);

        Vecd gradW_V_j = dW_ijV_j * e_ij - grad_corr_[index_i] * Vol_[index_j];
        mass_change_rate -= 2.0 * Vol_[index_i] *
                            (interface_state.rho_ * interface_state.vel_).dot(gradW_V_j);
        energy_change_rate -= 2.0 * Vol_[index_i] *
                              ((interface_state.E_ + interface_state.p_) * interface_state.vel_).dot(gradW_V_j);

        Real u_star_n = interface_state.vel_.dot(e_ij);
        alpha_change_rate += 2.0 * Vol_[index_i] * dW_ijV_j * u_star_n *
                             (alpha_[index_i] - interface_state.alpha_);
    }
    dmass_dt_[index_i] = mass_change_rate;
    dE_dt_[index_i] = energy_change_rate;
    dalpha_dt_[index_i] = alpha_change_rate / Vol_[index_i];
}
//=================================================================================================//
void EulerianMultiphaseIntegration2ndHalfMUSCL<Inner<>>::update(size_t index_i, Real dt)
{
    E_[index_i] += dE_dt_[index_i] * dt;
    mass_[index_i] += dmass_dt_[index_i] * dt;
    rho_[index_i] = mass_[index_i] / Vol_[index_i];

    alpha_[index_i] += dalpha_dt_[index_i] * dt;
    alpha_[index_i] = SMAX(0.0, SMIN(1.0, alpha_[index_i]));

    Real rho_e = E_[index_i] / Vol_[index_i]
               - 0.5 * (mom_[index_i] / mass_[index_i]).squaredNorm() * rho_[index_i];
    p_[index_i] = mixture_.MixturePressure(alpha_[index_i], rho_e);
}
//=================================================================================================//
//	MUSCL wall (reflective ghost) pass, mirroring the single-phase
//	EulerianCompressibleIntegration*MUSCL<Contact<Wall>> pattern: the wall-side
//	state is the fluid particle mirrored about the wall velocity, carrying the
//	fluid-side gradients, reconstructed to the wall midpoint by the same bridge.
//=================================================================================================//
using MultiphaseMUSCLWallBase = InteractionWithWall<BaseIntegrationInMultiphaseForWall>;
//=================================================================================================//
EulerianMultiphaseIntegration1stHalfMUSCL<Contact<Wall>>::EulerianMultiphaseIntegration1stHalfMUSCL(
    BaseContactRelation &contact_relation, MultiphaseMixture &mixture,
    const SecondOrderConfig &cfg)
    : MultiphaseMUSCLWallBase(contact_relation),
      bridge_(mixture, cfg),
      rho_grad_(particles_->getVariableDataByName<Vecd>("DensityGradient")),
      p_grad_(particles_->getVariableDataByName<Vecd>("PressureGradient")),
      alpha_grad_(particles_->getVariableDataByName<Vecd>("VolumeFractionGradient")),
      vel_grad_(particles_->getVariableDataByName<Matd>("VelocityGradient")) {}
//=================================================================================================//
void EulerianMultiphaseIntegration1stHalfMUSCL<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    Real energy_per_volume_i = this->E_[index_i] / this->Vol_[index_i];
    MultiphaseFluidState state_i(this->rho_[index_i], this->vel_[index_i], this->p_[index_i],
                                   energy_per_volume_i, this->alpha_[index_i]);
    // Prior ownership contract (same as the single-phase wall pass): Inner<>
    // seeds force_ with ForcePrior and assigns; this wall pass starts from zero
    // and only appends the wall flux, so the prior is counted exactly once.
    Vecd momentum_change_rate = Vecd::Zero();

    Matd vg_i = vel_grad_[index_i];
    Vecd grad_u_i = grad_row(vg_i, 0);
    Vecd grad_v_i = grad_row(vg_i, 1);
#if SPH_NDIM == 3
    Vecd grad_w_i = grad_row(vg_i, 2);
#endif
    // rho and alpha stay piecewise constant (zero slope), as in the inner pass.
    const Vecd zero = Vecd::Zero();

    const Vecd &xi = this->pos_[index_i];
    for (size_t k = 0; k != this->contact_configuration_.size(); ++k)
    {
        Real *Vol_k = this->wall_Vol_[k];
        Vecd *vel_ave_k = this->wall_vel_ave_[k];
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];
            Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * Vol_k[index_j];

            Vecd vel_reflect = 2.0 * vel_ave_k[index_j] - this->vel_[index_i];
            MultiphaseFluidState state_g(this->rho_[index_i], vel_reflect, this->p_[index_i],
                                           energy_per_volume_i, this->alpha_[index_i]);

            Vecd xj = xi - contact_neighborhood.r_ij_[n] * e_ij;
            Vecd xf = 0.5 * (xi + xj);

            MultiphaseFluidStarState interface_state = bridge_.getInterfaceState(
                state_i, state_g, xi, xj, xf, e_ij,
                zero, zero,
                grad_u_i, grad_u_i,
                grad_v_i, grad_v_i,
#if SPH_NDIM == 3
                grad_w_i, grad_w_i,
#endif
                p_grad_[index_i], p_grad_[index_i],
                zero, zero);

            Vecd gradW_V_j = dW_ijV_j * e_ij - grad_corr_[index_i] * Vol_k[index_j];
            Matd convect_flux = interface_state.rho_ * interface_state.vel_ * interface_state.vel_.transpose();
            momentum_change_rate -= 2.0 * this->Vol_[index_i] *
                                    (convect_flux + interface_state.p_ * Matd::Identity()) * gradW_V_j;
        }
    }
    this->force_[index_i] += momentum_change_rate;
}
//=================================================================================================//
EulerianMultiphaseIntegration2ndHalfMUSCL<Contact<Wall>>::EulerianMultiphaseIntegration2ndHalfMUSCL(
    BaseContactRelation &contact_relation, MultiphaseMixture &mixture,
    const SecondOrderConfig &cfg)
    : MultiphaseMUSCLWallBase(contact_relation),
      bridge_(mixture, cfg),
      rho_grad_(particles_->getVariableDataByName<Vecd>("DensityGradient")),
      p_grad_(particles_->getVariableDataByName<Vecd>("PressureGradient")),
      alpha_grad_(particles_->getVariableDataByName<Vecd>("VolumeFractionGradient")),
      vel_grad_(particles_->getVariableDataByName<Matd>("VelocityGradient")) {}
//=================================================================================================//
void EulerianMultiphaseIntegration2ndHalfMUSCL<Contact<Wall>>::interaction(size_t index_i, Real dt)
{
    Real energy_per_volume_i = this->E_[index_i] / this->Vol_[index_i];
    MultiphaseFluidState state_i(this->rho_[index_i], this->vel_[index_i], this->p_[index_i],
                                   energy_per_volume_i, this->alpha_[index_i]);
    // Prior ownership contract: Inner<> seeds dE_dt_ with the prior work and
    // assigns; this wall pass starts from zero and only appends the wall flux.
    Real mass_change_rate = 0.0;
    Real energy_change_rate = 0.0;
    Real alpha_change_rate = 0.0;

    Matd vg_i = vel_grad_[index_i];
    Vecd grad_u_i = grad_row(vg_i, 0);
    Vecd grad_v_i = grad_row(vg_i, 1);
#if SPH_NDIM == 3
    Vecd grad_w_i = grad_row(vg_i, 2);
#endif
    const Vecd zero = Vecd::Zero();

    const Vecd &xi = this->pos_[index_i];
    for (size_t k = 0; k != this->contact_configuration_.size(); ++k)
    {
        Real *Vol_k = this->wall_Vol_[k];
        Vecd *vel_ave_k = this->wall_vel_ave_[k];
        Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];
            Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * Vol_k[index_j];

            Vecd vel_reflect = 2.0 * vel_ave_k[index_j] - this->vel_[index_i];
            MultiphaseFluidState state_g(this->rho_[index_i], vel_reflect, this->p_[index_i],
                                           energy_per_volume_i, this->alpha_[index_i]);

            Vecd xj = xi - contact_neighborhood.r_ij_[n] * e_ij;
            Vecd xf = 0.5 * (xi + xj);

            MultiphaseFluidStarState interface_state = bridge_.getInterfaceState(
                state_i, state_g, xi, xj, xf, e_ij,
                zero, zero,
                grad_u_i, grad_u_i,
                grad_v_i, grad_v_i,
#if SPH_NDIM == 3
                grad_w_i, grad_w_i,
#endif
                p_grad_[index_i], p_grad_[index_i],
                zero, zero);

            Vecd gradW_V_j = dW_ijV_j * e_ij - grad_corr_[index_i] * Vol_k[index_j];
            mass_change_rate -= 2.0 * this->Vol_[index_i] *
                                (interface_state.rho_ * interface_state.vel_).dot(gradW_V_j);
            energy_change_rate -= 2.0 * this->Vol_[index_i] *
                                  ((interface_state.E_ + interface_state.p_) * interface_state.vel_).dot(gradW_V_j);

            // Identically zero for the mirrored ghost (alpha_g = alpha_i), kept
            // for structural parity with the first-order wall pass. Difference
            // form, so it keeps the raw dW_ijV_j e_ij (first-moment immune).
            Real u_star_n = interface_state.vel_.dot(e_ij);
            alpha_change_rate += 2.0 * this->Vol_[index_i] * dW_ijV_j * u_star_n *
                                 (this->alpha_[index_i] - interface_state.alpha_);
        }
    }
    this->dmass_dt_[index_i] += mass_change_rate;
    this->dE_dt_[index_i] += energy_change_rate;
    this->dalpha_dt_[index_i] += alpha_change_rate / this->Vol_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
