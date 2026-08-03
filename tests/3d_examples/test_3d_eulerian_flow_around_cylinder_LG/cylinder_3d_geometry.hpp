#ifndef CYLINDER_3D_GEOMETRY_HPP
#define CYLINDER_3D_GEOMETRY_HPP

#include "cylinder_3d_data.hpp"
#include "eulerian_open_boundary.h"
#include "sphinxsys.h"

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace SPH
{
namespace cylinder_3d
{

class Cylinder3DFluidBlock : public ComplexShape
{
  public:
    explicit Cylinder3DFluidBlock(const std::string &shape_name, const Cylinder3DConfig &cfg)
        : ComplexShape(shape_name)
    {
        const Vecd halfsize(0.5 * (cfg.DL + 2.0 * cfg.sponge_width),
                            0.5 * (cfg.DH + 2.0 * cfg.sponge_width),
                            0.5 * cfg.DW);
        const Vecd translation(0.5 * cfg.DL, 0.5 * cfg.DH, 0.5 * cfg.DW);
        add<GeometricShapeBox>(Transform(translation), halfsize, "OuterBoundary");

        const Vecd cylinder_center(cfg.cylinder_center_x, cfg.cylinder_center_y, 0.5 * cfg.DW);
        subtract<TriangleMeshShapeCylinder>(Vecd(0.0, 0.0, 1.0), cfg.cylinder_radius,
                                            0.5 * cfg.DW, 48,
                                            cylinder_center, "CylinderVoid");
    }
};

class Cylinder3DWallBlock : public ComplexShape
{
  public:
    explicit Cylinder3DWallBlock(const std::string &shape_name, const Cylinder3DConfig &cfg)
        : ComplexShape(shape_name)
    {
        const Vecd cylinder_center(cfg.cylinder_center_x, cfg.cylinder_center_y, 0.5 * cfg.DW);
        add<TriangleMeshShapeCylinder>(Vecd(0.0, 0.0, 1.0), cfg.cylinder_radius,
                                       0.5 * cfg.DW, 48,
                                       cylinder_center, "CylinderWall");
    }
};

class Cylinder3DPeriodicWallNormal : public NormalDirectionFromBodyShape
{
  public:
    Cylinder3DPeriodicWallNormal(SPHBody &sph_body, const Cylinder3DConfig &cfg)
        : NormalDirectionFromBodyShape(sph_body), cfg_(cfg) {}

    void update(size_t index_i, Real dt = 0.0)
    {
        Vecd radial(pos_[index_i][0] - cfg_.cylinder_center_x,
                    pos_[index_i][1] - cfg_.cylinder_center_y, 0.0);
        const Real radial_distance = radial.norm();
        const Vecd normal = radial_distance > TinyReal ? radial / radial_distance : Vecd(1.0, 0.0, 0.0);
        const Real signed_distance = radial_distance - cfg_.cylinder_radius;
        n_[index_i] = normal;
        n0_[index_i] = normal;
        phi_[index_i] = signed_distance;
        phi0_[index_i] = signed_distance;
    }

  private:
    const Cylinder3DConfig &cfg_;
};

class Cylinder3DInitialCondition : public fluid_dynamics::FluidInitialCondition
{
  public:
    explicit Cylinder3DInitialCondition(SPHBody &sph_body, const Cylinder3DConfig &cfg)
        : fluid_dynamics::FluidInitialCondition(sph_body),
          cfg_(cfg),
          rho_(particles_->registerStateVariableData<Real>("Density")),
          p_(particles_->registerStateVariableData<Real>("Pressure")),
          mass_(particles_->registerStateVariableData<Real>("Mass")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          mom_(particles_->registerStateVariableData<Vecd>("Momentum")) {};

    void update(size_t index_i, Real dt = 0.0)
    {
        rho_[index_i] = cfg_.rho0_f;
        p_[index_i] = 0.0;
        Real streamwise_velocity = cfg_.u_f;
        const Vecd &pos_i = pos_[index_i];
        if (pos_i[0] >= 0.0 && pos_i[0] <= cfg_.DL &&
            pos_i[1] >= 0.0 && pos_i[1] <= cfg_.DH)
        {
            streamwise_velocity *= 1.0 + 0.15 * signedIndexNoise(index_i);
        }
        vel_[index_i] = Vecd(streamwise_velocity, 0.0, 0.0);
        mass_[index_i] = rho_[index_i] * Vol_[index_i];
        mom_[index_i] = mass_[index_i] * vel_[index_i];
    }

  private:
    static Real signedIndexNoise(size_t index)
    {
        uint64_t bits = static_cast<uint64_t>(index) + 0x9e3779b97f4a7c15ULL;
        bits = (bits ^ (bits >> 30)) * 0xbf58476d1ce4e5b9ULL;
        bits = (bits ^ (bits >> 27)) * 0x94d049bb133111ebULL;
        bits ^= bits >> 31;
        return 2.0 * static_cast<Real>(bits >> 11) * (1.0 / 9007199254740992.0) - 1.0;
    }

    const Cylinder3DConfig &cfg_;
    Real *rho_;
    Real *p_;
    Real *mass_;
    Real *Vol_;
    Vecd *mom_;
};

inline bool isFiniteReal(Real value)
{
    return std::isfinite(static_cast<double>(value));
}

inline bool isFiniteVec(const Vecd &value)
{
    for (int axis = 0; axis != Dimensions; ++axis)
    {
        if (!isFiniteReal(value[axis]))
        {
            return false;
        }
    }
    return true;
}

inline bool insideCylinderCore(const Vecd &pos, const Cylinder3DConfig &cfg, Real radius_offset = 0.0)
{
    const Real dx = pos[0] - cfg.cylinder_center_x;
    const Real dy = pos[1] - cfg.cylinder_center_y;
    const Real r = std::sqrt(dx * dx + dy * dy);
    return r <= cfg.cylinder_radius + radius_offset;
}

inline Vecd orientedWallNormal(const Vecd &e_ij, const Vecd &wall_normal)
{
    return SGN(e_ij.dot(wall_normal)) * wall_normal;
}

inline Vecd projectForceToWallNormal(const Vecd &force, const Vecd &wall_normal)
{
    const Real normal_norm = wall_normal.norm();
    if (normal_norm <= Real(1.0e-12))
    {
        return force;
    }
    const Vecd unit_normal = wall_normal / normal_norm;
    return force.dot(unit_normal) * unit_normal;
}

template <class RiemannSolverType>
class EulerianPressureIntegration1stHalfWallContactNormalOnly
    : public fluid_dynamics::BaseEulerianIntegrationWithWall
{
  public:
    explicit EulerianPressureIntegration1stHalfWallContactNormalOnly(
        BaseContactRelation &wall_contact_relation, Real limiter_parameter = 15.0)
        : fluid_dynamics::BaseEulerianIntegrationWithWall(wall_contact_relation),
          riemann_solver_(this->fluid_, this->fluid_, limiter_parameter) {}

    void interaction(size_t index_i, Real dt = 0.0)
    {
        FluidStateIn state_i(this->rho_[index_i], this->vel_[index_i], this->p_[index_i]);
        Vecd momentum_change_rate = Vecd::Zero();
        for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
        {
            Vecd *n_k = this->wall_n_[k];
            Real *Vol_k = this->wall_Vol_[k];
            Vecd *vel_ave_k = this->wall_vel_ave_[k];
            Neighborhood &wall_neighborhood = (*this->contact_configuration_[k])[index_i];
            for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
            {
                const size_t index_j = wall_neighborhood.j_[n];
                const Vecd &e_ij = wall_neighborhood.e_ij_[n];
                const Real dW_ijV_j = wall_neighborhood.dW_ij_[n] * Vol_k[index_j];

                Real rho_j_in_wall = state_i.rho_;
                Real p_j_in_wall = state_i.p_;
                Vecd vel_j_in_wall = 2.0 * vel_ave_k[index_j] - state_i.vel_;
                FluidStateIn state_j(rho_j_in_wall, vel_j_in_wall, p_j_in_wall);
                FluidStateOut interface_state = riemann_solver_.InterfaceState(state_i, state_j, n_k[index_j]);
                const Matd convect_flux = interface_state.rho_ * interface_state.vel_ * interface_state.vel_.transpose();
                const Vecd convective_force =
                    -2.0 * this->Vol_[index_i] * convect_flux * e_ij * dW_ijV_j;
                const Vecd raw_pressure_force =
                    -2.0 * this->Vol_[index_i] * interface_state.p_ * e_ij * dW_ijV_j;
                const Vecd face_to_fluid_n = orientedWallNormal(e_ij, n_k[index_j]);
                momentum_change_rate += convective_force + projectForceToWallNormal(raw_pressure_force, face_to_fluid_n);
            }
        }
        this->dmom_dt_[index_i] += momentum_change_rate;
    }

  protected:
    RiemannSolverType riemann_solver_;
};

template <class RiemannSolverType>
class EulerianPressureIntegration1stHalfWithWallNormalOnly
    : public fluid_dynamics::EulerianIntegration1stHalf<Inner<>, RiemannSolverType>
{
    using InnerIntegration = fluid_dynamics::EulerianIntegration1stHalf<Inner<>, RiemannSolverType>;

  public:
    explicit EulerianPressureIntegration1stHalfWithWallNormalOnly(
        BaseInnerRelation &inner_relation, BaseContactRelation &wall_contact_relation,
        Real limiter_parameter = 15.0)
        : InnerIntegration(inner_relation, limiter_parameter),
          wall_contact_integration_(wall_contact_relation, limiter_parameter) {}

    using InnerIntegration::update;

    void initialization(size_t index_i, Real dt = 0.0) {}

    void interaction(size_t index_i, Real dt = 0.0)
    {
        InnerIntegration::interaction(index_i, dt);
        wall_contact_integration_.interaction(index_i, dt);
    }

  protected:
    EulerianPressureIntegration1stHalfWallContactNormalOnly<RiemannSolverType>
        wall_contact_integration_;
};

using EulerianPressureRelaxationWithWallNormalOnly =
    EulerianPressureIntegration1stHalfWithWallNormalOnly<AcousticRiemannSolver>;

class Cylinder3DFarFieldBoundary : public fluid_dynamics::NonReflectiveBoundaryCorrection
{
  public:
    explicit Cylinder3DFarFieldBoundary(BaseInnerRelation &inner_relation, const Cylinder3DConfig &cfg)
        : fluid_dynamics::NonReflectiveBoundaryCorrection(inner_relation),
          cfg_(cfg),
          boundary_state_(makeEulerianWeaklyCompressibleBoundaryState(*particles_)),
          wall_skip_band_(static_cast<Real>(cfg.boundary_n_layers) * cfg.global_resolution)
    {
        rho_farfield_ = cfg.rho0_f;
        sound_speed_ = cfg.c_f;
        vel_farfield_ = freestreamVelocity(cfg);
    };
    virtual ~Cylinder3DFarFieldBoundary() {};

    Vecd farfieldVelocity(size_t index_i) const
    {
        return freestreamVelocity(cfg_);
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        const Vecd &pos_i = pos_[index_i];
        if (isNearCylinderWall(pos_i))
        {
            return;
        }

        if (isOutsidePhysicalOpenBox(pos_i))
        {
            rho_[index_i] = rho_farfield_;
            vel_[index_i] = vel_farfield_;
            syncEulerianWeaklyCompressibleState(boundary_state_, index_i);
            return;
        }

        if (indicator_[index_i] != 1 && smeared_surface_[index_i] != 1)
        {
            return;
        }

        const Vecd &n_i = n_[index_i];
        const int face = dominantOpenFace(n_i);
        if (face == 2)
        {
            return;
        }

        const Vecd vel_ff = farfieldVelocity(index_i);
        const Real velocity_farfield_normal = vel_ff.dot(n_i);
        const Real velocity_boundary_normal = vel_[index_i].dot(n_i);
        Real w = inner_weight_summation_[index_i];
        if (!isFiniteReal(w))
        {
            w = 0.0;
        }
        w = SMIN(Real(1.0), SMAX(Real(0.0), w));
        const Real one_minus_w = 1.0 - w;
        const bool inflow = face < 0;

        if (std::fabs(velocity_boundary_normal) >= sound_speed_)
        {
            if (inflow)
            {
                rho_[index_i] = rho_farfield_;
                vel_[index_i] = vel_ff;
            }
            else
            {
                rho_[index_i] = rho_average_[index_i] + TinyReal;
                vel_[index_i] = vel_average_[index_i];
            }
            syncEulerianWeaklyCompressibleState(boundary_state_, index_i);
            return;
        }

        rho_[index_i] = rho_average_[index_i] * w + rho_farfield_ * one_minus_w;
        const Real vel_normal =
            vel_normal_average_[index_i] * w + velocity_farfield_normal * one_minus_w;
        if (inflow)
        {
            vel_[index_i] = vel_normal * n_i + (vel_ff - velocity_farfield_normal * n_i);
        }
        else
        {
            vel_[index_i] = vel_normal * n_i + vel_tangential_average_[index_i];
        }
        syncEulerianWeaklyCompressibleState(boundary_state_, index_i);
    }

  private:
    int dominantOpenFace(const Vecd &normal) const
    {
        const Real ax = std::fabs(normal[0]);
        const Real ay = std::fabs(normal[1]);
        const Real az = std::fabs(normal[2]);
        if (az > ax && az > ay)
        {
            return 2;
        }
        if (ax >= ay)
        {
            return normal[0] < 0.0 ? -1 : 1;
        }
        return -1;
    }

    bool isOutsidePhysicalOpenBox(const Vecd &pos) const
    {
        return pos[0] < 0.0 || pos[0] > cfg_.DL || pos[1] < 0.0 || pos[1] > cfg_.DH;
    }

    bool isNearCylinderWall(const Vecd &pos) const
    {
        return insideCylinderCore(pos, cfg_, wall_skip_band_);
    }

    const Cylinder3DConfig &cfg_;
    EulerianWeaklyCompressibleBoundaryState boundary_state_;
    Real wall_skip_band_;
};

inline bool reportFluidFiniteState(BaseParticles &particles, const std::string &stage)
{
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    Real *pressure = particles.getVariableDataByName<Real>("Pressure");
    Vecd *velocity = particles.getVariableDataByName<Vecd>("Velocity");
    Vecd *momentum = particles.getVariableDataByName<Vecd>("Momentum");
    const size_t total = particles.TotalRealParticles();

    Real rho_min = std::numeric_limits<Real>::max();
    Real rho_max = -std::numeric_limits<Real>::max();
    Real p_min = std::numeric_limits<Real>::max();
    Real p_max = -std::numeric_limits<Real>::max();
    Real speed_max = 0.0;
    size_t first_bad = total;
    size_t rho_min_index = total;
    size_t rho_max_index = total;
    size_t speed_max_index = total;

    for (size_t i = 0; i != total; ++i)
    {
        const bool finite = isFiniteVec(pos[i]) && isFiniteReal(vol[i]) &&
                            isFiniteReal(rho[i]) && isFiniteReal(mass[i]) &&
                            isFiniteReal(pressure[i]) && isFiniteVec(velocity[i]) &&
                            isFiniteVec(momentum[i]);
        if (!finite && first_bad == total)
        {
            first_bad = i;
        }
        if (isFiniteReal(rho[i]))
        {
            if (rho[i] < rho_min)
            {
                rho_min = rho[i];
                rho_min_index = i;
            }
            if (rho[i] > rho_max)
            {
                rho_max = rho[i];
                rho_max_index = i;
            }
        }
        if (isFiniteReal(pressure[i]))
        {
            p_min = SMIN(p_min, pressure[i]);
            p_max = SMAX(p_max, pressure[i]);
        }
        if (isFiniteVec(velocity[i]))
        {
            const Real speed = velocity[i].norm();
            if (speed > speed_max)
            {
                speed_max = speed;
                speed_max_index = i;
            }
        }
    }

    std::cout << "[Cylinder3D][State] " << stage
              << " total=" << total
              << " rho=[" << rho_min << ", " << rho_max << "]"
              << " p=[" << p_min << ", " << p_max << "]"
              << " max|u|=" << speed_max
              << " finite=" << (first_bad == total ? "yes" : "NO") << std::endl;
    if (first_bad != total)
    {
        std::cout << "[Cylinder3D][BadParticle] i=" << first_bad
                  << " pos=(" << pos[first_bad][0] << ", " << pos[first_bad][1]
                  << ", " << pos[first_bad][2] << ")"
                  << " rho=" << rho[first_bad]
                  << " p=" << pressure[first_bad]
                  << " vel=(" << velocity[first_bad][0] << ", " << velocity[first_bad][1]
                  << ", " << velocity[first_bad][2] << ")" << std::endl;
    }
    if (speed_max_index != total)
    {
        std::cout << "[Cylinder3D][MaxSpeedParticle] i=" << speed_max_index
                  << " pos=(" << pos[speed_max_index][0] << ", " << pos[speed_max_index][1]
                  << ", " << pos[speed_max_index][2] << ")"
                  << " rho=" << rho[speed_max_index]
                  << " p=" << pressure[speed_max_index]
                  << " vel=(" << velocity[speed_max_index][0] << ", " << velocity[speed_max_index][1]
                  << ", " << velocity[speed_max_index][2] << ")" << std::endl;
    }
    if (rho_min_index != total && rho_max_index != total)
    {
        std::cout << "[Cylinder3D][RhoExtrema]"
                  << " min_i=" << rho_min_index
                  << " min_pos=(" << pos[rho_min_index][0] << ", " << pos[rho_min_index][1]
                  << ", " << pos[rho_min_index][2] << ")"
                  << " max_i=" << rho_max_index
                  << " max_pos=(" << pos[rho_max_index][0] << ", " << pos[rho_max_index][1]
                  << ", " << pos[rho_max_index][2] << ")" << std::endl;
    }
    return first_bad == total;
}

struct GeometryGateResult
{
    size_t fluid_particles = 0;
    size_t wall_particles = 0;
    size_t cylinder_internal_fluid = 0;
    size_t inlet = 0;
    size_t outlet = 0;
    size_t top_bottom = 0;
    size_t z_face = 0;
    size_t cylinder_contact_particles = 0;
    Real z_end_inner_mean = 0.0;
    Real z_mid_inner_mean = 0.0;
    Real z_end_contact_mean = 0.0;
    Real z_mid_contact_mean = 0.0;
    Real wall_z_min = std::numeric_limits<Real>::infinity();
    Real wall_z_max = -std::numeric_limits<Real>::infinity();
    Real wall_max_abs_normal_z = 0.0;
    bool pass = false;
};

inline GeometryGateResult reportGeometryGate(FluidBody &fluid_body, SolidBody &wall_body,
                                             BaseInnerRelation &fluid_inner,
                                             ContactRelation &fluid_wall_contact,
                                             const Cylinder3DConfig &cfg)
{
    GeometryGateResult gate;
    BaseParticles &fluid_particles = fluid_body.getBaseParticles();
    BaseParticles &wall_particles = wall_body.getBaseParticles();
    Vecd *pos = fluid_particles.getVariableDataByName<Vecd>("Position");
    Vecd *normal = fluid_particles.getVariableDataByName<Vecd>("NormalDirection");
    Vecd *wall_pos = wall_particles.getVariableDataByName<Vecd>("Position");
    Vecd *wall_normal = wall_particles.getVariableDataByName<Vecd>("NormalDirection");
    gate.fluid_particles = fluid_particles.TotalRealParticles();
    gate.wall_particles = wall_particles.TotalRealParticles();

    size_t z_end_count = 0;
    size_t z_mid_count = 0;
    Real z_end_inner_sum = 0.0;
    Real z_mid_inner_sum = 0.0;
    Real z_end_sum = 0.0;
    Real z_mid_sum = 0.0;
    const Real z_band = static_cast<Real>(cfg.boundary_n_layers) * cfg.dp;
    const Real mid_low = 0.5 * cfg.DW - z_band;
    const Real mid_high = 0.5 * cfg.DW + z_band;

    for (size_t i = 0; i != gate.wall_particles; ++i)
    {
        gate.wall_z_min = SMIN(gate.wall_z_min, wall_pos[i][2]);
        gate.wall_z_max = SMAX(gate.wall_z_max, wall_pos[i][2]);
        gate.wall_max_abs_normal_z = SMAX(gate.wall_max_abs_normal_z, std::fabs(wall_normal[i][2]));
    }

    for (size_t i = 0; i != gate.fluid_particles; ++i)
    {
        if (insideCylinderCore(pos[i], cfg, -0.25 * cfg.dp))
        {
            ++gate.cylinder_internal_fluid;
        }
        const Real ax = std::fabs(normal[i][0]);
        const Real ay = std::fabs(normal[i][1]);
        const Real az = std::fabs(normal[i][2]);
        if (az > ax && az > ay)
        {
            ++gate.z_face;
        }
        else if (ax >= ay && normal[i][0] < -0.5)
        {
            ++gate.inlet;
        }
        else if (ax >= ay && normal[i][0] > 0.5)
        {
            ++gate.outlet;
        }
        else if (ay > ax && std::fabs(normal[i][1]) > 0.5)
        {
            ++gate.top_bottom;
        }

        const size_t contact_neighbors = fluid_wall_contact.contact_configuration_[0][i].current_size_;
        if (contact_neighbors > 0)
        {
            ++gate.cylinder_contact_particles;
        }
        if (insideCylinderCore(pos[i], cfg, cfg.dp) && !insideCylinderCore(pos[i], cfg, -cfg.dp))
        {
            const size_t inner_neighbors = fluid_inner.inner_configuration_[i].current_size_;
            if (pos[i][2] <= z_band || pos[i][2] >= cfg.DW - z_band)
            {
                z_end_inner_sum += static_cast<Real>(inner_neighbors);
                z_end_sum += static_cast<Real>(contact_neighbors);
                ++z_end_count;
            }
            if (pos[i][2] >= mid_low && pos[i][2] <= mid_high)
            {
                z_mid_inner_sum += static_cast<Real>(inner_neighbors);
                z_mid_sum += static_cast<Real>(contact_neighbors);
                ++z_mid_count;
            }
        }
    }
    gate.z_end_inner_mean = z_end_count > 0 ? z_end_inner_sum / static_cast<Real>(z_end_count) : 0.0;
    gate.z_mid_inner_mean = z_mid_count > 0 ? z_mid_inner_sum / static_cast<Real>(z_mid_count) : 0.0;
    gate.z_end_contact_mean = z_end_count > 0 ? z_end_sum / static_cast<Real>(z_end_count) : 0.0;
    gate.z_mid_contact_mean = z_mid_count > 0 ? z_mid_sum / static_cast<Real>(z_mid_count) : 0.0;
    const bool z_contact_ok = gate.z_mid_contact_mean <= 0.0 ||
                              gate.z_end_contact_mean >= 0.8 * gate.z_mid_contact_mean;
    const bool wall_periodic_geometry_ok = gate.wall_z_min >= -TinyReal &&
                                           gate.wall_z_max <= cfg.DW + TinyReal &&
                                           gate.wall_max_abs_normal_z <= TinyReal;
    gate.pass = gate.cylinder_internal_fluid == 0 && gate.inlet > 0 && gate.outlet > 0 &&
                gate.top_bottom > 0 && gate.cylinder_contact_particles > 0 && z_contact_ok &&
                wall_periodic_geometry_ok;

    std::cout << "[Cylinder3D][GeometryGate] fluid=" << gate.fluid_particles
              << " wall=" << gate.wall_particles
              << " internal_fluid=" << gate.cylinder_internal_fluid
              << " inlet=" << gate.inlet
              << " outlet=" << gate.outlet
              << " top_bottom=" << gate.top_bottom
              << " z_face=" << gate.z_face
              << " contact_particles=" << gate.cylinder_contact_particles
              << " z_end_inner_mean=" << gate.z_end_inner_mean
              << " z_mid_inner_mean=" << gate.z_mid_inner_mean
              << " z_end_contact_mean=" << gate.z_end_contact_mean
              << " z_mid_contact_mean=" << gate.z_mid_contact_mean
              << " wall_z=[" << gate.wall_z_min << ", " << gate.wall_z_max << "]"
              << " wall_max_abs_normal_z=" << gate.wall_max_abs_normal_z
              << " pass=" << (gate.pass ? "yes" : "NO") << std::endl;
    return gate;
}

struct BoundaryGateResult
{
    size_t open_particles = 0;
    Real outlet_pressure_mean = 0.0;
    Real outlet_density_mean = 0.0;
    Real min_weight = 1.0;
    Real max_weight = 0.0;
    bool pass = false;
};

inline BoundaryGateResult reportBoundaryGate(BaseParticles &particles, const Cylinder3DConfig &cfg,
                                             const std::string &stage)
{
    BoundaryGateResult gate;
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Vecd *normal = particles.getVariableDataByName<Vecd>("NormalDirection");
    Real *pressure = particles.getVariableDataByName<Real>("Pressure");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *weight = particles.getVariableDataByName<Real>("InnerWeightSummation");
    const size_t total = particles.TotalRealParticles();
    Real outlet_p_sum = 0.0;
    Real outlet_rho_sum = 0.0;
    size_t outlet_count = 0;
    bool all_weight_ok = true;

    for (size_t i = 0; i != total; ++i)
    {
        const Real ax = std::fabs(normal[i][0]);
        const Real ay = std::fabs(normal[i][1]);
        const Real az = std::fabs(normal[i][2]);
        if (az > ax && az > ay)
        {
            continue;
        }
        const bool open_face = ax > 0.5 || ay > 0.5 || pos[i][0] < 0.0 || pos[i][0] > cfg.DL ||
                               pos[i][1] < 0.0 || pos[i][1] > cfg.DH;
        if (!open_face || insideCylinderCore(pos[i], cfg, cfg.boundary_n_layers * cfg.dp))
        {
            continue;
        }
        ++gate.open_particles;
        if (weight != nullptr)
        {
            gate.min_weight = SMIN(gate.min_weight, weight[i]);
            gate.max_weight = SMAX(gate.max_weight, weight[i]);
            all_weight_ok = all_weight_ok && isFiniteReal(weight[i]) &&
                            weight[i] >= -1.0e-12 && weight[i] <= 1.0 + 1.0e-12;
        }
        if (normal[i][0] > 0.5 || pos[i][0] > cfg.DL)
        {
            outlet_p_sum += pressure[i];
            outlet_rho_sum += rho[i];
            ++outlet_count;
        }
    }
    if (outlet_count > 0)
    {
        gate.outlet_pressure_mean = outlet_p_sum / static_cast<Real>(outlet_count);
        gate.outlet_density_mean = outlet_rho_sum / static_cast<Real>(outlet_count);
    }
    gate.pass = gate.open_particles > 0 && outlet_count > 0 &&
                isFiniteReal(gate.outlet_pressure_mean) &&
                isFiniteReal(gate.outlet_density_mean) && all_weight_ok;
    std::cout << "[Cylinder3D][BoundaryGate] " << stage
              << " open=" << gate.open_particles
              << " outlet_p_mean=" << gate.outlet_pressure_mean
              << " outlet_rho_mean=" << gate.outlet_density_mean
              << " weight=[" << gate.min_weight << ", " << gate.max_weight << "]"
              << " pass=" << (gate.pass ? "yes" : "NO") << std::endl;
    return gate;
}

struct PeriodicGateResult
{
    Real pressure_jump = 0.0;
    Real velocity_jump = 0.0;
    bool computed = false;
    bool pass = false;
};

inline PeriodicGateResult reportPeriodicGate(BaseParticles &particles, const Cylinder3DConfig &cfg,
                                             const std::string &stage)
{
    PeriodicGateResult gate;
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    Real *pressure = particles.getVariableDataByName<Real>("Pressure");
    const Real z_band = static_cast<Real>(cfg.boundary_n_layers) * cfg.dp;
    Vecd v_low = Vecd::Zero();
    Vecd v_high = Vecd::Zero();
    Real p_low = 0.0;
    Real p_high = 0.0;
    size_t n_low = 0;
    size_t n_high = 0;

    for (size_t i = 0; i != particles.TotalRealParticles(); ++i)
    {
        if (pos[i][2] <= z_band)
        {
            v_low += vel[i];
            p_low += pressure[i];
            ++n_low;
        }
        else if (pos[i][2] >= cfg.DW - z_band)
        {
            v_high += vel[i];
            p_high += pressure[i];
            ++n_high;
        }
    }
    if (n_low > 0 && n_high > 0)
    {
        v_low /= static_cast<Real>(n_low);
        v_high /= static_cast<Real>(n_high);
        p_low /= static_cast<Real>(n_low);
        p_high /= static_cast<Real>(n_high);
        gate.velocity_jump = (v_low - v_high).norm();
        gate.pressure_jump = std::fabs(p_low - p_high);
        gate.computed = true;
        gate.pass = isFiniteReal(gate.velocity_jump) && isFiniteReal(gate.pressure_jump);
    }
    std::cout << "[Cylinder3D][PeriodicGate] " << stage
              << " velocity_jump=" << gate.velocity_jump
              << " pressure_jump=" << gate.pressure_jump
              << " computed=" << (gate.computed ? "yes" : "NO")
              << " pass=" << (gate.pass ? "yes" : "NO") << std::endl;
    return gate;
}

} // namespace cylinder_3d
} // namespace SPH

#endif // CYLINDER_3D_GEOMETRY_HPP
