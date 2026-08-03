#ifndef CYLINDER_3D_COMPRESSIBLE_DIAGNOSTICS_HPP
#define CYLINDER_3D_COMPRESSIBLE_DIAGNOSTICS_HPP

/**
 * Diagnostics for the compressible cylinder case: wall pair direction contract,
 * thermodynamic positivity, wall mass/energy flux and load direction.
 *
 * This header only *verifies* the shared MUSCL/HLLC wall specialisation -- it
 * must not re-solve star states or provide a second flux implementation.
 * Failures are fail-fast; abs(area), clamping negative areas and silent
 * fallbacks are banned.
 *
 * Normal contract:
 *   n_wall  cylinder outward normal, solid -> fluid
 *   e_pair  fluid owner -> wall neighbour, gated by dot(e_pair, n_wall) > 0
 *   n_flux  fluid control-surface outward normal, n_flux = -e_pair
 */

#include "cylinder_3d_compressible_state.hpp"
#include "sphinxsys.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

namespace SPH
{
namespace cylinder_3d_compressible
{

/**
 * Dead zone for dot(e_pair, n_wall). e_pair and n_wall are both unit vectors, so
 * a purely tangential pair gives exactly 0 analytically and lands within a few
 * machine eps of it after the neighbour search. Treating that as a direction
 * failure would reject a large fraction of a curved wall's legitimate stencil.
 */
inline constexpr Real kPairAlignmentDeadZone = 1.0e-12;

struct PeriodicTopologyGateResult
{
    std::array<size_t, 2> fluid_x_seam_owners{{0, 0}};
    std::array<size_t, 4> fluid_xz_corner_owners{{0, 0, 0, 0}};
    std::array<size_t, 4> y_wall_x_seam_owners{{0, 0, 0, 0}};
    size_t nonfinite_inner_pairs = 0;
    size_t nonfinite_y_wall_pairs = 0;
    size_t inconsistent_inner_pairs = 0;
    size_t inconsistent_y_wall_pairs = 0;
    bool pass = false;
};

/**
 * Spatial A/B operator for the physical y-wall kernel-support band.
 *
 * The production gradients are retained everywhere by default. In the opt-in
 * diagnostic mode this zeros only the reconstruction gradients for fluid
 * owners that have a physical y wall in their support, so their MUSCL bridge
 * reduces to its cell-centred (first-order HLLC) state. It does not modify
 * conservative variables, wall states, flux areas, or the time step.
 */
class YWallFirstOrderReconstruction : public LocalDynamics
{
  public:
    YWallFirstOrderReconstruction(SPHBody &sph_body, const CompressibleCylinderConfig &cfg)
        : LocalDynamics(sph_body), cfg_(cfg),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          rho_grad_(particles_->getVariableDataByName<Vecd>("DensityGradient")),
          p_grad_(particles_->getVariableDataByName<Vecd>("PressureGradient")),
          vel_grad_(particles_->getVariableDataByName<Matd>("VelocityGradient")),
          mask_(particles_->registerStateVariableData<int>("FirstOrderReconstructionMask"))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        const Real distance_to_y_wall = std::min(pos_[index_i][1], cfg_.DH - pos_[index_i][1]);
        const bool selected = distance_to_y_wall <= cfg_.y_wall_first_order_band_dp * cfg_.dp;
        mask_[index_i] = selected ? 1 : 0;
        if (selected)
        {
            rho_grad_[index_i] = Vecd::Zero();
            vel_grad_[index_i] = Matd::Zero();
            p_grad_[index_i] = Vecd::Zero();
        }
    }

  private:
    const CompressibleCylinderConfig &cfg_;
    Vecd *pos_, *rho_grad_, *p_grad_;
    Matd *vel_grad_;
    int *mask_;
};

/**
 * First-order reconstruction only where the x-max open boundary intersects
 * either physical y wall. The x-min corners and y-wall mid-span remain at the
 * production MUSCL order, making this a direct ablation of the original first
 * bad owner's mixed-boundary neighbourhood.
 */
class OutletYWallCornerFirstOrderReconstruction : public LocalDynamics
{
  public:
    OutletYWallCornerFirstOrderReconstruction(
        SPHBody &sph_body, const CompressibleCylinderConfig &cfg)
        : LocalDynamics(sph_body), cfg_(cfg),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          rho_grad_(particles_->getVariableDataByName<Vecd>("DensityGradient")),
          p_grad_(particles_->getVariableDataByName<Vecd>("PressureGradient")),
          vel_grad_(particles_->getVariableDataByName<Matd>("VelocityGradient")),
          mask_(particles_->registerStateVariableData<int>("FirstOrderReconstructionMask"))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        const Real distance_to_y_wall = std::min(pos_[index_i][1], cfg_.DH - pos_[index_i][1]);
        const Real distance_to_x_outlet = cfg_.DL - pos_[index_i][0];
        const bool selected =
            distance_to_y_wall <= cfg_.y_wall_first_order_band_dp * cfg_.dp &&
            distance_to_x_outlet <= cfg_.outlet_corner_first_order_band_dp * cfg_.dp;
        mask_[index_i] = selected ? 1 : 0;
        if (selected)
        {
            rho_grad_[index_i] = Vecd::Zero();
            vel_grad_[index_i] = Matd::Zero();
            p_grad_[index_i] = Vecd::Zero();
        }
    }

  private:
    const CompressibleCylinderConfig &cfg_;
    Vecd *pos_, *rho_grad_, *p_grad_;
    Matd *vel_grad_;
    int *mask_;
};

/**
 * First-order reconstruction for the exact owner set used by wall-contact
 * interactions. Selection is topology-based: a fluid owner is marked when any
 * configured wall body's current contact neighborhood is non-empty.
 */
class WallContactFirstOrderReconstruction : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit WallContactFirstOrderReconstruction(BaseContactRelation &wall_contact_relation)
        : LocalDynamics(wall_contact_relation.getSPHBody()),
          DataDelegateContact(wall_contact_relation),
          rho_grad_(particles_->getVariableDataByName<Vecd>("DensityGradient")),
          p_grad_(particles_->getVariableDataByName<Vecd>("PressureGradient")),
          vel_grad_(particles_->getVariableDataByName<Matd>("VelocityGradient")),
          mask_(particles_->registerStateVariableData<int>("FirstOrderReconstructionMask"))
    {
        if (contact_configuration_.empty())
        {
            throw std::runtime_error(
                "WallContactFirstOrderReconstruction requires at least one wall contact body.");
        }
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        bool selected = false;
        for (ParticleConfiguration *configuration : contact_configuration_)
        {
            if ((*configuration)[index_i].current_size_ > 0)
            {
                selected = true;
                break;
            }
        }
        mask_[index_i] = selected ? 1 : 0;
        if (selected)
        {
            rho_grad_[index_i] = Vecd::Zero();
            vel_grad_[index_i] = Matd::Zero();
            p_grad_[index_i] = Vecd::Zero();
        }
    }

  private:
    Vecd *rho_grad_, *p_grad_;
    Matd *vel_grad_;
    int *mask_;
};

/**
 * Positivity floor (Zhang--Shu style), executed once per step after the
 * conservative second-half update.
 *
 * The A3 dissipation-limited HLLC removed the catastrophic Ma=2 failure
 * (wall through-flow + unchecked shock undershoot), but the second-order
 * MUSCL update can still drive isolated bow-shock cells to a slightly
 * negative EOS-recovered pressure (observed: a single cell at
 * p ~ -0.1*p_inf, bounded and non-propagating). This floor clamps rho and
 * the EOS-recovered pressure at factor x freestream and rewrites
 * TotalEnergy consistently in the same cell -- p and E must move together,
 * otherwise the state gate's recoveredPressure == p invariant would trip
 * on the clamped cell.
 *
 * Case-local by design: it touches only the seven registered state fields
 * through this body's own arrays, so the shared
 * EulerianCompressibleIntegration2ndHalfMUSCL update is unchanged for
 * every other MUSCL user. The floor energy shows up in the conservation
 * report's (ungated) state-minus-applied-rate imbalance, keeping the
 * correction visible instead of hidden inside the solver's rates.
 */
class Cylinder3DPositivityFloor : public LocalDynamics
{
  public:
    Cylinder3DPositivityFloor(SPHBody &sph_body, const CompressibleCylinderConfig &cfg)
        : LocalDynamics(sph_body), gamma_(cfg.gamma),
          rho_floor_(cfg.positivity_floor_factor * cfg.rho_inf),
          p_floor_(cfg.positivity_floor_factor * cfg.p_inf),
          rho_(particles_->getVariableDataByName<Real>("Density")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          mass_(particles_->getVariableDataByName<Real>("Mass")),
          p_(particles_->getVariableDataByName<Real>("Pressure")),
          E_(particles_->getVariableDataByName<Real>("TotalEnergy")),
          mom_(particles_->getVariableDataByName<Vecd>("Momentum"))
    {
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        // Degenerate-mass guard: mom/mass below assumes mass > 0. A cell that
        // lost all mass is reset to a static floor state (zero velocity,
        // floor pressure) instead of dividing by a non-positive mass.
        if (mass_[index_i] <= 0.0)
        {
            rho_[index_i] = rho_floor_;
            mass_[index_i] = rho_floor_ * Vol_[index_i];
            mom_[index_i] = Vecd::Zero();
            p_[index_i] = p_floor_;
            E_[index_i] = p_floor_ / (gamma_ - 1.0) * Vol_[index_i];
            return;
        }
        // ANY correction rewrites p and E together from the same post-floor
        // split: flooring rho alone would leave the solver-written p
        // inconsistent with the new mass, tripping the state gate's
        // recoveredPressure == p invariant even though p itself is positive.
        bool corrected = false;
        if (rho_[index_i] < rho_floor_)
        {
            rho_[index_i] = rho_floor_;
            mass_[index_i] = rho_floor_ * Vol_[index_i];
            corrected = true;
        }
        // The conservative velocity is the one the solver's own EOS
        // inversion uses (Momentum/Mass), not the stored lagging Velocity.
        const Vecd vel_conservative = mom_[index_i] / mass_[index_i];
        const Real kinetic_per_volume = 0.5 * rho_[index_i] * vel_conservative.squaredNorm();
        const Real rho_e = E_[index_i] / Vol_[index_i] - kinetic_per_volume;
        Real p_new = (gamma_ - 1.0) * rho_e;
        if (p_new < p_floor_)
        {
            p_new = p_floor_;
            corrected = true;
        }
        if (corrected)
        {
            // Keep the internal energy rho_e (possibly floored to the
            // pressure floor) and rewrite E against the post-floor kinetic
            // part, so Mass / TotalEnergy / Pressure stay EOS-consistent.
            p_[index_i] = p_new;
            E_[index_i] = (p_new / (gamma_ - 1.0) + kinetic_per_volume) * Vol_[index_i];
        }
    }

  private:
    const Real gamma_, rho_floor_, p_floor_;
    Real *rho_, *Vol_, *mass_, *p_, *E_;
    Vecd *mom_;
};

struct FirstOrderReconstructionMaskReport
{
    size_t marked = 0;
    size_t expected = 0;
    size_t lower_wall = 0;
    size_t upper_wall = 0;
    size_t mismatched = 0;
    size_t x_outlet_band_leak = 0;
    size_t y_interior_leak = 0;
    size_t cylinder_only = 0;
    size_t y_wall_only = 0;
    size_t multiple_wall_bodies = 0;
    size_t unmarked_contact = 0;
    size_t marked_without_contact = 0;
    bool pass = false;
};

inline FirstOrderReconstructionMaskReport reportFirstOrderReconstructionMask(
    BaseParticles &particles, const CompressibleCylinderConfig &cfg,
    BaseContactRelation &fluid_wall_contact, const std::string &stage)
{
    FirstOrderReconstructionMaskReport report;
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    int *mask = particles.getVariableDataByName<int>("FirstOrderReconstructionMask");
    const Real y_band = cfg.y_wall_first_order_band_dp * cfg.dp;
    const Real x_band = cfg.outlet_corner_first_order_band_dp * cfg.dp;
    const bool outlet_corner_mode = cfg.outlet_y_wall_corner_first_order_reconstruction;
    const bool wall_contact_mode = cfg.wall_contact_first_order_reconstruction;
    if (wall_contact_mode && fluid_wall_contact.contact_configuration_.empty())
    {
        throw std::runtime_error(
            "Wall-contact first-order mask gate requires at least one contact configuration.");
    }

    for (size_t i = 0; i != particles.TotalRealParticles(); ++i)
    {
        const Real distance_to_y_wall = std::min(pos[i][1], cfg.DH - pos[i][1]);
        const Real distance_to_x_outlet = cfg.DL - pos[i][0];
        const bool near_y_wall = distance_to_y_wall <= y_band;
        bool cylinder_contact = false;
        bool y_wall_contact = false;
        bool any_wall_contact = false;
        if (wall_contact_mode)
        {
            for (size_t k = 0; k != fluid_wall_contact.contact_configuration_.size(); ++k)
            {
                const bool has_contact =
                    fluid_wall_contact.contact_configuration_[k][i].current_size_ > 0;
                any_wall_contact = any_wall_contact || has_contact;
                cylinder_contact = cylinder_contact || (k == 0 && has_contact);
                y_wall_contact = y_wall_contact || (k > 0 && has_contact);
            }
        }
        const bool expected = wall_contact_mode
                                  ? any_wall_contact
                                  : near_y_wall &&
                                        (!outlet_corner_mode || distance_to_x_outlet <= x_band);
        const bool marked = mask[i] == 1;
        report.expected += expected ? 1 : 0;
        report.marked += marked ? 1 : 0;
        report.mismatched += marked != expected ? 1 : 0;
        report.unmarked_contact += wall_contact_mode && expected && !marked ? 1 : 0;
        report.marked_without_contact += wall_contact_mode && marked && !expected ? 1 : 0;
        if (wall_contact_mode && expected)
        {
            report.cylinder_only += cylinder_contact && !y_wall_contact ? 1 : 0;
            report.y_wall_only += y_wall_contact && !cylinder_contact ? 1 : 0;
            report.multiple_wall_bodies += cylinder_contact && y_wall_contact ? 1 : 0;
        }
        if (marked)
        {
            if (pos[i][1] < 0.5 * cfg.DH)
            {
                ++report.lower_wall;
            }
            else
            {
                ++report.upper_wall;
            }
            report.y_interior_leak += !wall_contact_mode && !near_y_wall ? 1 : 0;
            report.x_outlet_band_leak +=
                !wall_contact_mode && outlet_corner_mode && distance_to_x_outlet > x_band ? 1 : 0;
        }
    }

    if (wall_contact_mode)
    {
        const size_t cylinder_contact_owners = report.cylinder_only + report.multiple_wall_bodies;
        const size_t y_wall_contact_owners = report.y_wall_only + report.multiple_wall_bodies;
        report.pass = report.marked == report.expected && report.marked > 0 &&
                      report.mismatched == 0 && report.unmarked_contact == 0 &&
                      report.marked_without_contact == 0 && cylinder_contact_owners > 0 &&
                      (cfg.y_boundary_mode != YBoundaryMode::Wall || y_wall_contact_owners > 0);
    }
    else
    {
        report.pass = report.marked == report.expected && report.marked > 0 &&
                      report.lower_wall > 0 && report.upper_wall > 0 &&
                      report.mismatched == 0 && report.x_outlet_band_leak == 0 &&
                      report.y_interior_leak == 0;
    }
    std::cout << "[Cylinder3DCompressible][FirstOrderMaskGate] " << stage
              << " mode=" << (wall_contact_mode
                                    ? "wall-contact-owners"
                                    : (outlet_corner_mode ? "x-max-y-wall-corners" : "full-y-wall-band"))
              << " marked=" << report.marked
              << " expected=" << report.expected
              << " lower=" << report.lower_wall
              << " upper=" << report.upper_wall
              << " mismatched=" << report.mismatched
              << " x_outlet_band_leak=" << report.x_outlet_band_leak
              << " y_interior_leak=" << report.y_interior_leak
              << " cylinder_only=" << report.cylinder_only
              << " y_wall_only=" << report.y_wall_only
              << " multiple_wall_bodies=" << report.multiple_wall_bodies
              << " unmarked_contact=" << report.unmarked_contact
              << " marked_without_contact=" << report.marked_without_contact
              << " pass=" << (report.pass ? "yes" : "NO") << std::endl;
    return report;
}

/**
 * Verify the periodic links actually used by this case-local topology.
 *
 * Cell-linked-list periodicity inserts translated entries but preserves the real
 * particle index. The raw positions therefore identify which pair crossed a
 * seam, while r/e/dW verify the translated neighbour data passed to the HLLC
 * reconstruction. The four x/z corners catch a missing diagonal insertion;
 * the four y-wall x-seam bands independently prove fluid-wall contact closure.
 */
inline PeriodicTopologyGateResult reportPeriodicTopologyGate(
    FluidBody &fluid_body, BaseInnerRelation &fluid_inner,
    ContactRelation &fluid_wall_contact, BaseParticles &y_wall_particles,
    const CompressibleCylinderConfig &cfg, const std::string &stage)
{
    if (cfg.x_boundary_mode != XBoundaryMode::Periodic ||
        cfg.y_boundary_mode != YBoundaryMode::Wall)
    {
        throw std::runtime_error("reportPeriodicTopologyGate requires x=periodic and y=wall.");
    }
    if (fluid_wall_contact.contact_configuration_.size() < 2)
    {
        throw std::runtime_error("reportPeriodicTopologyGate requires the y-wall contact configuration.");
    }

    PeriodicTopologyGateResult gate;
    Vecd *fluid_pos = fluid_body.getBaseParticles().getVariableDataByName<Vecd>("Position");
    Vecd *y_wall_pos = y_wall_particles.getVariableDataByName<Vecd>("Position");
    const Real support = fluid_body.getSPHAdaptation().getKernel()->CutOffRadius();
    const size_t total = fluid_body.getBaseParticles().TotalRealParticles();

    const auto seam_side = [](Real coordinate, Real length, Real band)
    {
        return coordinate < band ? 0 : (coordinate > length - band ? 1 : -1);
    };
    const auto opposite = [](int side)
    { return side == 0 ? 1 : 0; };
    const auto minimum_image_displacement = [&cfg](const Vecd &pos_i, const Vecd &pos_j)
    {
        Vecd displacement = pos_i - pos_j;
        const auto minimum_image_component = [](Real delta, Real period)
        {
            if (delta > 0.5 * period)
            {
                return delta - period;
            }
            if (delta < -0.5 * period)
            {
                return delta + period;
            }
            return delta;
        };
        displacement[0] = minimum_image_component(displacement[0], cfg.DL);
        displacement[2] = minimum_image_component(displacement[2], cfg.DW);
        return displacement;
    };
    const auto has_consistent_periodic_geometry = [&minimum_image_displacement, support](
                                                      const Vecd &pos_i, const Vecd &pos_j,
                                                      const Neighborhood &neighborhood, size_t n)
    {
        const Vecd expected_displacement = minimum_image_displacement(pos_i, pos_j);
        const Real expected_distance = expected_displacement.norm();
        const Real tolerance = 1.0e-10 * SMAX(Real(1.0), expected_distance);
        if (!(expected_distance > TinyReal && expected_distance <= support + tolerance) ||
            std::fabs(neighborhood.r_ij_[n] - expected_distance) > tolerance)
        {
            return false;
        }
        const Vecd expected_direction = expected_displacement / (expected_distance + TinyReal);
        return (neighborhood.e_ij_[n] - expected_direction).norm() <= 1.0e-10;
    };

    for (size_t i = 0; i != total; ++i)
    {
        const int x_side = seam_side(fluid_pos[i][0], cfg.DL, support);
        const int z_side = seam_side(fluid_pos[i][2], cfg.DW, support);
        bool x_cross_seen = false;
        bool xz_cross_seen = false;
        if (x_side >= 0)
        {
            Neighborhood &neighborhood = fluid_inner.inner_configuration_[i];
            for (size_t n = 0; n != neighborhood.current_size_; ++n)
            {
                const size_t j = neighborhood.j_[n];
                if (j >= total || !isFiniteRealValue(neighborhood.r_ij_[n]) ||
                    !isFiniteRealValue(neighborhood.dW_ij_[n]) ||
                    !isFiniteVecValue(neighborhood.e_ij_[n]))
                {
                    ++gate.nonfinite_inner_pairs;
                    continue;
                }
                const int x_neighbor_side = seam_side(fluid_pos[j][0], cfg.DL, support);
                const bool crosses_x = x_neighbor_side == opposite(x_side);
                if (crosses_x &&
                    !has_consistent_periodic_geometry(fluid_pos[i], fluid_pos[j], neighborhood, n))
                {
                    ++gate.inconsistent_inner_pairs;
                    continue;
                }
                x_cross_seen = x_cross_seen || crosses_x;
                if (z_side >= 0)
                {
                    const int z_neighbor_side = seam_side(fluid_pos[j][2], cfg.DW, support);
                    xz_cross_seen = xz_cross_seen ||
                                    (crosses_x && z_neighbor_side == opposite(z_side));
                }
            }
            if (x_cross_seen)
            {
                ++gate.fluid_x_seam_owners[static_cast<size_t>(x_side)];
            }
            if (z_side >= 0 && xz_cross_seen)
            {
                ++gate.fluid_xz_corner_owners[static_cast<size_t>(2 * x_side + z_side)];
            }
        }

        const int y_side = seam_side(fluid_pos[i][1], cfg.DH, support);
        if (x_side < 0 || y_side < 0)
        {
            continue;
        }
        bool y_wall_cross_seen = false;
        Neighborhood &y_wall_neighborhood = fluid_wall_contact.contact_configuration_[1][i];
        for (size_t n = 0; n != y_wall_neighborhood.current_size_; ++n)
        {
            const size_t j = y_wall_neighborhood.j_[n];
            if (j >= y_wall_particles.TotalRealParticles() ||
                !isFiniteRealValue(y_wall_neighborhood.r_ij_[n]) ||
                !isFiniteRealValue(y_wall_neighborhood.dW_ij_[n]) ||
                !isFiniteVecValue(y_wall_neighborhood.e_ij_[n]))
            {
                ++gate.nonfinite_y_wall_pairs;
                continue;
            }
            const int x_neighbor_side = seam_side(y_wall_pos[j][0], cfg.DL, support);
            const bool crosses_x = x_neighbor_side == opposite(x_side);
            if (crosses_x && !has_consistent_periodic_geometry(fluid_pos[i], y_wall_pos[j],
                                                               y_wall_neighborhood, n))
            {
                ++gate.inconsistent_y_wall_pairs;
                continue;
            }
            y_wall_cross_seen = y_wall_cross_seen || crosses_x;
        }
        if (y_wall_cross_seen)
        {
            ++gate.y_wall_x_seam_owners[static_cast<size_t>(2 * y_side + x_side)];
        }
    }

    const auto all_nonzero = [](const auto &counts)
    {
        for (size_t count : counts)
        {
            if (count == 0)
            {
                return false;
            }
        }
        return true;
    };
    gate.pass = all_nonzero(gate.fluid_x_seam_owners) &&
                all_nonzero(gate.fluid_xz_corner_owners) &&
                all_nonzero(gate.y_wall_x_seam_owners) &&
                gate.nonfinite_inner_pairs == 0 && gate.nonfinite_y_wall_pairs == 0 &&
                gate.inconsistent_inner_pairs == 0 && gate.inconsistent_y_wall_pairs == 0;

    std::cout << "[Cylinder3DCompressible][PeriodicTopologyGate] " << stage
              << " fluid_x=[" << gate.fluid_x_seam_owners[0] << ", "
              << gate.fluid_x_seam_owners[1] << "]"
              << " fluid_xz=[" << gate.fluid_xz_corner_owners[0] << ", "
              << gate.fluid_xz_corner_owners[1] << ", "
              << gate.fluid_xz_corner_owners[2] << ", "
              << gate.fluid_xz_corner_owners[3] << "]"
              << " y_wall_x=[" << gate.y_wall_x_seam_owners[0] << ", "
              << gate.y_wall_x_seam_owners[1] << ", "
              << gate.y_wall_x_seam_owners[2] << ", "
              << gate.y_wall_x_seam_owners[3] << "]"
              << " nonfinite_inner=" << gate.nonfinite_inner_pairs
              << " nonfinite_y_wall=" << gate.nonfinite_y_wall_pairs
              << " inconsistent_inner=" << gate.inconsistent_inner_pairs
              << " inconsistent_y_wall=" << gate.inconsistent_y_wall_pairs
              << " pass=" << (gate.pass ? "yes" : "NO") << std::endl;
    return gate;
}

struct BulkFlowReport
{
    Real total_mass = 0.0;
    Real total_x_momentum = 0.0;
    Real u_bulk = 0.0;
    Real p_min = std::numeric_limits<Real>::infinity();
    size_t p_min_index = 0;
    Real max_mach = 0.0;
    bool pass = false;
};

/**
 * Periodic-channel diagnostic only: the conservative velocity Momentum/Mass is
 * mass-weighted to avoid the stored Velocity's intentional second-half lag.
 * No forcing is added in this first boundary-isolation experiment, so U_bulk
 * distinguishes a stable topology from a wake that merely decayed away.
 */
inline BulkFlowReport reportBulkFlow(const CompressibleConservativeStateView &state,
                                     BaseParticles &particles, const std::string &stage)
{
    BulkFlowReport report;
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t total = particles.TotalRealParticles();
    for (size_t i = 0; i != total; ++i)
    {
        const Real mass = state.mass[i];
        const Real rho = state.rho[i];
        const Real pressure = state.p[i];
        const Vecd &momentum = state.mom[i];
        if (!isFiniteRealValue(mass) || !isFiniteRealValue(rho) ||
            !isFiniteRealValue(pressure) || !isFiniteVecValue(momentum) ||
            mass <= 0.0 || rho <= 0.0 || pressure <= 0.0)
        {
            std::cout << "[Cylinder3DCompressible][BulkFlow] " << stage
                      << " invalid_state_i=" << i
                      << " pos=(" << pos[i][0] << ", " << pos[i][1] << ", " << pos[i][2] << ")"
                      << " mass=" << mass
                      << " rho=" << rho
                      << " p=" << pressure
                      << " E=" << state.E[i]
                      << " mom=(" << momentum[0] << ", " << momentum[1] << ", "
                      << momentum[2] << ")" << std::endl;
            return report;
        }
        report.total_mass += mass;
        report.total_x_momentum += momentum[0];
        if (pressure < report.p_min)
        {
            report.p_min = pressure;
            report.p_min_index = i;
        }
        const Real local_sound_speed = std::sqrt(state.gamma * pressure / rho);
        report.max_mach = SMAX(report.max_mach, (momentum / mass).norm() / local_sound_speed);
    }
    report.pass = report.total_mass > 0.0 && isFiniteRealValue(report.total_mass) &&
                  isFiniteRealValue(report.total_x_momentum) && isFiniteRealValue(report.max_mach);
    if (report.pass)
    {
        report.u_bulk = report.total_x_momentum / report.total_mass;
    }
    std::cout << "[Cylinder3DCompressible][BulkFlow] " << stage
              << " mass=" << report.total_mass
              << " px=" << report.total_x_momentum
              << " u_bulk=" << report.u_bulk
              << " p_min=" << report.p_min
              << " p_min_i=" << report.p_min_index;
    if (pos != nullptr && report.p_min_index < total)
    {
        std::cout << " p_min_pos=(" << pos[report.p_min_index][0] << ", "
                  << pos[report.p_min_index][1] << ", " << pos[report.p_min_index][2] << ")";
    }
    std::cout << " max_mach=" << report.max_mach
              << " pass=" << (report.pass ? "yes" : "NO") << std::endl;
    return report;
}

struct WallPairGateResult
{
    size_t pairs = 0;
    size_t owners_with_pairs = 0;
    size_t tangential_pairs = 0;
    size_t direction_failures = 0;
    size_t nonfinite_weights = 0;
    size_t near_zero_normals = 0;
    size_t nonpositive_effective_area = 0;
    Real min_dot_e_pair_n_wall = std::numeric_limits<Real>::max();
    Real min_pair_distance = std::numeric_limits<Real>::infinity();
    bool pass = false;
};

/**
 * Pair-level gate over the fluid -> cylinder contact configuration.
 *
 * The outward-oriented effective area of a pair is -dW_ij * Vol_j (dW_ij is
 * negative by SPH convention); a non-positive value means the discrete surface
 * element points the wrong way and is reported, never repaired.
 */
inline WallPairGateResult reportWallPairGate(BaseParticles &fluid_particles,
                                             BaseParticles &wall_particles,
                                             ContactRelation &fluid_wall_contact,
                                             const std::string &stage,
                                             size_t wall_contact_index = 0)
{
    WallPairGateResult gate;
    if (wall_contact_index >= fluid_wall_contact.contact_configuration_.size())
    {
        throw std::runtime_error("reportWallPairGate: wall contact index is out of range.");
    }
    Vecd *pos = fluid_particles.getVariableDataByName<Vecd>("Position");
    Vecd *wall_normal = wall_particles.getVariableDataByName<Vecd>("NormalDirection");
    Real *wall_Vol = wall_particles.getVariableDataByName<Real>("VolumetricMeasure");

    for (size_t i = 0; i != fluid_particles.TotalRealParticles(); ++i)
    {
        Neighborhood &neighborhood = fluid_wall_contact.contact_configuration_[wall_contact_index][i];
        if (neighborhood.current_size_ > 0)
        {
            ++gate.owners_with_pairs;
        }
        for (size_t n = 0; n != neighborhood.current_size_; ++n)
        {
            ++gate.pairs;
            const size_t index_j = neighborhood.j_[n];
            const Vecd &e_pair = neighborhood.e_ij_[n];
            const Vecd &n_wall = wall_normal[index_j];

            if (!isFiniteVecValue(e_pair) || !isFiniteVecValue(n_wall) ||
                !isFiniteRealValue(neighborhood.dW_ij_[n]) ||
                !isFiniteRealValue(neighborhood.r_ij_[n]) ||
                !isFiniteRealValue(wall_Vol[index_j]))
            {
                ++gate.nonfinite_weights;
                continue;
            }
            if (n_wall.norm() <= 1.0e-12)
            {
                ++gate.near_zero_normals;
                continue;
            }

            const Real alignment = e_pair.dot(n_wall);
            gate.min_dot_e_pair_n_wall = SMIN(gate.min_dot_e_pair_n_wall, alignment);
            if (std::fabs(alignment) <= kPairAlignmentDeadZone)
            {
                // Purely tangential pair: dot == 0 in exact arithmetic, within a
                // few eps of it in practice. Both vectors are unit length, so the
                // dead zone is absolute. Counted, not failed.
                ++gate.tangential_pairs;
            }
            else if (alignment < 0.0)
            {
                if (gate.direction_failures == 0)
                {
                    std::cout << "[Cylinder3DCompressible][WallPair] first direction failure"
                              << " stage=" << stage
                              << " fluid_i=" << i
                              << " pos=(" << pos[i][0] << ", " << pos[i][1] << ", " << pos[i][2] << ")"
                              << " wall_j=" << index_j
                              << " e_pair=(" << e_pair[0] << ", " << e_pair[1] << ", " << e_pair[2] << ")"
                              << " n_wall=(" << n_wall[0] << ", " << n_wall[1] << ", " << n_wall[2] << ")"
                              << " dot=" << alignment << std::endl;
                }
                ++gate.direction_failures;
            }

            const Real effective_area = -neighborhood.dW_ij_[n] * wall_Vol[index_j];
            if (effective_area <= 0.0)
            {
                ++gate.nonpositive_effective_area;
            }
            gate.min_pair_distance = SMIN(gate.min_pair_distance, neighborhood.r_ij_[n]);
        }
    }

    gate.pass = gate.pairs > 0 && gate.owners_with_pairs > 0 &&
                gate.direction_failures == 0 && gate.nonfinite_weights == 0 &&
                gate.near_zero_normals == 0 && gate.nonpositive_effective_area == 0;

    std::cout << "[Cylinder3DCompressible][WallPairGate] " << stage
              << " pairs=" << gate.pairs
              << " owners=" << gate.owners_with_pairs
              << " tangential=" << gate.tangential_pairs
              << " direction_failures=" << gate.direction_failures
              << " nonfinite=" << gate.nonfinite_weights
              << " near_zero_normals=" << gate.near_zero_normals
              << " nonpositive_area=" << gate.nonpositive_effective_area
              << " min_dot=" << gate.min_dot_e_pair_n_wall
              << " min_r=" << gate.min_pair_distance
              << " pass=" << (gate.pass ? "yes" : "NO") << std::endl;
    return gate;
}

struct WallFluxBalanceResult
{
    size_t wall_owners = 0;
    Real net_mass_flux = 0.0;
    Real net_energy_flux = 0.0;
    /** Sum of |per-owner rate|, i.e. the scale the net is compared against. */
    Real gross_mass_flux = 0.0;
    Real gross_energy_flux = 0.0;
    Real relative_mass_imbalance = 0.0;
    Real relative_energy_imbalance = 0.0;
    /** True when the net flux is already tight in absolute terms. */
    bool absolute_floor_met = false;
    bool pass = false;
};

/**
 * Relative mass / energy imbalance of the cylinder wall flux.
 *
 * PRECONDITION: the caller must have zeroed MassChangeRate / TotalEnergyChangeRate
 * and then executed the *wall-contact-only* second half step, so the rates hold
 * the wall flux alone. Running this after the combined inner+wall pass would
 * measure the interior flux too and the judgement would be meaningless.
 *
 * This is NOT a zero-penetration proof. The plan's absolute "net penetration
 * < 1e-12" criterion applies to a flat wall under normal incidence, where the
 * mirrored ghost state makes the normal mass flux vanish exactly. On a
 * lattice-generated curved cylinder at D/dp = 5 the discrete surface is not
 * symmetric, so each owner carries an O(dp) tangential residual and the signed
 * sum is O(dp), not O(1e-12). What remains meaningful is that the residual stays
 * small *relative to the gross flux through the same surface*: a genuinely
 * leaking or mis-oriented wall gives an imbalance of order one.
 *
 * The absolute flat-wall contract belongs in a dedicated unit-test fixture, not
 * on this geometry.
 */
inline WallFluxBalanceResult reportWallFluxBalance(BaseParticles &fluid_particles,
                                                   ContactRelation &fluid_wall_contact,
                                                   const std::string &stage,
                                                   Real relative_tolerance,
                                                   Real absolute_floor)
{
    WallFluxBalanceResult gate;
    Real *dmass_dt = fluid_particles.getVariableDataByName<Real>("MassChangeRate");
    Real *dE_dt = fluid_particles.getVariableDataByName<Real>("TotalEnergyChangeRate");
    if (dmass_dt == nullptr || dE_dt == nullptr)
    {
        throw std::runtime_error("reportWallFluxBalance: MassChangeRate / TotalEnergyChangeRate not registered.");
    }

    for (size_t i = 0; i != fluid_particles.TotalRealParticles(); ++i)
    {
        bool has_wall_contact = false;
        for (const auto &configuration : fluid_wall_contact.contact_configuration_)
        {
            has_wall_contact = has_wall_contact || configuration[i].current_size_ > 0;
        }
        if (!has_wall_contact)
        {
            continue;
        }
        ++gate.wall_owners;
        if (!isFiniteRealValue(dmass_dt[i]) || !isFiniteRealValue(dE_dt[i]))
        {
            throw std::runtime_error("reportWallFluxBalance: non-finite wall-adjacent change rate.");
        }
        gate.net_mass_flux += dmass_dt[i];
        gate.net_energy_flux += dE_dt[i];
        gate.gross_mass_flux += std::fabs(dmass_dt[i]);
        gate.gross_energy_flux += std::fabs(dE_dt[i]);
    }

    // A wall that produces no flux at all would give 0/0; treat it as a failure
    // rather than a vacuous pass, since the stencil is then not exercised.
    const bool flux_exercised = gate.gross_mass_flux > 0.0 && gate.gross_energy_flux > 0.0;
    if (flux_exercised)
    {
        gate.relative_mass_imbalance = std::fabs(gate.net_mass_flux) / gate.gross_mass_flux;
        gate.relative_energy_imbalance = std::fabs(gate.net_energy_flux) / gate.gross_energy_flux;
    }
    // Absolute floor. On an exactly uniform freestream the mirrored ghost state
    // (vel_reflect = 2*vel_wall_ave - vel_fluid) cancels the wall flux down to
    // machine precision, so BOTH net and gross land around 1e-18 and their ratio
    // is pure rounding noise -- a relative judgement is meaningless there. Below
    // the floor the wall is provably tight in absolute terms, which is the
    // stronger statement; above it the relative criterion takes over and stays
    // meaningful once the wake makes the flux O(1).
    //
    // NOTE on the premise: that machine-precision cancellation requires the
    // initial state to be uniform. With velocity_noise_amplitude > 0 the per-
    // particle scatter means vel_reflect no longer cancels at t = 0, the floor is
    // not met, and the verdict falls through to the relative criterion -- which is
    // the branch that carries the O(1) wake case anyway, but was calibrated on the
    // unperturbed state. Both quantities are printed below, so which branch
    // decided the verdict is visible in the log rather than implicit.
    gate.absolute_floor_met = std::fabs(gate.net_mass_flux) <= absolute_floor &&
                              std::fabs(gate.net_energy_flux) <= absolute_floor;
    gate.pass = gate.wall_owners > 0 && flux_exercised &&
                (gate.absolute_floor_met ||
                 (gate.relative_mass_imbalance <= relative_tolerance &&
                  gate.relative_energy_imbalance <= relative_tolerance));

    std::cout << "[Cylinder3DCompressible][WallFluxBalanceGate] " << stage
              << " owners=" << gate.wall_owners
              << " net_mass=" << gate.net_mass_flux
              << " gross_mass=" << gate.gross_mass_flux
              << " rel_mass=" << gate.relative_mass_imbalance
              << " net_energy=" << gate.net_energy_flux
              << " gross_energy=" << gate.gross_energy_flux
              << " rel_energy=" << gate.relative_energy_imbalance
              << " (rel tol " << relative_tolerance
              << ", abs floor " << absolute_floor << ")"
              << " abs_floor_met=" << (gate.absolute_floor_met ? "yes" : "no")
              << " pass=" << (gate.pass ? "yes" : "NO") << std::endl;
    return gate;
}

struct CylinderLoadReport
{
    Vecd pressure_force = Vecd::Zero();
    Vecd viscous_force = Vecd::Zero();
    Vecd total_force = Vecd::Zero();
    Real cd = 0.0;
    Real cl = 0.0;
    bool finite = false;
    /** Drag must oppose the streamwise freestream, i.e. Cd > 0. */
    bool drag_direction_ok = false;
};

/**
 * Summed cylinder load from the existing solid_dynamics load fields.
 *
 * Pressure and viscous contributions stay separated, and Cd/Cl are normalised by
 * 0.5*rho_inf*u_inf^2*D*DW. Nothing here re-derives the surface traction.
 */
inline CylinderLoadReport reportCylinderLoad(BaseParticles &wall_particles,
                                             const CompressibleCylinderConfig &cfg,
                                             const std::string &stage)
{
    CylinderLoadReport report;
    Vecd *pressure_force = wall_particles.getVariableDataByName<Vecd>("PressureForceFromFluid");
    Vecd *viscous_force = wall_particles.getVariableDataByName<Vecd>("ViscousForceFromFluid");
    if (pressure_force == nullptr || viscous_force == nullptr)
    {
        throw std::runtime_error("reportCylinderLoad: cylinder load fields are not registered.");
    }

    for (size_t i = 0; i != wall_particles.TotalRealParticles(); ++i)
    {
        report.pressure_force += pressure_force[i];
        report.viscous_force += viscous_force[i];
    }
    report.total_force = report.pressure_force + report.viscous_force;
    report.finite = isFiniteVecValue(report.pressure_force) && isFiniteVecValue(report.viscous_force);

    const Real scale = forceScale(cfg);
    if (scale > 0.0 && report.finite)
    {
        report.cd = report.total_force[0] / scale;
        report.cl = report.total_force[1] / scale;
        // PressureForceFromFluid accumulates the force ON the solid, so a body in
        // a +x freestream must see a positive streamwise load. A negative Cd
        // means the surface traction sign or the wall normal is inverted -- not a
        // physically small drag.
        report.drag_direction_ok = report.cd > 0.0;
    }

    std::cout << "[Cylinder3DCompressible][Load] " << stage
              << " Fp=(" << report.pressure_force[0] << ", " << report.pressure_force[1]
              << ", " << report.pressure_force[2] << ")"
              << " Fv=(" << report.viscous_force[0] << ", " << report.viscous_force[1]
              << ", " << report.viscous_force[2] << ")"
              << " Cd=" << report.cd
              << " Cl=" << report.cl
              << " finite=" << (report.finite ? "yes" : "NO")
              << " drag_direction=" << (report.drag_direction_ok ? "ok" : "WRONG") << std::endl;
    return report;
}

struct ConservationDriftReport
{
    /** Raw change of the totals. Reported only: an open domain may change. */
    Real mass_change = 0.0;
    Real energy_change = 0.0;
    /** Budget residual: change minus the flux that crossed the boundary. */
    Real mass_imbalance = 0.0;
    Real energy_imbalance = 0.0;
    bool pass = false;
};

/**
 * Bookkeeping report for an OPEN domain. REPORTED, NOT JUDGED -- read on.
 *
 * With four open x/y faces the totals are physically free to move: the dp=0.04
 * run gains mass while the startup compression wave is convected out and the
 * inlet keeps feeding rho_inf*u_inf in. So "sum(Mass) stays constant" is the
 * wrong statement -- it only holds for a closed domain, and judging it reports a
 * growing "drift" (1.9e-3 at 600 steps) on a scheme that is in fact conservative.
 *
 * The obvious replacement -- require the change of the totals to match the
 * accumulated net flux -- is NOT a test either, and this is worth being explicit
 * about because it looks like one. The shared solver advances the state with
 *     mass_[i] += dmass_dt_[i] * dt        (2nd half update)
 * and accumulateAppliedRates() re-sums THE SAME dmass_dt_ array over THE SAME index
 * range with THE SAME dt. Hence
 *     sum(Mass) - sum(Mass)_ref  ==  integral(net flux) dt
 * identically, to summation round-off, for a conservative scheme, a leaking one,
 * and one with a sign error alike: mass created inside dmass_dt_ appears on both
 * sides of the subtraction. The residual only verifies that the diagnostic can
 * re-add the solver's own increments.
 *
 * A real conservation test needs the boundary flux computed INDEPENDENTLY of the
 * state update -- a surface integral over the ghost pairs, i.e. its own stencil.
 * Until that exists, both numbers are printed and neither is gated, so this
 * cannot report a false pass. pass is therefore finiteness only.
 */
inline ConservationDriftReport reportConservationDrift(const ConservationBudget &reference,
                                                       const ConservationBudget &current,
                                                       const CompressibleCylinderConfig &cfg,
                                                       const std::string &stage)
{
    ConservationDriftReport report;
    if (reference.total_mass <= 0.0 || reference.total_energy <= 0.0)
    {
        throw std::runtime_error("reportConservationDrift: reference budget must be positive.");
    }
    const Real mass_delta = current.total_mass - reference.total_mass;
    const Real energy_delta = current.total_energy - reference.total_energy;
    report.mass_change = std::fabs(mass_delta) / reference.total_mass;
    report.energy_change = std::fabs(energy_delta) / reference.total_energy;
    report.mass_imbalance =
        std::fabs(mass_delta - current.accumulated_mass_rate) / reference.total_mass;
    report.energy_imbalance =
        std::fabs(energy_delta - current.accumulated_energy_rate) / reference.total_energy;
    // Finiteness only. The imbalance terms are identities (see the header
    // comment) and the raw changes are legitimately non-zero on an open domain,
    // so neither may enter the verdict. A NaN total, however, is a real defect.
    report.pass = isFiniteRealValue(report.mass_change) && isFiniteRealValue(report.energy_change) &&
                  isFiniteRealValue(current.total_mass) && isFiniteRealValue(current.total_energy) &&
                  current.total_mass > 0.0 && current.total_energy > 0.0;

    std::cout << "[Cylinder3DCompressible][ConservationReport] " << stage
              << " mass_change=" << report.mass_change
              << " energy_change=" << report.energy_change
              << " mass=" << current.total_mass
              << " energy=" << current.total_energy
              << " | not a conservation test:"
              << " state_minus_accumulated_rate_mass=" << report.mass_imbalance
              << " energy=" << report.energy_imbalance
              << " (identically zero by construction; needs an independent"
              << " boundary surface integral to have content)"
              << " finite_and_positive=" << (report.pass ? "yes" : "NO") << std::endl;
    return report;
}

} // namespace cylinder_3d_compressible
} // namespace SPH

#endif // CYLINDER_3D_COMPRESSIBLE_DIAGNOSTICS_HPP
