#ifndef CYLINDER_3D_COMPRESSIBLE_BOUNDARY_HPP
#define CYLINDER_3D_COMPRESSIBLE_BOUNDARY_HPP

/**
 * Fully compressible far-field open boundary built on the shared Eulerian ghost
 * framework (GhostCreationInESPH / GhostBoundaryConditionSetupInESPH /
 * GhostKernelGradientUpdate).
 *
 * The x faces are cfg-selectable: Mach-aware far field (the default) or periodic.
 * The y faces are cfg-selectable: Mach-aware far field (the default) or physical
 * static walls. Open owners use
 * BoundaryType = 9 and a single applyFarFieldBoundary() branch, mirroring the
 * existing 2D fully compressible example. Solid-contact faces are masked before
 * GhostCreationInESPH is constructed: that constructor immediately creates the
 * ghosts, so post-filtering is impossible.
 */

#include "cylinder_3d_compressible_geometry.hpp"
#include "eulerian_ghost_boundary.h"
#include "sphinxsys.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace SPH
{
namespace cylinder_3d_compressible
{

/** The four non-periodic open faces, in fixed corner tie-break priority order. */
enum class OpenFace
{
    XMin = 0,
    XMax = 1,
    YMin = 2,
    YMax = 3,
    None = 4
};

inline const char *openFaceName(OpenFace face)
{
    switch (face)
    {
    case OpenFace::XMin:
        return "x-min";
    case OpenFace::XMax:
        return "x-max";
    case OpenFace::YMin:
        return "y-min";
    case OpenFace::YMax:
        return "y-max";
    default:
        return "none";
    }
}

inline Vecd openFaceOutwardNormal(OpenFace face)
{
    switch (face)
    {
    case OpenFace::XMin:
        return Vecd(-1.0, 0.0, 0.0);
    case OpenFace::XMax:
        return Vecd(1.0, 0.0, 0.0);
    case OpenFace::YMin:
        return Vecd(0.0, -1.0, 0.0);
    case OpenFace::YMax:
        return Vecd(0.0, 1.0, 0.0);
    default:
        return Vecd::Zero();
    }
}

inline bool isActiveOpenFace(OpenFace face, const CompressibleCylinderConfig &cfg)
{
    return (cfg.x_boundary_mode == XBoundaryMode::FarField &&
            (face == OpenFace::XMin || face == OpenFace::XMax)) ||
           (cfg.y_boundary_mode == YBoundaryMode::FarField &&
            (face == OpenFace::YMin || face == OpenFace::YMax));
}

inline bool hasActiveOpenFace(const CompressibleCylinderConfig &cfg)
{
    return cfg.x_boundary_mode == XBoundaryMode::FarField ||
           cfg.y_boundary_mode == YBoundaryMode::FarField;
}

/**
 * Classify a position onto exactly one open face: the *nearest* one.
 *
 * Picking by nearest distance rather than by a fixed scan order matters at the
 * corners. A particle at (0.05, 0.02) is genuinely a y-min surface particle near
 * the inlet corner; a priority scan that tests x first would label it XMin and
 * hand it n_out = (-1,0,0), up to 90 degrees away from the direction along which
 * GhostCreationInESPH actually placed its ghost (that direction comes from the
 * kernel-gradient deficiency). The face priority below is only a tie-break for
 * exactly equal distances, which keeps the assignment deterministic.
 *
 */
inline OpenFace classifyOpenFace(const Vecd &pos, const CompressibleCylinderConfig &cfg)
{
    const Real band = static_cast<Real>(cfg.boundary_n_layers) * cfg.dp;
    const Real distance[4] = {
        pos[0],          // x-min
        cfg.DL - pos[0], // x-max
        pos[1],          // y-min
        cfg.DH - pos[1], // y-max
    };

    OpenFace nearest = OpenFace::None;
    Real nearest_distance = std::numeric_limits<Real>::max();
    for (size_t face = 0; face != 4; ++face)
    {
        // Strict '<' keeps the declaration order (x-min > x-max > y-min > y-max)
        // as the tie-break for exactly equal distances.
        const OpenFace candidate = static_cast<OpenFace>(face);
        if (isActiveOpenFace(candidate, cfg) && distance[face] <= band &&
            distance[face] < nearest_distance)
        {
            nearest_distance = distance[face];
            nearest = candidate;
        }
    }
    return nearest;
}

/**
 * Keep Indicator = 1 only on active non-periodic open faces.
 *
 * FreeSurfaceIndicationComplex must run first: it uses the cylinder solid
 * contact to restore the position divergence around the cylinder, which is what
 * makes the cylinder neighbourhood separable from a genuine open face here.
 */
class OpenBoundaryIndicatorMask : public LocalDynamics
{
  public:
    OpenBoundaryIndicatorMask(SPHBody &sph_body, const CompressibleCylinderConfig &cfg)
        : LocalDynamics(sph_body),
          cfg_(cfg),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          indicator_(particles_->getVariableDataByName<int>("Indicator")) {}

    void update(size_t index_i, Real dt = 0.0)
    {
        if (indicator_[index_i] != 1)
        {
            return;
        }
        // Only the actual cylinder void is a solid-contact boundary. The z
        // stencil is closed by the periodic cell-linked-list topology, so this
        // mask must never manufacture a z far-field owner.
        if (insideCylinderCore(pos_[index_i], cfg_) ||
            classifyOpenFace(pos_[index_i], cfg_) == OpenFace::None)
        {
            indicator_[index_i] = 0;
        }
    }

  private:
    const CompressibleCylinderConfig &cfg_;
    Vecd *pos_;
    int *indicator_;
};

struct OpenBoundaryMaskResult
{
    std::array<size_t, 4> owners_per_face{{0, 0, 0, 0}};
    size_t leaked_cylinder = 0;
    size_t leaked_interior = 0;
    size_t kept_indicator = 0;
    bool pass = false;
};

/** Report the mask outcome; call after the mask has been executed. */
inline OpenBoundaryMaskResult reportOpenBoundaryMask(BaseParticles &particles,
                                                     const CompressibleCylinderConfig &cfg,
                                                     const std::string &stage)
{
    OpenBoundaryMaskResult result;
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    int *indicator = particles.getVariableDataByName<int>("Indicator");

    for (size_t i = 0; i != particles.TotalRealParticles(); ++i)
    {
        if (indicator[i] != 1)
        {
            continue;
        }
        ++result.kept_indicator;
        const OpenFace face = classifyOpenFace(pos[i], cfg);
        if (insideCylinderCore(pos[i], cfg))
        {
            ++result.leaked_cylinder;
        }
        else if (face == OpenFace::None)
        {
            ++result.leaked_interior;
        }
        else
        {
            ++result.owners_per_face[static_cast<size_t>(face)];
        }
    }

    const size_t face_total = std::accumulate(result.owners_per_face.begin(),
                                              result.owners_per_face.end(), size_t(0));
    bool active_face_owners_ok = true;
    for (size_t face = 0; face != result.owners_per_face.size(); ++face)
    {
        const bool active = isActiveOpenFace(static_cast<OpenFace>(face), cfg);
        active_face_owners_ok = active_face_owners_ok &&
                                (active ? result.owners_per_face[face] > 0
                                        : result.owners_per_face[face] == 0);
    }
    result.pass = result.leaked_cylinder == 0 && result.leaked_interior == 0 &&
                  active_face_owners_ok && result.kept_indicator == face_total;

    std::cout << "[Cylinder3DCompressible][OpenMaskGate] " << stage
              << " kept=" << result.kept_indicator
              << " x_min=" << result.owners_per_face[0]
              << " x_max=" << result.owners_per_face[1]
              << " y_min=" << result.owners_per_face[2]
              << " y_max=" << result.owners_per_face[3]
              << " leaked_cylinder=" << result.leaked_cylinder
              << " leaked_interior=" << result.leaked_interior
              << " pass=" << (result.pass ? "yes" : "NO") << std::endl;
    return result;
}

/** Per-face ghost counts plus real -> ghost mapping uniqueness. */
struct GhostMapReport
{
    std::array<size_t, 4> ghosts_per_face{{0, 0, 0, 0}};
    size_t unclassified = 0;
    size_t duplicate_real_owners = 0;
    size_t duplicate_ghost_indices = 0;
    size_t inadmissible_ghosts = 0;
    bool pass = false;
};

struct FaceDirectedGhostReport
{
    size_t corrected = 0;
    size_t missing_link = 0;
    size_t invalid_projection = 0;
    bool pass = false;
};

/**
 * In y-wall mode, an x/y corner has an x far-field state and a y solid-contact
 * stencil. GhostCreationInESPH is intentionally geometry-agnostic: it derives
 * both its ghost direction and distance from the *fluid* shape, so its generic
 * path can incorrectly place that x-state ghost by the closer y face. This
 * correction retains the shared ghost allocation and configuration, but makes
 * the one remaining open-face contract explicit: x owner, x distance, x link
 * direction, and x far-field state all agree. It is never used in the default
 * all-far-field mode.
 */
class XFaceDirectedGhostCorrection : public LocalDynamics, public DataDelegateInner
{
  public:
    XFaceDirectedGhostCorrection(BaseInnerRelation &inner_relation,
                                 GhostCreationInESPH &ghost_creation,
                                 const CompressibleCylinderConfig &cfg)
        : LocalDynamics(inner_relation.getSPHBody()),
          DataDelegateInner(inner_relation),
          cfg_(cfg),
          ghost_data_(ghost_creation.real_and_ghost_particle_data_),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure"))
    {
        if (cfg_.y_boundary_mode != YBoundaryMode::Wall)
        {
            throw std::runtime_error("XFaceDirectedGhostCorrection is valid only for y_boundary=wall.");
        }
    }

    FaceDirectedGhostReport exec(bool emit_report = false)
    {
        FaceDirectedGhostReport report;
        const size_t total_real_particles = particles_->TotalRealParticles();
        for (auto &entry : ghost_data_)
        {
            const OpenFace face = classifyOpenFace(pos_[entry.real_index_], cfg_);
            if (face != OpenFace::XMin && face != OpenFace::XMax)
            {
                ++report.invalid_projection;
                continue;
            }
            const Vecd e_face = -openFaceOutwardNormal(face); // ghost -> owner convention.
            const Real distance_to_face =
                face == OpenFace::XMin ? pos_[entry.real_index_][0]
                                       : cfg_.DL - pos_[entry.real_index_][0];
            if (distance_to_face <= TinyReal || !isFiniteRealValue(distance_to_face))
            {
                ++report.invalid_projection;
                continue;
            }

            // Keep only the kernel deficit projected onto the single active
            // face. The y component is closed by physical y-wall contact, not
            // by an x far-field ghost.
            Vecd real_neighbor_gradient = Vecd::Zero();
            Neighborhood &neighborhood = inner_configuration_[entry.real_index_];
            size_t ghost_slot = neighborhood.current_size_;
            for (size_t n = 0; n != neighborhood.current_size_; ++n)
            {
                const size_t index_j = neighborhood.j_[n];
                if (index_j < total_real_particles)
                {
                    real_neighbor_gradient +=
                        neighborhood.dW_ij_[n] * Vol_[index_j] * neighborhood.e_ij_[n];
                }
                else if (index_j == entry.ghost_index_)
                {
                    ghost_slot = n;
                }
            }
            const Real face_projection = real_neighbor_gradient.dot(e_face);
            if (ghost_slot == neighborhood.current_size_ || face_projection <= TinyReal ||
                !isFiniteRealValue(face_projection))
            {
                ++report.missing_link;
                continue;
            }

            // Mirror-link distance is twice the owner-to-face distance, matching
            // GhostCreationInESPH's established HLLC relation convention.
            neighborhood.r_ij_[ghost_slot] = 2.0 * distance_to_face;
            neighborhood.e_ij_[ghost_slot] = e_face;
            neighborhood.dW_ij_[ghost_slot] = -face_projection / (Vol_[entry.ghost_index_] + TinyReal);
            pos_[entry.ghost_index_] = pos_[entry.real_index_] - distance_to_face * e_face;
            entry.e_ij_ghost_ = e_face;
            ++report.corrected;
        }
        report.pass = report.corrected == ghost_data_.size() && report.missing_link == 0 &&
                      report.invalid_projection == 0;
        if (emit_report)
        {
            std::cout << "[Cylinder3DCompressible][FaceDirectedGhostGate] corrected=" << report.corrected
                      << " missing_link=" << report.missing_link
                      << " invalid_projection=" << report.invalid_projection
                      << " pass=" << (report.pass ? "yes" : "NO") << std::endl;
        }
        return report;
    }

  private:
    const CompressibleCylinderConfig &cfg_;
    std::vector<RealAndGhostParticleData> &ghost_data_;
    Vecd *pos_;
    Real *Vol_;
};

/**
 * Mach-aware far-field ghost state for active non-periodic open faces.
 *
 * The ghost primitive state blends the inner-neighbourhood average with the
 * freestream using the shared kernel weight summation, as the existing 2D fully
 * compressible example does, and is written through the single seven-field entry
 * point so Mass / Momentum / TotalEnergy stay consistent.
 */
class Cylinder3DCompressibleBoundaryCondition : public GhostBoundaryConditionSetupInESPH
{
  public:
    Cylinder3DCompressibleBoundaryCondition(BaseInnerRelation &inner_relation,
                                            GhostCreationInESPH &ghost_creation,
                                            const CompressibleCylinderConfig &cfg)
        : GhostBoundaryConditionSetupInESPH(inner_relation, ghost_creation),
          cfg_(cfg),
          freestream_(makeFreestreamState(cfg)),
          state_(makeCompressibleStateView(*particles_, cfg.gamma))
    {
        // The base constructor already dispatched setupBoundaryTypes() through the
        // base vtable, i.e. as a no-op. Classify here, where cfg_ is initialised.
        assignBoundaryTypes();
    }
    virtual ~Cylinder3DCompressibleBoundaryCondition() {}

    /**
     * Every ghost owner must be an active open face, so all of them get
     * BoundaryType 9. An owner that the mask should have removed is a hard
     * error: continuing would silently create a far-field ghost on the cylinder
     * rather than on the cylinder solid-contact boundary.
     */
    void assignBoundaryTypes()
    {
        for (const auto &entry : real_and_ghost_particle_data_)
        {
            const size_t index_i = entry.real_index_;
            const OpenFace face = classifyOpenFace(pos_[index_i], cfg_);
            if (face == OpenFace::None || !isActiveOpenFace(face, cfg_) ||
                insideCylinderCore(pos_[index_i], cfg_))
            {
                throw std::runtime_error(
                    "Cylinder3DCompressibleBoundaryCondition: ghost owner is not on an active open face; "
                    "OpenBoundaryIndicatorMask must run before GhostCreationInESPH is constructed.");
            }
            boundary_type_[index_i] = 9;
        }
    }

    void applyFarFieldBoundary(size_t ghost_index, size_t index_i) override
    {
        const OpenFace face = classifyOpenFace(pos_[index_i], cfg_);
        const Vecd n_out = openFaceOutwardNormal(face);
        const Real velocity_farfield_normal = freestream_.vel.dot(n_out);
        const Real velocity_boundary_normal = state_.vel[index_i].dot(n_out);

        // Inner-neighbourhood averages over real particles only: a ghost must not
        // be averaged into the state that defines it.
        Real inner_weight_summation = W0_ * state_.Vol[index_i];
        Real rho_summation = 0.0;
        Real p_summation = 0.0;
        Real vel_normal_summation = 0.0;
        Vecd vel_tangential_summation = Vecd::Zero();
        size_t neighbor_count = 0;
        Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            const size_t index_j = inner_neighborhood.j_[n];
            if (index_j >= particles_->TotalRealParticles())
            {
                continue;
            }
            inner_weight_summation += inner_neighborhood.W_ij_[n] * state_.Vol[index_j];
            rho_summation += state_.rho[index_j];
            p_summation += state_.p[index_j];
            vel_normal_summation += state_.vel[index_j].dot(n_out);
            vel_tangential_summation += state_.vel[index_j] - state_.vel[index_j].dot(n_out) * n_out;
            ++neighbor_count;
        }

        if (neighbor_count == 0)
        {
            // No real neighbour to average: use the pure freestream rather than
            // dividing by an empty stencil.
            setFreestreamState(state_, ghost_index, freestream_);
            return;
        }

        const Real count = static_cast<Real>(neighbor_count);
        const Real rho_average = rho_summation / count;
        const Real p_average = p_summation / count;
        const Real vel_normal_average = vel_normal_summation / count;
        const Vecd vel_tangential_average = vel_tangential_summation / count;

        const Real w = SMIN(Real(1.0), SMAX(Real(0.0), inner_weight_summation));
        const Real one_minus_w = 1.0 - w;
        const bool inflow = velocity_boundary_normal <= 0.0;

        Real rho_ghost = 0.0;
        Real p_ghost = 0.0;
        Vecd vel_ghost = Vecd::Zero();

        const Real local_sound_speed = localSoundSpeed(state_, index_i);
        if (std::fabs(velocity_boundary_normal) >= local_sound_speed)
        {
            // Supersonic characteristic contract: prescribe the complete state
            // at inflow and extrapolate the complete owner state at outflow.
            if (inflow)
            {
                rho_ghost = freestream_.rho;
                p_ghost = freestream_.p;
                vel_ghost = freestream_.vel;
            }
            else
            {
                rho_ghost = state_.rho[index_i];
                p_ghost = state_.p[index_i];
                vel_ghost = state_.vel[index_i];
            }
        }
        else
        {
            rho_ghost = rho_average * w + freestream_.rho * one_minus_w;
            p_ghost = p_average * w + freestream_.p * one_minus_w;
            const Real vel_normal = vel_normal_average * w + velocity_farfield_normal * one_minus_w;
            // Inflow: the tangential part is prescribed by the freestream.
            // Outflow: the tangential part is convected out.
            // Both branches are materialised into Vecd -- an Eigen expression
            // ternary would need both arms to share one expression type.
            const Vecd vel_tangential =
                inflow ? Vecd(freestream_.vel - velocity_farfield_normal * n_out)
                       : vel_tangential_average;
            vel_ghost = vel_normal * n_out + vel_tangential;
        }

        if (!isFiniteRealValue(rho_ghost) || !isFiniteRealValue(p_ghost) || rho_ghost <= 0.0 ||
            p_ghost <= 0.0)
        {
            std::cout << "[Cylinder3DCompressible][FarFieldGhostStateError]"
                      << " face=" << static_cast<size_t>(face)
                      << " owner=" << index_i
                      << " ghost=" << ghost_index
                      << " branch=" << (std::fabs(velocity_boundary_normal) >= local_sound_speed ? "supersonic" : "subsonic")
                      << " rho_owner=" << state_.rho[index_i]
                      << " p_owner=" << state_.p[index_i]
                      << " rho_average=" << rho_average
                      << " p_average=" << p_average
                      << " rho_ghost=" << rho_ghost
                      << " p_ghost=" << p_ghost
                      << " normal_velocity=" << velocity_boundary_normal
                      << " sound_speed=" << local_sound_speed
                      << " neighbors=" << neighbor_count << std::endl;
        }
        setCompressiblePrimitiveState(state_, ghost_index, rho_ghost, vel_ghost, p_ghost);
    }

    GhostMapReport reportGhostMap(const std::string &stage)
    {
        GhostMapReport report;
        std::vector<size_t> real_owners;
        std::vector<size_t> ghost_indices;
        real_owners.reserve(real_and_ghost_particle_data_.size());
        ghost_indices.reserve(real_and_ghost_particle_data_.size());

        for (const auto &entry : real_and_ghost_particle_data_)
        {
            const OpenFace face = classifyOpenFace(pos_[entry.real_index_], cfg_);
            if (face == OpenFace::None)
            {
                ++report.unclassified;
            }
            else
            {
                ++report.ghosts_per_face[static_cast<size_t>(face)];
            }
            if (!isThermodynamicallyAdmissible(state_, entry.ghost_index_))
            {
                if (report.inadmissible_ghosts == 0)
                {
                    std::cout << "[Cylinder3DCompressible][GhostState] first inadmissible"
                              << " face=" << openFaceName(face)
                              << " real=" << entry.real_index_
                              << " ghost=" << entry.ghost_index_
                              << " rho_real=" << state_.rho[entry.real_index_]
                              << " p_real=" << state_.p[entry.real_index_]
                              << " rho_ghost=" << state_.rho[entry.ghost_index_]
                              << " p_ghost=" << state_.p[entry.ghost_index_]
                              << " E_ghost=" << state_.E[entry.ghost_index_] << std::endl;
                }
                ++report.inadmissible_ghosts;
            }
            real_owners.push_back(entry.real_index_);
            ghost_indices.push_back(entry.ghost_index_);
        }

        std::sort(real_owners.begin(), real_owners.end());
        std::sort(ghost_indices.begin(), ghost_indices.end());
        report.duplicate_real_owners =
            real_owners.size() -
            static_cast<size_t>(std::distance(real_owners.begin(),
                                              std::unique(real_owners.begin(), real_owners.end())));
        report.duplicate_ghost_indices =
            ghost_indices.size() -
            static_cast<size_t>(std::distance(ghost_indices.begin(),
                                              std::unique(ghost_indices.begin(), ghost_indices.end())));

        bool active_face_ghosts_ok = true;
        for (size_t face = 0; face != report.ghosts_per_face.size(); ++face)
        {
            const bool active = isActiveOpenFace(static_cast<OpenFace>(face), cfg_);
            active_face_ghosts_ok = active_face_ghosts_ok &&
                                    (active ? report.ghosts_per_face[face] > 0
                                            : report.ghosts_per_face[face] == 0);
        }
        report.pass = report.unclassified == 0 && report.duplicate_real_owners == 0 &&
                      report.duplicate_ghost_indices == 0 && report.inadmissible_ghosts == 0 &&
                      active_face_ghosts_ok;

        std::cout << "[Cylinder3DCompressible][GhostMapGate] " << stage
                  << " total=" << real_and_ghost_particle_data_.size()
                  << " x_min=" << report.ghosts_per_face[0]
                  << " x_max=" << report.ghosts_per_face[1]
                  << " y_min=" << report.ghosts_per_face[2]
                  << " y_max=" << report.ghosts_per_face[3]
                  << " unclassified=" << report.unclassified
                  << " dup_real=" << report.duplicate_real_owners
                  << " dup_ghost=" << report.duplicate_ghost_indices
                  << " inadmissible=" << report.inadmissible_ghosts
                  << " pass=" << (report.pass ? "yes" : "NO") << std::endl;
        return report;
    }

  private:
    const CompressibleCylinderConfig &cfg_;
    CompressibleFreestreamState freestream_;
    CompressibleConservativeStateView state_;
};

} // namespace cylinder_3d_compressible
} // namespace SPH

#endif // CYLINDER_3D_COMPRESSIBLE_BOUNDARY_HPP
