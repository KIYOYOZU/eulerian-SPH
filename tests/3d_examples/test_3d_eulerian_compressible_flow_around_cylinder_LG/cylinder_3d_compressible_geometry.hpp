#ifndef CYLINDER_3D_COMPRESSIBLE_GEOMETRY_HPP
#define CYLINDER_3D_COMPRESSIBLE_GEOMETRY_HPP

/**
 * Geometry and cylinder solid for the fully
 * compressible cylinder case.
 *
 * The fluid shape is exactly the physical box [0,DL] x [0,DH] x [0,DW] minus
 * the cylinder: no sponge / padding layer of real particles outside the physical
 * x/y faces, because such layers would combine an unvalidated relaxation
 * boundary with the ghost flux boundary.
 *
 * n_wall points from the solid to the fluid and has a zero z component; the
 * fluid-owner -> wall-neighbour direction e_pair must satisfy
 * dot(e_pair, n_wall) > 0. The flux normal is n_flux = -e_pair.
 */

#include "cylinder_3d_compressible_state.hpp"
#include "sphinxsys.h"

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace SPH
{
namespace cylinder_3d_compressible
{

/**
 * Half-height of the cylinder mesh used for both the fluid void and the solid
 * wall body, extended a few dp beyond the physical z-span.
 *
 * The z-min/z-max planes are far fields, while the cylinder continues beyond
 * the fluid span so its radial void remains well-defined at those planes. The
 * mesh height is finite only because the shape API needs one. Using exactly
 * 0.5*DW puts the two flat end caps exactly on z=0 and z=DW. TriangleMeshShapeCylinder's
 * containment test is a mesh ray cast, which is markedly less robust near a
 * flat cap than on the cylindrical side, and at dp=0.02 a lattice sample point
 * lands EXACTLY on that plane (DW/dp = 12.5, an exact half-integer -- see
 * BW = 4*dp in the lattice bounding box). The two effects together silently
 * misclassified 76 cylinder-interior fluid particles as inside the fluid
 * (verified: pos=(0.41,0.47,0.25), radial distance 0.0949 < radius 0.10, at
 * z = DW exactly), which then produced fluid/wall contact pairs pointing INTO
 * the solid (dot(e_pair,n_wall) = -1, not noise).
 *
 * Extending the half-height moves both caps outside [0, DW] without changing
 * the radial containment test anywhere inside that range -- the cylinder's
 * cross-section is z-independent by construction -- so no cap plane can ever
 * coincide with a sampled z, at this or any other resolution.
 */
inline Real cylinderMeshHalfHeight(const CompressibleCylinderConfig &cfg)
{
    return 0.5 * cfg.DW + 4.0 * cfg.dp;
}

class Cylinder3DFluidBlock : public ComplexShape
{
  public:
    explicit Cylinder3DFluidBlock(const std::string &shape_name, const CompressibleCylinderConfig &cfg)
        : ComplexShape(shape_name)
    {
        // Physical box only -- x/y faces use far-field ghosts and z is periodic.
        const Vecd halfsize(0.5 * cfg.DL, 0.5 * cfg.DH, 0.5 * cfg.DW);
        const Vecd translation(0.5 * cfg.DL, 0.5 * cfg.DH, 0.5 * cfg.DW);
        add<GeometricShapeBox>(Transform(translation), halfsize, "OuterBoundary");

        const Vecd cylinder_center(cfg.cylinder_center_x, cfg.cylinder_center_y, 0.5 * cfg.DW);
        subtract<TriangleMeshShapeCylinder>(Vecd(0.0, 0.0, 1.0), cfg.cylinder_radius,
                                            cylinderMeshHalfHeight(cfg), 48,
                                            cylinder_center, "CylinderVoid");
    }
};

/**
 * Solid cylinder wall for a spanwise-periodic domain.
 *
 * Its z-particle placement is supplied by PeriodicCylinderWallLattice below.
 * That generator gives the finite physical span an exact periodic lattice even
 * when DW/dp is not an integer.
 */
class Cylinder3DWallBlock : public ComplexShape
{
  public:
    explicit Cylinder3DWallBlock(const std::string &shape_name, const CompressibleCylinderConfig &cfg)
        : ComplexShape(shape_name)
    {
        const Vecd cylinder_center(cfg.cylinder_center_x, cfg.cylinder_center_y, 0.5 * cfg.DW);
        add<TriangleMeshShapeCylinder>(Vecd(0.0, 0.0, 1.0), cfg.cylinder_radius,
                                       0.5 * cfg.DW, 48,
                                       cylinder_center, "CylinderWall");
    }
};

/** Tags for fluid/wall generators whose z centres close exactly under one DW shift. */
class PeriodicCylinderFluidLattice;
class PeriodicCylinderWallLattice;
class YWallVolumeLattice;

/**
 * Radial cylinder normal with an exactly zero z component.
 *
 * The body-shape normal would pick up finite-mesh end caps. The periodic wall
 * particle generator excludes those caps, and its remaining solid contact is
 * radial across the z seam.
 */
class Cylinder3DRadialWallNormal : public NormalDirectionFromBodyShape
{
  public:
    Cylinder3DRadialWallNormal(SPHBody &sph_body, const CompressibleCylinderConfig &cfg)
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
    const CompressibleCylinderConfig &cfg_;
};

/** Constant solid-to-fluid normals for the two filled y-wall slabs. */
class YWallNormal : public LocalDynamics
{
  public:
    explicit YWallNormal(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          n_(particles_->registerStateVariableData<Vecd>("NormalDirection")) {}

    void update(size_t index_i, Real dt = 0.0)
    {
        n_[index_i] = pos_[index_i][1] < 0.0 ? Vecd(0.0, 1.0, 0.0) : Vecd(0.0, -1.0, 0.0);
    }

  private:
    Vecd *pos_, *n_;
};

/**
 * Analytic fluid-domain normal, replacing NormalDirectionFromBodyShape.
 *
 * WHY NOT the shared version: it calls Shape::findNormalDirection() plus
 * findSignedDistance() per particle (general_geometric.cpp:20-25), i.e. two full
 * boolean-tree queries against ComplexShape = box MINUS TriangleMeshShapeCylinder.
 * The cylinder leg resolves through Simbody's mesh closest-point search, and a
 * component level set cannot help: defineComponentLevelSetShape() does not
 * replace initial_shape_ (base_body.hpp:23-27), so the query still walks the
 * boolean tree. Measured on BSCC-M9 at dp=0.02 (62.5k particles):
 * 37 minutes without finishing, ~7 cores busy, before a single time step.
 *
 * The fluid domain boundary relevant to the ghost path is analytically known --
 * four x/y planar faces and one cylindrical surface -- so the nearest one can be
 * selected directly. The z direction is periodic and must not enter this normal.
 * This is the
 * same radial-normal substitution Cylinder3DRadialWallNormal already makes for the wall, and
 * it matches the upstream 2D ESPH prototype (test_2d_eulerian_supersonic_flow_new_BC),
 * which never runs a normal-direction pass over the fluid at all and instead
 * evaluates the normal only for the few ghost owners that need it.
 *
 * Consumers: the fluid NormalDirection is written to VTP for visualisation only.
 * Every gate reads the WALL normal (analytic radial), and the open faces use
 * openFaceOutwardNormal(classifyOpenFace(...)). Sign convention follows the shared
 * class: the normal points out of the fluid, so it is -radial on the cylinder
 * surface (into the void) and +axis on the outer box faces.
 */
class Cylinder3DFluidNormal : public NormalDirectionFromBodyShape
{
  public:
    Cylinder3DFluidNormal(SPHBody &sph_body, const CompressibleCylinderConfig &cfg)
        : NormalDirectionFromBodyShape(sph_body), cfg_(cfg) {}

    void update(size_t index_i, Real dt = 0.0)
    {
        const Vecd &pos = pos_[index_i];
        // Distance to each non-periodic box face, with its outward normal.
        const Real face_distance[4] = {
            pos[0], cfg_.DL - pos[0],
            pos[1], cfg_.DH - pos[1]};
        static const Vecd face_normal[4] = {
            Vecd(-1.0, 0.0, 0.0), Vecd(1.0, 0.0, 0.0),
            Vecd(0.0, -1.0, 0.0), Vecd(0.0, 1.0, 0.0)};

        Real nearest_distance = face_distance[0];
        Vecd nearest_normal = face_normal[0];
        for (size_t face = 1; face != 4; ++face)
        {
            if (face_distance[face] < nearest_distance)
            {
                nearest_distance = face_distance[face];
                nearest_normal = face_normal[face];
            }
        }

        // The cylinder void, competing on the same nearest-surface basis. Its
        // outward-from-fluid normal points toward the cylinder axis.
        const Vecd radial(pos[0] - cfg_.cylinder_center_x, pos[1] - cfg_.cylinder_center_y, 0.0);
        const Real radial_distance = radial.norm();
        const Real cylinder_distance = radial_distance - cfg_.cylinder_radius;
        if (cylinder_distance < nearest_distance)
        {
            nearest_distance = cylinder_distance;
            nearest_normal = radial_distance > TinyReal ? Vecd(-radial / radial_distance)
                                                        : Vecd(-1.0, 0.0, 0.0);
        }

        n_[index_i] = nearest_normal;
        n0_[index_i] = nearest_normal;
        // Signed distance to the fluid boundary: positive inside the fluid.
        phi_[index_i] = nearest_distance;
        phi0_[index_i] = nearest_distance;
    }

  private:
    const CompressibleCylinderConfig &cfg_;
};

/**
 * Freestream initial condition written through the single seven-field entry point.
 *
 * The perturbation SCALES the streamwise component per particle -- the direction
 * stays exactly +x, so the y-symmetry is broken by the magnitude scatter, not by
 * a tilt. Density and pressure are written at their freestream values, and
 * setCompressiblePrimitiveState() derives TotalEnergy from the perturbed speed,
 * so the recovered internal energy is exactly p_inf/(gamma-1) for every particle
 * and the state is thermodynamically admissible at any amplitude in [0, 1).
 */
class Cylinder3DCompressibleInitialCondition : public fluid_dynamics::FluidInitialCondition
{
  public:
    Cylinder3DCompressibleInitialCondition(SPHBody &sph_body, const CompressibleCylinderConfig &cfg)
        : fluid_dynamics::FluidInitialCondition(sph_body),
          cfg_(cfg),
          freestream_(makeFreestreamState(cfg)),
          state_(makeCompressibleStateView(*particles_, cfg.gamma)),
          noise_amplitude_(cfg.velocity_noise_amplitude) {}

    void update(size_t index_i, Real dt = 0.0)
    {
        Vecd velocity = freestream_.vel;
        velocity[0] *= 1.0 + noise_amplitude_ * signedIndexNoise(index_i);
        setCompressiblePrimitiveState(state_, index_i, freestream_.rho, velocity, freestream_.p);
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

    const CompressibleCylinderConfig &cfg_;
    CompressibleFreestreamState freestream_;
    CompressibleConservativeStateView state_;
    Real noise_amplitude_;
};

inline bool insideCylinderCore(const Vecd &pos, const CompressibleCylinderConfig &cfg,
                               Real radius_offset = 0.0)
{
    const Real dx = pos[0] - cfg.cylinder_center_x;
    const Real dy = pos[1] - cfg.cylinder_center_y;
    const Real r = std::sqrt(dx * dx + dy * dy);
    return r <= cfg.cylinder_radius + radius_offset;
}

inline bool outsidePhysicalBox(const Vecd &pos, const CompressibleCylinderConfig &cfg, Real slack)
{
    return pos[0] < -slack || pos[0] > cfg.DL + slack ||
           pos[1] < -slack || pos[1] > cfg.DH + slack ||
           pos[2] < -slack || pos[2] > cfg.DW + slack;
}

inline bool reportFluidFiniteState(BaseParticles &particles, const std::string &stage)
{
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    Real *pressure = particles.getVariableDataByName<Real>("Pressure");
    Real *energy = particles.getVariableDataByName<Real>("TotalEnergy");
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
    size_t speed_max_index = total;

    for (size_t i = 0; i != total; ++i)
    {
        const bool finite = isFiniteVecValue(pos[i]) && isFiniteRealValue(vol[i]) &&
                            isFiniteRealValue(rho[i]) && isFiniteRealValue(mass[i]) &&
                            isFiniteRealValue(pressure[i]) && isFiniteRealValue(energy[i]) &&
                            isFiniteVecValue(velocity[i]) && isFiniteVecValue(momentum[i]);
        if (!finite && first_bad == total)
        {
            first_bad = i;
        }
        if (isFiniteRealValue(rho[i]))
        {
            if (rho[i] < rho_min)
            {
                rho_min = rho[i];
                rho_min_index = i;
            }
            rho_max = SMAX(rho_max, rho[i]);
        }
        if (isFiniteRealValue(pressure[i]))
        {
            p_min = SMIN(p_min, pressure[i]);
            p_max = SMAX(p_max, pressure[i]);
        }
        if (isFiniteVecValue(velocity[i]))
        {
            const Real speed = velocity[i].norm();
            if (speed > speed_max)
            {
                speed_max = speed;
                speed_max_index = i;
            }
        }
    }

    std::cout << "[Cylinder3DCompressible][State] " << stage
              << " total=" << total
              << " rho=[" << rho_min << ", " << rho_max << "]"
              << " p=[" << p_min << ", " << p_max << "]"
              << " max|u|=" << speed_max
              << " finite=" << (first_bad == total ? "yes" : "NO") << std::endl;
    if (first_bad != total)
    {
        std::cout << "[Cylinder3DCompressible][BadParticle] i=" << first_bad
                  << " pos=(" << pos[first_bad][0] << ", " << pos[first_bad][1]
                  << ", " << pos[first_bad][2] << ")"
                  << " rho=" << rho[first_bad]
                  << " p=" << pressure[first_bad]
                  << " E=" << energy[first_bad]
                  << " vel=(" << velocity[first_bad][0] << ", " << velocity[first_bad][1]
                  << ", " << velocity[first_bad][2] << ")" << std::endl;
    }
    if (speed_max_index != total && rho_min_index != total)
    {
        std::cout << "[Cylinder3DCompressible][Extrema]"
                  << " max_speed_i=" << speed_max_index
                  << " max_speed_pos=(" << pos[speed_max_index][0] << ", " << pos[speed_max_index][1]
                  << ", " << pos[speed_max_index][2] << ")"
                  << " rho_min_i=" << rho_min_index
                  << " rho_min_pos=(" << pos[rho_min_index][0] << ", " << pos[rho_min_index][1]
                  << ", " << pos[rho_min_index][2] << ")" << std::endl;
    }
    return first_bad == total;
}

struct GeometryGateResult
{
    size_t fluid_particles = 0;
    size_t wall_particles = 0;
    size_t fluid_outside_box = 0;
    size_t cylinder_internal_fluid = 0;
    size_t x_min_surface = 0;
    size_t x_max_surface = 0;
    size_t y_surface = 0;
    size_t z_surface = 0;
    size_t cylinder_contact_particles = 0;
    size_t pair_direction_failures = 0;
    Real min_fluid_wall_distance = std::numeric_limits<Real>::infinity();
    Real z_end_contact_mean = 0.0;
    Real z_mid_contact_mean = 0.0;
    Real wall_z_min = std::numeric_limits<Real>::infinity();
    Real wall_z_max = -std::numeric_limits<Real>::infinity();
    Real wall_max_abs_normal_z = 0.0;
    bool pass = false;
};

/**
 * Geometry gate: no solid particles outside the physical box, no fluid inside
 * the cylinder core, spanwise-consistent contact counts, and a valid
 * dot(e_pair, n_wall) > 0 for every effective contact pair.
 */
inline GeometryGateResult reportGeometryGate(FluidBody &fluid_body, SolidBody &wall_body,
                                             BaseInnerRelation &fluid_inner,
                                             ContactRelation &fluid_wall_contact,
                                             const CompressibleCylinderConfig &cfg)
{
    GeometryGateResult gate;
    BaseParticles &fluid_particles = fluid_body.getBaseParticles();
    BaseParticles &wall_particles = wall_body.getBaseParticles();
    Vecd *pos = fluid_particles.getVariableDataByName<Vecd>("Position");
    Vecd *wall_pos = wall_particles.getVariableDataByName<Vecd>("Position");
    Vecd *wall_normal = wall_particles.getVariableDataByName<Vecd>("NormalDirection");
    gate.fluid_particles = fluid_particles.TotalRealParticles();
    gate.wall_particles = wall_particles.TotalRealParticles();

    size_t z_end_count = 0;
    size_t z_mid_count = 0;
    Real z_end_sum = 0.0;
    Real z_mid_sum = 0.0;
    const Real z_band = static_cast<Real>(cfg.boundary_n_layers) * cfg.dp;
    const Real mid_low = 0.5 * cfg.DW - z_band;
    const Real mid_high = 0.5 * cfg.DW + z_band;
    // One dp of slack absorbs the lattice half-spacing offset at the faces.
    const Real box_slack = cfg.dp;

    for (size_t i = 0; i != gate.wall_particles; ++i)
    {
        gate.wall_z_min = SMIN(gate.wall_z_min, wall_pos[i][2]);
        gate.wall_z_max = SMAX(gate.wall_z_max, wall_pos[i][2]);
        gate.wall_max_abs_normal_z = SMAX(gate.wall_max_abs_normal_z, std::fabs(wall_normal[i][2]));
    }

    for (size_t i = 0; i != gate.fluid_particles; ++i)
    {
        if (outsidePhysicalBox(pos[i], cfg, box_slack))
        {
            ++gate.fluid_outside_box;
        }
        if (insideCylinderCore(pos[i], cfg, -0.25 * cfg.dp))
        {
            ++gate.cylinder_internal_fluid;
        }

        // Face candidates are counted by POSITION, not by normal direction. The
        // cylinder-adjacent particles carry radial normals whose x/y mix would
        // otherwise land them in the x-min/x-max/y buckets depending on angle,
        // so a normal-based count can be satisfied without any particle actually
        // sitting on a physical outer face.
        const Real face_band = cfg.dp;
        const bool near_cylinder = insideCylinderCore(pos[i], cfg, face_band);
        if (!near_cylinder)
        {
            if (pos[i][0] <= face_band)
            {
                ++gate.x_min_surface;
            }
            if (pos[i][0] >= cfg.DL - face_band)
            {
                ++gate.x_max_surface;
            }
            if (pos[i][1] <= face_band || pos[i][1] >= cfg.DH - face_band)
            {
                ++gate.y_surface;
            }
            if (pos[i][2] <= face_band || pos[i][2] >= cfg.DW - face_band)
            {
                ++gate.z_surface;
            }
        }

        Neighborhood &contact_neighborhood = fluid_wall_contact.contact_configuration_[0][i];
        const size_t contact_neighbors = contact_neighborhood.current_size_;
        if (contact_neighbors > 0)
        {
            ++gate.cylinder_contact_particles;
        }
        for (size_t n = 0; n != contact_neighbors; ++n)
        {
            const size_t index_j = contact_neighborhood.j_[n];
            // e_pair points from the fluid owner towards the wall neighbour and
            // must agree with the outward (solid -> fluid) wall normal.
            const Vecd &e_pair = contact_neighborhood.e_ij_[n];
            // A purely tangential pair gives dot == 0 analytically and lands
            // within a few machine eps of it; only a genuinely negative
            // alignment breaks the direction contract. Same dead zone as
            // reportWallPairGate().
            if (e_pair.dot(wall_normal[index_j]) < -1.0e-12)
            {
                if (gate.pair_direction_failures == 0)
                {
                    std::cout << "[Cylinder3DCompressible][PairDirection] first failure"
                              << " fluid_i=" << i
                              << " pos=(" << pos[i][0] << ", " << pos[i][1] << ", " << pos[i][2] << ")"
                              << " wall_j=" << index_j
                              << " e_pair=(" << e_pair[0] << ", " << e_pair[1] << ", " << e_pair[2] << ")"
                              << " n_wall=(" << wall_normal[index_j][0] << ", " << wall_normal[index_j][1]
                              << ", " << wall_normal[index_j][2] << ")" << std::endl;
                }
                ++gate.pair_direction_failures;
            }
            gate.min_fluid_wall_distance = SMIN(gate.min_fluid_wall_distance, contact_neighborhood.r_ij_[n]);
        }

        if (insideCylinderCore(pos[i], cfg, cfg.dp) && !insideCylinderCore(pos[i], cfg, -cfg.dp))
        {
            if (pos[i][2] <= z_band || pos[i][2] >= cfg.DW - z_band)
            {
                z_end_sum += static_cast<Real>(contact_neighbors);
                ++z_end_count;
            }
            if (pos[i][2] >= mid_low && pos[i][2] <= mid_high)
            {
                z_mid_sum += static_cast<Real>(contact_neighbors);
                ++z_mid_count;
            }
        }
    }
    gate.z_end_contact_mean = z_end_count > 0 ? z_end_sum / static_cast<Real>(z_end_count) : 0.0;
    gate.z_mid_contact_mean = z_mid_count > 0 ? z_mid_sum / static_cast<Real>(z_mid_count) : 0.0;

    const bool z_contact_ok = gate.z_mid_contact_mean <= 0.0 ||
                              gate.z_end_contact_mean >= 0.8 * gate.z_mid_contact_mean;
    const bool wall_radial_geometry_ok = gate.wall_z_min >= -TinyReal &&
                                         gate.wall_z_max <= cfg.DW + TinyReal &&
                                         gate.wall_max_abs_normal_z <= TinyReal;
    gate.pass = gate.fluid_outside_box == 0 && gate.cylinder_internal_fluid == 0 &&
                gate.x_min_surface > 0 && gate.x_max_surface > 0 && gate.y_surface > 0 &&
                gate.cylinder_contact_particles > 0 && gate.pair_direction_failures == 0 &&
                z_contact_ok && wall_radial_geometry_ok;

    std::cout << "[Cylinder3DCompressible][GeometryGate] fluid=" << gate.fluid_particles
              << " wall=" << gate.wall_particles
              << " outside_box=" << gate.fluid_outside_box
              << " internal_fluid=" << gate.cylinder_internal_fluid
              << " x_min_surface=" << gate.x_min_surface
              << " x_max_surface=" << gate.x_max_surface
              << " y_surface=" << gate.y_surface
              << " z_surface_raw=" << gate.z_surface
              << " contact_particles=" << gate.cylinder_contact_particles
              << " pair_direction_failures=" << gate.pair_direction_failures
              << " min_fluid_wall_r=" << gate.min_fluid_wall_distance
              << " z_end_contact_mean=" << gate.z_end_contact_mean
              << " z_mid_contact_mean=" << gate.z_mid_contact_mean
              << " wall_z=[" << gate.wall_z_min << ", " << gate.wall_z_max << "]"
              << " wall_max_abs_normal_z=" << gate.wall_max_abs_normal_z
              << " pass=" << (gate.pass ? "yes" : "NO") << std::endl;
    return gate;
}

} // namespace cylinder_3d_compressible

/**
 * Generate the solid cylinder on a tensor-product lattice with a wall-only
 * spanwise spacing. The fluid retains cfg.dp. For the target DW=0.25 and
 * dp=0.02, Nz=13 and dz=0.25/13, so the first and last wall layers have the
 * same periodic separation as every interior pair. The volume measure uses
 * dx*dy*dz, consistently accounting for the anisotropic particle spacing.
 */
template <>
class ParticleGenerator<BaseParticles, cylinder_3d_compressible::PeriodicCylinderFluidLattice>
    : public ParticleGenerator<BaseParticles>
{
  public:
    ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles,
                      const cylinder_3d_compressible::CompressibleCylinderConfig &cfg)
        : ParticleGenerator<BaseParticles>(sph_body, base_particles),
          initial_shape_(sph_body.getInitialShape()), cfg_(cfg) {}

    void prepareGeometricData() override
    {
        const size_t nx = static_cast<size_t>(std::round(cfg_.DL / cfg_.dp));
        const size_t ny = static_cast<size_t>(std::round(cfg_.DH / cfg_.dp));
        const size_t nz = SMAX(size_t(1), static_cast<size_t>(std::round(cfg_.DW / cfg_.dp)));
        const Real dx = cfg_.DL / static_cast<Real>(nx);
        const Real dy = cfg_.DH / static_cast<Real>(ny);
        const Real dz = cfg_.DW / static_cast<Real>(nz);
        const Real volume = dx * dy * dz;

        for (size_t i = 0; i < nx; ++i)
        {
            const Real x = (static_cast<Real>(i) + 0.5) * dx;
            for (size_t j = 0; j < ny; ++j)
            {
                const Real y = (static_cast<Real>(j) + 0.5) * dy;
                for (size_t k = 0; k < nz; ++k)
                {
                    const Vecd position(x, y, (static_cast<Real>(k) + 0.5) * dz);
                    if (initial_shape_.checkContain(position))
                    {
                        addPositionAndVolumetricMeasure(position, volume);
                    }
                }
            }
        }
    }

  private:
    Shape &initial_shape_;
    const cylinder_3d_compressible::CompressibleCylinderConfig &cfg_;
};

/**
 * The wall uses the same periodic z lattice as the fluid but an analytic radial
 * inclusion. This deliberately bypasses the finite level-set cylinder caps:
 * every layer has the same 2-D disk, so the seam cannot acquire cap-only
 * particles or a z-dependent radial contact stencil.
 */
template <>
class ParticleGenerator<BaseParticles, cylinder_3d_compressible::PeriodicCylinderWallLattice>
    : public ParticleGenerator<BaseParticles>
{
  public:
    ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles,
                      const cylinder_3d_compressible::CompressibleCylinderConfig &cfg)
        : ParticleGenerator<BaseParticles>(sph_body, base_particles), cfg_(cfg) {}

    void prepareGeometricData() override
    {
        const size_t nx = static_cast<size_t>(std::round(cfg_.DL / cfg_.dp));
        const size_t ny = static_cast<size_t>(std::round(cfg_.DH / cfg_.dp));
        const size_t nz = SMAX(size_t(1), static_cast<size_t>(std::round(cfg_.DW / cfg_.dp)));
        const Real dx = cfg_.DL / static_cast<Real>(nx);
        const Real dy = cfg_.DH / static_cast<Real>(ny);
        const Real dz = cfg_.DW / static_cast<Real>(nz);
        const Real volume = dx * dy * dz;
        const Real radius_squared = cfg_.cylinder_radius * cfg_.cylinder_radius;

        for (size_t i = 0; i < nx; ++i)
        {
            const Real x = (static_cast<Real>(i) + 0.5) * dx;
            for (size_t j = 0; j < ny; ++j)
            {
                const Real y = (static_cast<Real>(j) + 0.5) * dy;
                const Real dx_center = x - cfg_.cylinder_center_x;
                const Real dy_center = y - cfg_.cylinder_center_y;
                if (dx_center * dx_center + dy_center * dy_center >= radius_squared)
                {
                    continue;
                }
                for (size_t k = 0; k < nz; ++k)
                {
                    addPositionAndVolumetricMeasure(
                        Vecd(x, y, (static_cast<Real>(k) + 0.5) * dz), volume);
                }
            }
        }
    }

  private:
    const cylinder_3d_compressible::CompressibleCylinderConfig &cfg_;
};

/**
 * Upper and lower physical y walls. Eulerian HLLC wall contact weighs every
 * neighbour by its VolumetricMeasure, so the wall must be a filled solid slab
 * (as the cylinder is), not a surface particle whose measure is an area. Three
 * layers reach the 2.6*dp fluid-kernel support while remaining outside the
 * physical fluid domain. In x-farfield mode the slabs continue by three layers
 * beyond both x-open faces: an x ghost closes only the fluid-fluid stencil,
 * whereas a fluid near an x/y corner still needs the missing tangential half of
 * its y-wall contact stencil. In x-periodic mode they contain physical-domain
 * cells only and their contact stencil closes through the x/z periodic cell
 * lists; retaining the extension there would wrap duplicate wall entries into
 * the domain.
 */
template <>
class ParticleGenerator<BaseParticles, cylinder_3d_compressible::YWallVolumeLattice>
    : public ParticleGenerator<BaseParticles>
{
  public:
    ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles,
                      const cylinder_3d_compressible::CompressibleCylinderConfig &cfg)
        : ParticleGenerator<BaseParticles>(sph_body, base_particles), cfg_(cfg) {}

    void prepareGeometricData() override
    {
        const size_t nx = static_cast<size_t>(std::round(cfg_.DL / cfg_.dp));
        const size_t nz = SMAX(size_t(1), static_cast<size_t>(std::round(cfg_.DW / cfg_.dp)));
        const Real dx = cfg_.DL / static_cast<Real>(nx);
        const Real dz = cfg_.DW / static_cast<Real>(nz);
        const Real volume = dx * cfg_.dp * dz;
        constexpr size_t wall_layers = 3;
        constexpr int x_extension_layers = 3;
        const int x_begin = cfg_.x_boundary_mode == cylinder_3d_compressible::XBoundaryMode::Periodic
                                ? 0
                                : -x_extension_layers;
        const int x_end = cfg_.x_boundary_mode == cylinder_3d_compressible::XBoundaryMode::Periodic
                              ? static_cast<int>(nx)
                              : static_cast<int>(nx) + x_extension_layers;

        for (int i = x_begin; i < x_end; ++i)
        {
            const Real x = (static_cast<Real>(i) + 0.5) * dx;
            for (size_t k = 0; k < nz; ++k)
            {
                const Real z = (static_cast<Real>(k) + 0.5) * dz;
                for (size_t layer = 0; layer < wall_layers; ++layer)
                {
                    const Real offset = (static_cast<Real>(layer) + 0.5) * cfg_.dp;
                    addPositionAndVolumetricMeasure(Vecd(x, -offset, z), volume);
                    addPositionAndVolumetricMeasure(Vecd(x, cfg_.DH + offset, z), volume);
                }
            }
        }
    }

  private:
    const cylinder_3d_compressible::CompressibleCylinderConfig &cfg_;
};

} // namespace SPH

#endif // CYLINDER_3D_COMPRESSIBLE_GEOMETRY_HPP
