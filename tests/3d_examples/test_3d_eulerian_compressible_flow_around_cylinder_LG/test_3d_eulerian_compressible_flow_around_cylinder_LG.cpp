/**
 * @file test_3d_eulerian_compressible_flow_around_cylinder_LG.cpp
 * @brief 3D configurable-Mach, Re=100 ideal-gas fully compressible Eulerian SPH flow around
 *        a cylinder with configurable x far-field/periodic topology, z-periodic
 *        topology and configurable y far-field/wall faces.
 *
 * The fluid advances Mass / Momentum / TotalEnergy and recovers pressure from the
 * ideal-gas EOS. Inviscid fluxes reuse the shared MUSCL-HLLC inner + wall contact
 * integration; active open faces reuse the shared Eulerian ghost framework
 * with a fully compressible far-field state. In the x/z-periodic, y-wall
 * isolation mode no open-boundary ghost object exists at all.
 */
#include "sphinxsys.h"

#include "cylinder_3d_compressible_boundary.hpp"
#include "cylinder_3d_compressible_data.hpp"
#include "cylinder_3d_compressible_diagnostics.hpp"
#include "cylinder_3d_compressible_geometry.hpp"
#include "cylinder_3d_compressible_state.hpp"

#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

using namespace SPH;
using namespace SPH::cylinder_3d_compressible;

struct SmokeRunResult
{
    bool geometry_gate = false;
    bool open_mask_gate = false;
    bool ghost_map_gate = false;
    bool face_ghost_gate = false;
    bool wall_pair_gate = false;
    bool y_wall_gate = false;
    bool periodic_topology_gate = false;
    bool no_open_ghost_gate = false;
    bool wall_flux_gate = false;
    bool first_order_mask_gate = false;
    bool init_state_gate = false;
    bool fifty_step_finite = false;
    bool final_finite = false;
    bool conservation_gate = false;
    bool force_gate = false;
    bool dt_not_collapsed = false;
    size_t iterations = 0;
};

//----------------------------------------------------------------------
//  Configuration contract tests (--contract-tests).
//----------------------------------------------------------------------
namespace
{
bool expectThrow(const std::function<void()> &action)
{
    try
    {
        action();
        return false;
    }
    catch (const std::runtime_error &)
    {
        return true;
    }
}
} // namespace

bool runContractTests()
{
    bool pass = true;
    const auto expect_true = [&](bool condition, const std::string &name)
    {
        std::cout << "[Contract] " << name << " : " << (condition ? "PASS" : "FAIL") << std::endl;
        pass = pass && condition;
    };

    CompressibleCylinderConfig cfg = loadConfig(resolveDefaultConfigPath().string());
    printConfigSummary(cfg);

    // ---- derived freestream values ----
    const auto relative_error = [](Real value, Real reference)
    { return std::fabs(value - reference) / std::fabs(reference); };

    // Contract tests must stay config-agnostic: they pin the derivation and
    // parser logic, never the case's current operating point (mach, boundary
    // modes, stability switches). Anything that changes with config.ini
    // belongs to the run, not to this contract.
    expect_true(relative_error(cfg.u_inf, cfg.mach_inf * cfg.c_inf) < 1.0e-12, "u_inf = Ma*c_inf");
    expect_true(relative_error(cfg.p_inf, cfg.rho_inf * cfg.c_inf * cfg.c_inf / cfg.gamma) < 1.0e-12,
                "p_inf = rho*c^2/gamma");
    expect_true(relative_error(cfg.mu_inf, cfg.rho_inf * cfg.u_inf * cfg.D / cfg.re) < 1.0e-12,
                "mu = rho*u*D/Re");
    expect_true(relative_error(cfg.D, 2.0 * cfg.cylinder_radius) < 1.0e-12, "D = 2*radius");
    expect_true(relative_error(cfg.E_inf_per_volume,
                               cfg.p_inf / (cfg.gamma - 1.0) + 0.5 * cfg.rho_inf * cfg.u_inf * cfg.u_inf) < 1.0e-12,
                "E_inf per volume matches the ideal-gas definition");

    const CompressibleFreestreamState freestream = makeFreestreamState(cfg);
    expect_true(relative_error(freestream.vel[0], cfg.u_inf) < 1.0e-12 &&
                    std::fabs(freestream.vel[1]) < 1.0e-14 && std::fabs(freestream.vel[2]) < 1.0e-14,
                "freestream velocity is purely streamwise");
    expect_true(relative_error(std::sqrt(cfg.gamma * cfg.p_inf / cfg.rho_inf), cfg.c_inf) < 1.0e-12,
                "sound speed from p_inf recovers c_inf");
    expect_true(relative_error(cfg.u_inf / cfg.c_inf, cfg.mach_inf) < 1.0e-12, "Ma = u_inf / c_inf");
    expect_true(parseXBoundaryMode("farfield") == XBoundaryMode::FarField &&
                    parseXBoundaryMode("periodic") == XBoundaryMode::Periodic,
                "x boundary parser accepts both modes");
    expect_true(expectThrow([]
                            { parseXBoundaryMode("reflective"); }),
                "reject unsupported x boundary mode");
    expect_true(parseYBoundaryMode("farfield") == YBoundaryMode::FarField &&
                    parseYBoundaryMode("wall") == YBoundaryMode::Wall,
                "y boundary parser accepts both modes");
    expect_true(expectThrow([]
                            { parseYBoundaryMode("reflective"); }),
                "reject unsupported y boundary mode");

    // ---- rejection cases ----
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.mach_inf = 0.0;
                                deriveAndValidate(bad); }),
                "reject non-positive mach_inf");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.mach_inf = std::numeric_limits<Real>::infinity();
                                deriveAndValidate(bad); }),
                "reject non-finite mach_inf");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.gamma = 1.667;
                                deriveAndValidate(bad); }),
                "reject gamma != 1.4");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.rho_inf = 1.2;
                                deriveAndValidate(bad); }),
                "reject rho_inf != 1");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.c_inf = 2.0;
                                deriveAndValidate(bad); }),
                "reject c_inf != 1");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.global_resolution = 0.0;
                                deriveAndValidate(bad); }),
                "reject non-positive resolution");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.global_resolution = 3.0;
                                deriveAndValidate(bad); }),
                "reject a resolution that yields zero lattice cells");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.re = 0.0;
                                deriveAndValidate(bad); }),
                "reject non-positive Re");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.cylinder_center_x = 0.05;
                                deriveAndValidate(bad); }),
                "reject cylinder touching the boundary");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.acoustic_cfl = 0.0;
                                deriveAndValidate(bad); }),
                "reject non-positive acoustic CFL");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.y_wall_first_order_reconstruction = true;
                                bad.outlet_y_wall_corner_first_order_reconstruction = false;
                                bad.wall_contact_first_order_reconstruction = false;
                                bad.y_wall_first_order_band_dp = 0.0;
                                deriveAndValidate(bad); }),
                "reject non-positive y-wall first-order A/B band");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.y_wall_first_order_reconstruction = true;
                                bad.outlet_y_wall_corner_first_order_reconstruction = false;
                                bad.wall_contact_first_order_reconstruction = false;
                                bad.y_boundary_mode = YBoundaryMode::FarField;
                                deriveAndValidate(bad); }),
                "reject y-wall first-order A/B without y walls");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.wall_contact_first_order_reconstruction = false;
                                bad.outlet_y_wall_corner_first_order_reconstruction = true;
                                bad.outlet_corner_first_order_band_dp = 0.0;
                                deriveAndValidate(bad); }),
                "reject non-positive outlet-corner first-order A/B band");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.wall_contact_first_order_reconstruction = false;
                                bad.outlet_y_wall_corner_first_order_reconstruction = true;
                                bad.y_wall_first_order_band_dp = 0.0;
                                deriveAndValidate(bad); }),
                "reject non-positive y-wall band in outlet-corner first-order A/B mode");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.wall_contact_first_order_reconstruction = false;
                                bad.outlet_y_wall_corner_first_order_reconstruction = true;
                                bad.outlet_corner_first_order_band_dp = 0.5 * bad.DL / bad.dp;
                                deriveAndValidate(bad); }),
                "reject an outlet-corner x band that reaches the upstream half-domain");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.wall_contact_first_order_reconstruction = false;
                                bad.outlet_y_wall_corner_first_order_reconstruction = true;
                                bad.y_wall_first_order_band_dp = 0.5 * bad.DH / bad.dp;
                                deriveAndValidate(bad); }),
                "reject outlet-corner y bands that merge across the channel mid-plane");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.wall_contact_first_order_reconstruction = false;
                                bad.outlet_y_wall_corner_first_order_reconstruction = true;
                                bad.x_boundary_mode = XBoundaryMode::Periodic;
                                deriveAndValidate(bad); }),
                "reject outlet-corner first-order A/B without an x far-field outlet");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.wall_contact_first_order_reconstruction = false;
                                bad.outlet_y_wall_corner_first_order_reconstruction = true;
                                bad.y_boundary_mode = YBoundaryMode::FarField;
                                deriveAndValidate(bad); }),
                "reject outlet-corner first-order A/B without y walls");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.wall_contact_first_order_reconstruction = false;
                                bad.y_wall_first_order_reconstruction = true;
                                bad.outlet_y_wall_corner_first_order_reconstruction = true;
                                deriveAndValidate(bad); }),
                "reject simultaneous full-y-wall and outlet-corner first-order A/B modes");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.y_wall_first_order_reconstruction = true;
                                deriveAndValidate(bad); }),
                "reject simultaneous full-y-wall and wall-contact first-order A/B modes");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.outlet_y_wall_corner_first_order_reconstruction = true;
                                deriveAndValidate(bad); }),
                "reject simultaneous outlet-corner and wall-contact first-order A/B modes");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.DW = 0.0;
                                deriveAndValidate(bad); }),
                "reject non-positive spanwise period");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.x_boundary_mode = XBoundaryMode::Periodic;
                                bad.y_boundary_mode = YBoundaryMode::FarField;
                                deriveAndValidate(bad); }),
                "reject x periodic with y farfield");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.muscl_hllc_dissipation_limiter = true;
                                bad.muscl_hllc_limiter_parameter = 0.0;
                                deriveAndValidate(bad); }),
                "reject a non-positive HLLC limiter parameter when the limiter is enabled");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.muscl_hllc_dissipation_limiter = true;
                                bad.muscl_hllc_limiter_parameter =
                                    std::numeric_limits<Real>::infinity();
                                deriveAndValidate(bad); }),
                "reject a non-finite HLLC limiter parameter when the limiter is enabled");
    expect_true(([&]()
                 {
                     // Symmetric contract: the parameter is only read on the
                     // enabled path, so a disabled limiter must accept any
                     // value. Locks the conditional-validation design against
                     // a future unconditional check.
                     try
                     {
                         CompressibleCylinderConfig relaxed = cfg;
                         relaxed.muscl_hllc_dissipation_limiter = false;
                         relaxed.muscl_hllc_limiter_parameter = -1.0;
                         deriveAndValidate(relaxed);
                         return true;
                     }
                     catch (const std::runtime_error &)
                     {
                         return false;
                     } })(),
                "accept an arbitrary HLLC limiter parameter when the limiter is disabled");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.positivity_floor = true;
                                bad.positivity_floor_factor = 0.0;
                                deriveAndValidate(bad); }),
                "reject a zero positivity floor factor when the floor is enabled");
    expect_true(expectThrow([&]
                            {
                                CompressibleCylinderConfig bad = cfg;
                                bad.positivity_floor = true;
                                bad.positivity_floor_factor = 1.0;
                                deriveAndValidate(bad); }),
                "reject a positivity floor factor >= 1 when the floor is enabled");
    expect_true(([&]()
                 {
                     // Symmetric contract: a disabled floor must accept any
                     // factor; the value is only read on the enabled path.
                     try
                     {
                         CompressibleCylinderConfig relaxed = cfg;
                         relaxed.positivity_floor = false;
                         relaxed.positivity_floor_factor = -1.0;
                         deriveAndValidate(relaxed);
                         return true;
                     }
                     catch (const std::runtime_error &)
                     {
                         return false;
                     } })(),
                "accept an arbitrary positivity floor factor when the floor is disabled");

    // ---- legacy weakly-compressible keys must be rejected by loadConfig ----
    const std::filesystem::path temp_dir =
        std::filesystem::path(__FILE__).parent_path() / "tmp_contract";
    std::filesystem::create_directories(temp_dir);
    const auto write_legacy_config = [&](const std::string &extra_section_line,
                                         const std::string &extra_key_line)
    {
        const std::filesystem::path path = temp_dir / "legacy_config.ini";
        std::ofstream out(path, std::ios::out | std::ios::trunc);
        out << "[geometry]\ndl = 2.0\ndh = 1.0\ndw = 0.24\nglobal_resolution = 0.04\n"
            << "cylinder_center_x = 0.5\ncylinder_center_y = 0.5\ncylinder_radius = 0.10\n";
        if (extra_section_line == "geometry")
        {
            out << extra_key_line << "\n";
        }
        out << "[physical]\ngamma = 1.4\nrho_inf = 1.0\nc_inf = 1.0\nmach_inf = 0.3\nre = 100.0\n";
        if (extra_section_line == "physical")
        {
            out << extra_key_line << "\n";
        }
        out << "[simulation]\nend_time = 50.0\noutput_interval = 50\nscreen_output_interval = 1000\n"
            << "acoustic_cfl = 0.1\nboundary_n_layers = 3\n";
        out.close();
        return path;
    };

    expect_true(expectThrow([&]
                            { loadConfig(write_legacy_config("physical", "sound_speed_factor = 12.0").string()); }),
                "reject legacy sound_speed_factor");
    expect_true(expectThrow([&]
                            { loadConfig(write_legacy_config("physical", "outlet_pressure = 0.0").string()); }),
                "reject legacy gauge outlet_pressure");
    expect_true(expectThrow([&]
                            { loadConfig(write_legacy_config("geometry", "sponge_width_factor = 3.0").string()); }),
                "reject legacy sponge_width_factor");
    expect_true(([&]()
                 {
                     // The stripped-down config without legacy keys must still load.
                     try
                     {
                         loadConfig(write_legacy_config("none", "").string());
                         return true;
                     }
                     catch (const std::runtime_error &)
                     {
                         return false;
                     } })(),
                "accept a minimal compressible config");
    std::error_code ec;
    std::filesystem::remove_all(temp_dir, ec);

    return pass;
}

//----------------------------------------------------------------------
//  Main simulation / smoke run.
//----------------------------------------------------------------------
/**
 * @param geometry_only stop after the geometry / ghost gates.
 * @param smoke_only    stop right after the 50-step gate; this is interface
 *                      verification, never a physically resolved run.
 */
SmokeRunResult runCylinder3D(const CompressibleCylinderConfig &cfg, bool geometry_only = false,
                             bool smoke_only = false)
{
    SmokeRunResult result;
    // The freestream state is derived where it is consumed (initial condition,
    // ghost boundary condition), so it is deliberately not held here: GCC's
    // -Werror=unused-but-set-variable rejects a set-but-unread local.
    const Real BW = 4.0 * cfg.dp;
    BoundingBoxd system_domain_bounds(
        Vecd(-BW, -BW, -BW),
        Vecd(cfg.DL + BW, cfg.DH + BW, cfg.DW + BW));
    SPHSystem sph_system(system_domain_bounds, cfg.dp);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(false);

    //----------------------------------------------------------------------
    //  Bodies. The fluid uses reserve-aware generation so the shared ghost
    //  framework has capacity for all open-face owners.
    //----------------------------------------------------------------------
    FluidBody fluid_block(sph_system, makeShared<Cylinder3DFluidBlock>("Cylinder3DFluid", cfg));
    // No level set on the fluid. defineComponentLevelSetShape("OuterBoundary")
    // used to be here, but it does NOT replace initial_shape_ -- it only hangs a
    // level set on that one sub-shape (base_body.hpp:23-27), so the boolean tree
    // stayed on the query path and the tabulation itself was never read. Lattice
    // generation uses checkContain on the ComplexShape, and the normals are now
    // analytic (Cylinder3DFluidNormal), so nothing consumes it. Building it at
    // dp=0.02 is pure cost.
    fluid_block.defineMatterMaterial<CompressibleFluid>(cfg.gamma);
    fluid_block.addMaterialProperty<Viscosity>(cfg.mu_inf);
    Ghost<ReserveSizeFactor> ghost_boundary(cfg.ghost_reserve_factor);
    fluid_block.generateParticlesWithReserve<BaseParticles, PeriodicCylinderFluidLattice>(ghost_boundary, cfg);

    SolidBody cylinder_wall(sph_system, makeShared<Cylinder3DWallBlock>("Cylinder", cfg));
    cylinder_wall.defineAdaptationRatios(1.3, 1.0);
    cylinder_wall.defineBodyLevelSetShape();
    cylinder_wall.defineMatterMaterial<Solid>();
    cylinder_wall.generateParticles<BaseParticles, PeriodicCylinderWallLattice>(cfg);

    const bool use_x_periodic = cfg.x_boundary_mode == XBoundaryMode::Periodic;
    const bool use_y_walls = cfg.y_boundary_mode == YBoundaryMode::Wall;
    const bool has_open_faces = hasActiveOpenFace(cfg);
    SolidBody y_walls(sph_system, makeShared<DefaultShape>("YWalls"));
    y_walls.defineMatterMaterial<Solid>();
    if (use_y_walls)
    {
        y_walls.generateParticles<BaseParticles, YWallVolumeLattice>(cfg);
    }

    BaseParticles &fluid_particles = fluid_block.getBaseParticles();
    BaseParticles &wall_particles = cylinder_wall.getBaseParticles();
    BaseParticles &y_wall_particles = y_walls.getBaseParticles();

    // The six evolving flow fields must be registered before RestartIO is built.
    CompressibleConservativeStateView state = makeCompressibleStateView(fluid_particles, cfg.gamma);
    registerCompressibleEvolvingVariables(fluid_particles);

    size_t restart_step = 0;
    if (cfg.enable_restart)
    {
        if (cfg.restart_step == -1)
        {
            const int detected_step = detectLatestRestartStep();
            if (detected_step <= 0)
            {
                throw std::runtime_error("enable_restart=true and restart_step=-1, but no restart checkpoint was found.");
            }
            restart_step = static_cast<size_t>(detected_step);
        }
        else if (cfg.restart_step > 0)
        {
            restart_step = static_cast<size_t>(cfg.restart_step);
        }
    }
    sph_system.setRestartStep(restart_step);
    std::unique_ptr<RestartIO> restart_io;
    if (cfg.enable_restart)
    {
        restart_io = std::make_unique<RestartIO>(sph_system);
        std::cout << "[Cylinder3DCompressible][Restart] enable=true step=" << restart_step
                  << " output_interval_steps="
                  << static_cast<size_t>(cfg.screen_output_interval) * static_cast<size_t>(cfg.restart_output_factor)
                  << " keep_last_n=" << cfg.restart_keep_last_n << std::endl;
    }

    //----------------------------------------------------------------------
    //  Relations.
    //----------------------------------------------------------------------
    InnerRelation fluid_inner(fluid_block);
    RealBodyVector fluid_wall_bodies{&cylinder_wall};
    if (use_y_walls)
    {
        fluid_wall_bodies.push_back(&y_walls);
    }
    ContactRelation fluid_wall_contact(fluid_block, fluid_wall_bodies);
    ComplexRelation fluid_wall_complex(fluid_inner, fluid_wall_contact);
    std::unique_ptr<ContactRelation> y_wall_probe_contact;
    if (use_y_walls)
    {
        y_wall_probe_contact = std::make_unique<ContactRelation>(fluid_block, RealBodyVector{&y_walls});
    }
    ContactRelation cylinder_contact(cylinder_wall, RealBodyVector{&fluid_block});

    // The target is an infinite spanwise cylinder. z-periodic CLL closes the
    // spanwise fluid/cylinder-wall stencil. In the isolation mode x-periodic CLL
    // also closes both the fluid seam and the upper/lower y-wall contact seam.
    const BoundingBoxd periodic_x_bounds(
        Vecd(0.0, 0.0, 0.0), Vecd(cfg.DL, cfg.DH, cfg.DW));
    const BoundingBoxd periodic_z_bounds(
        Vecd(0.0, 0.0, 0.0), Vecd(cfg.DL, cfg.DH, cfg.DW));
    PeriodicAlongAxis periodic_along_x(periodic_x_bounds, xAxis);
    PeriodicAlongAxis periodic_along_z(periodic_z_bounds, zAxis);
    std::unique_ptr<PeriodicConditionUsingCellLinkedList> periodic_condition_x;
    std::unique_ptr<PeriodicConditionUsingCellLinkedList> periodic_condition_y_wall_x;
    if (use_x_periodic)
    {
        periodic_condition_x =
            std::make_unique<PeriodicConditionUsingCellLinkedList>(fluid_block, periodic_along_x);
        periodic_condition_y_wall_x =
            std::make_unique<PeriodicConditionUsingCellLinkedList>(y_walls, periodic_along_x);
    }
    PeriodicConditionUsingCellLinkedList periodic_condition_z(fluid_block, periodic_along_z);
    PeriodicConditionUsingCellLinkedList periodic_condition_wall_z(cylinder_wall, periodic_along_z);
    PeriodicConditionUsingCellLinkedList periodic_condition_y_wall_z(y_walls, periodic_along_z);

    //----------------------------------------------------------------------
    //  Dynamics. Construction order matters: kernel correction -> the three
    //  inner gradients -> the two MUSCL half steps.
    //----------------------------------------------------------------------
    SimpleDynamics<Cylinder3DRadialWallNormal> cylinder_normal_direction(cylinder_wall, cfg);
    // Analytic, not NormalDirectionFromBodyShape: the shared version does two
    // boolean-tree queries per particle against box-minus-triangle-mesh, which at
    // dp=0.02 did not finish in 37 minutes on 256 cores. See Cylinder3DFluidNormal.
    SimpleDynamics<Cylinder3DFluidNormal> fluid_normal_direction(fluid_block, cfg);
    std::unique_ptr<SimpleDynamics<YWallNormal>> y_wall_normal_direction;
    if (use_y_walls)
    {
        y_wall_normal_direction = std::make_unique<SimpleDynamics<YWallNormal>>(y_walls);
    }
    SimpleDynamics<Cylinder3DCompressibleInitialCondition> initial_condition(fluid_block, cfg);

    InteractionWithUpdate<FreeSurfaceIndicationComplex> surface_indicator(fluid_inner, fluid_wall_contact);
    SimpleDynamics<OpenBoundaryIndicatorMask> open_boundary_mask(fluid_block, cfg);

    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> kernel_correction_matrix(fluid_inner);
    InteractionDynamics<KernelGradientCorrectionInner> kernel_gradient_update(fluid_inner);
    InteractionWithUpdate<GhostKernelGradientUpdate> ghost_kernel_gradient_update(fluid_inner);

    InteractionWithUpdate<fluid_dynamics::DensityGradient<Inner<LinearGradientCorrection>>>
        density_gradient(fluid_inner);
    InteractionWithUpdate<fluid_dynamics::VelocityGradient<Inner<LinearGradientCorrection>>>
        velocity_gradient(fluid_inner);
    InteractionWithUpdate<fluid_dynamics::PressureGradient<Inner<LinearGradientCorrection>>>
        pressure_gradient(fluid_inner);
    std::unique_ptr<SimpleDynamics<YWallFirstOrderReconstruction>> y_wall_first_order_reconstruction;
    if (cfg.y_wall_first_order_reconstruction)
    {
        y_wall_first_order_reconstruction =
            std::make_unique<SimpleDynamics<YWallFirstOrderReconstruction>>(fluid_block, cfg);
    }
    std::unique_ptr<SimpleDynamics<OutletYWallCornerFirstOrderReconstruction>>
        outlet_y_wall_corner_first_order_reconstruction;
    if (cfg.outlet_y_wall_corner_first_order_reconstruction)
    {
        outlet_y_wall_corner_first_order_reconstruction = std::make_unique<
            SimpleDynamics<OutletYWallCornerFirstOrderReconstruction>>(fluid_block, cfg);
    }
    std::unique_ptr<SimpleDynamics<WallContactFirstOrderReconstruction>>
        wall_contact_first_order_reconstruction;
    if (cfg.wall_contact_first_order_reconstruction)
    {
        wall_contact_first_order_reconstruction =
            std::make_unique<SimpleDynamics<WallContactFirstOrderReconstruction>>(fluid_wall_contact);
    }
    const bool has_first_order_reconstruction =
        y_wall_first_order_reconstruction || outlet_y_wall_corner_first_order_reconstruction ||
        wall_contact_first_order_reconstruction;
    const auto exec_first_order_reconstruction = [&]()
    {
        if (y_wall_first_order_reconstruction)
        {
            y_wall_first_order_reconstruction->exec();
        }
        if (outlet_y_wall_corner_first_order_reconstruction)
        {
            outlet_y_wall_corner_first_order_reconstruction->exec();
        }
        if (wall_contact_first_order_reconstruction)
        {
            wall_contact_first_order_reconstruction->exec();
        }
    };
    result.first_order_mask_gate = !has_first_order_reconstruction;

    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(fluid_inner, fluid_wall_contact);
    // MUSCL-HLLC bridge configuration. Two separate configs:
    //
    //   muscl_bridge_cfg (inner): the dissipation-limited HLLC
    //   (HLLCWithLimiterRiemannSolver) is used for the inner interactions.
    //   Its Roe-averaged widened wave speeds prevent the bow-shock MUSCL
    //   undershoot that the plain Davis-speed HLLC leaves unchecked at Ma=2.
    //
    //   wall_bridge_cfg (wall contact): the PLAIN HLLC is always used for
    //   wall-contact interactions, regardless of the limiter switch. The
    //   limited solver's dissipation limiter (SMIN(param*SMAX((ul-ur)/clr,0),1))
    //   evaluates to zero when fluid flows toward the wall (ul < 0, ur = -ul),
    //   which suppresses the p_star momentum-flux correction p - rho*c*u_n.
    //   Without that correction the wall exerts no deceleration on the
    //   approaching fluid, and the no-slip constraint is lost -- the near-wall
    //   velocity stays at the freestream level. The plain HLLC applies the full
    //   p_star = p + rho*(s_l-ul)*(s_star-ul) and keeps the no-slip wall tight.
    //   This mirrors the 2D supersonic case, which has no wall contact at all.
    fluid_dynamics::MUSCLHLLCBridgeConfig muscl_bridge_cfg;
    muscl_bridge_cfg.use_hllc_dissipation_limiter = cfg.muscl_hllc_dissipation_limiter;
    muscl_bridge_cfg.hllc_limiter_parameter = cfg.muscl_hllc_limiter_parameter;
    fluid_dynamics::MUSCLHLLCBridgeConfig wall_bridge_cfg;
    wall_bridge_cfg.use_hllc_dissipation_limiter = false;
    InteractionWithUpdate<fluid_dynamics::EulerianCompressibleIntegration1stHalfMUSCLWithWall>
        momentum_relaxation(DynamicsArgs(fluid_inner, muscl_bridge_cfg),
                            DynamicsArgs(fluid_wall_contact, wall_bridge_cfg));
    InteractionWithUpdate<fluid_dynamics::EulerianCompressibleIntegration2ndHalfMUSCLWithWall>
        density_and_energy_relaxation(DynamicsArgs(fluid_inner, muscl_bridge_cfg),
                                      DynamicsArgs(fluid_wall_contact, wall_bridge_cfg));
    // Wall-contact-only second half step, used purely by the wall-flux contract
    // below. It shares the state arrays but never advances them (interaction
    // only writes the rate accumulators; no update() is invoked).
    InteractionDynamics<fluid_dynamics::EulerianCompressibleIntegration2ndHalfMUSCL<Contact<Wall>>>
        wall_flux_only(fluid_wall_contact, wall_bridge_cfg);
    std::unique_ptr<InteractionDynamics<fluid_dynamics::EulerianCompressibleIntegration2ndHalfMUSCL<Contact<Wall>>>>
        y_wall_flux_only;
    if (y_wall_probe_contact)
    {
        y_wall_flux_only = std::make_unique<InteractionDynamics<
            fluid_dynamics::EulerianCompressibleIntegration2ndHalfMUSCL<Contact<Wall>>>>(
            *y_wall_probe_contact, wall_bridge_cfg);
    }
    ReduceDynamics<fluid_dynamics::EulerianCompressibleAcousticTimeStepSize>
        get_acoustic_dt(fluid_block, cfg.acoustic_cfl);
    // Positivity floor: runs once per step after the conservative update
    // (and after the rate accounting, so the floor energy stays visible in
    // the ungated conservation imbalance instead of hiding in the rates).
    std::unique_ptr<SimpleDynamics<Cylinder3DPositivityFloor>> positivity_floor;
    if (cfg.positivity_floor)
    {
        positivity_floor =
            std::make_unique<SimpleDynamics<Cylinder3DPositivityFloor>>(fluid_block, cfg);
    }

    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(cylinder_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_and_energy_relaxation)>>
        pressure_force_from_fluid(cylinder_contact);
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");

    //----------------------------------------------------------------------
    //  Topology, restart and initial state. Both the fresh-start and the
    //  restart path rebuild the real-particle topology before the single
    //  ghost construction below.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    if (periodic_condition_x)
    {
        periodic_condition_x->update_cell_linked_list_.exec();
    }
    periodic_condition_z.update_cell_linked_list_.exec();
    periodic_condition_wall_z.update_cell_linked_list_.exec();
    if (use_y_walls)
    {
        if (periodic_condition_y_wall_x)
        {
            periodic_condition_y_wall_x->update_cell_linked_list_.exec();
        }
        periodic_condition_y_wall_z.update_cell_linked_list_.exec();
    }
    sph_system.initializeSystemConfigurations();
    cylinder_normal_direction.exec();
    if (y_wall_normal_direction)
    {
        y_wall_normal_direction->exec();
    }

    if (cfg.enable_restart && restart_step > 0)
    {
        physical_time = restart_io->readRestartFiles(restart_step);
        if (periodic_condition_x)
        {
            periodic_condition_x->bounding_.exec();
        }
        periodic_condition_z.bounding_.exec();
        fluid_block.updateCellLinkedList();
        if (periodic_condition_x)
        {
            periodic_condition_x->update_cell_linked_list_.exec();
        }
        periodic_condition_z.update_cell_linked_list_.exec();
        periodic_condition_wall_z.bounding_.exec();
        cylinder_wall.updateCellLinkedList();
        periodic_condition_wall_z.update_cell_linked_list_.exec();
        if (use_y_walls)
        {
            if (periodic_condition_y_wall_x)
            {
                periodic_condition_y_wall_x->bounding_.exec();
            }
            periodic_condition_y_wall_z.bounding_.exec();
            y_walls.updateCellLinkedList();
            if (periodic_condition_y_wall_x)
            {
                periodic_condition_y_wall_x->update_cell_linked_list_.exec();
            }
            periodic_condition_y_wall_z.update_cell_linked_list_.exec();
        }
        fluid_wall_complex.updateConfiguration();
        if (y_wall_probe_contact)
        {
            y_wall_probe_contact->updateConfiguration();
        }
        cylinder_contact.updateConfiguration();
        cylinder_normal_direction.exec();
        if (y_wall_normal_direction)
        {
            y_wall_normal_direction->exec();
        }
        std::cout << "[Cylinder3DCompressible][Restart] loaded step=" << restart_step
                  << " physical_time=" << physical_time << std::endl;
        // Force / ForcePrior and the change rates are recomputed every step and
        // must not be inherited from the checkpoint.
        Vecd *force = fluid_particles.getVariableDataByName<Vecd>("Force");
        Vecd *force_prior = fluid_particles.getVariableDataByName<Vecd>("ForcePrior");
        Real *dmass_dt = fluid_particles.getVariableDataByName<Real>("MassChangeRate");
        Real *dE_dt = fluid_particles.getVariableDataByName<Real>("TotalEnergyChangeRate");
        for (size_t i = 0; i != fluid_particles.TotalRealParticles(); ++i)
        {
            force[i] = Vecd::Zero();
            force_prior[i] = Vecd::Zero();
            dmass_dt[i] = 0.0;
            dE_dt[i] = 0.0;
        }
        // Never continue from an inconsistent checkpoint (plan Design Notes 1.3).
        if (!checkStateConsistency(state, 0, fluid_particles.TotalRealParticles(),
                                   "after restart read", cfg.u_inf)
                 .pass)
        {
            throw std::runtime_error(
                "Restart checkpoint failed the seven-field consistency gate; refusing to continue.");
        }
    }
    else
    {
        initial_condition.exec();
    }
    fluid_normal_direction.exec();

    result.init_state_gate =
        checkStateConsistency(state, 0, fluid_particles.TotalRealParticles(), "after init", cfg.u_inf).pass;
    result.geometry_gate =
        reportGeometryGate(fluid_block, cylinder_wall, fluid_inner, fluid_wall_contact, cfg).pass;
    result.wall_pair_gate =
        reportWallPairGate(fluid_particles, wall_particles, fluid_wall_contact, "after init").pass;
    result.y_wall_gate = !use_y_walls ||
                         reportWallPairGate(fluid_particles, y_wall_particles, fluid_wall_contact,
                                            "after init y walls", 1)
                             .pass;
    result.wall_pair_gate = result.wall_pair_gate && result.y_wall_gate;
    result.periodic_topology_gate =
        !use_x_periodic || reportPeriodicTopologyGate(fluid_block, fluid_inner, fluid_wall_contact,
                                                      y_wall_particles, cfg, "after init")
                               .pass;
    reportFluidFiniteState(fluid_particles, "after init");
    reportBulkFlow(state, fluid_particles, "after init");

    // Local particle count: the pre-generation estimate in deriveAndValidate()
    // works on the lattice grid, so the real count is checked here, where it
    // exists, before anything expensive runs.
    const size_t total_real_particles = fluid_particles.TotalRealParticles();
    std::cout << "[Cylinder3DCompressible][ParticleCount] fluid=" << total_real_particles
              << " wall=" << wall_particles.TotalRealParticles()
              << " y_wall=" << y_wall_particles.TotalRealParticles()
              << " limit=" << cfg.max_local_particles
              << " D/dp=" << cfg.D / cfg.dp << std::endl;
    if (total_real_particles > static_cast<size_t>(cfg.max_local_particles))
    {
        throw std::runtime_error("Actual fluid particle count exceeds max_local_particles; refusing local large case.");
    }

    //----------------------------------------------------------------------
    //  Open boundary: mask first, then construct ghosts only when an active open
    //  face exists. x=periodic/y=wall must not allocate any ghost object.
    //----------------------------------------------------------------------
    surface_indicator.exec();
    open_boundary_mask.exec();
    const OpenBoundaryMaskResult mask_result = reportOpenBoundaryMask(fluid_particles, cfg, "after mask");
    result.open_mask_gate = mask_result.pass;

    std::unique_ptr<GhostCreationInESPH> ghost_creation;
    std::unique_ptr<Cylinder3DCompressibleBoundaryCondition> farfield_boundary;
    std::unique_ptr<XFaceDirectedGhostCorrection> x_face_ghost_correction;
    if (has_open_faces)
    {
        // Fail with a message here rather than letting checkWithinGhostSize()
        // abort from inside the library when a finer dp overflows the reserve.
        if (mask_result.kept_indicator > ghost_boundary.getGhostSize())
        {
            throw std::runtime_error(
                "Ghost reserve is smaller than the number of open-face owners; raise ghost_reserve_factor.");
        }
        std::cout << "[Cylinder3DCompressible][GhostReserve] owners=" << mask_result.kept_indicator
                  << " reserve=" << ghost_boundary.getGhostSize() << std::endl;
        ghost_creation = std::make_unique<GhostCreationInESPH>(fluid_inner, ghost_boundary);
        farfield_boundary = std::make_unique<Cylinder3DCompressibleBoundaryCondition>(
            fluid_inner, *ghost_creation, cfg);
        if (use_y_walls)
        {
            x_face_ghost_correction = std::make_unique<XFaceDirectedGhostCorrection>(
                fluid_inner, *ghost_creation, cfg);
            result.face_ghost_gate = x_face_ghost_correction->exec().pass;
        }
        else
        {
            result.face_ghost_gate = true;
        }
        farfield_boundary->resetBoundaryConditions();
        result.ghost_map_gate = farfield_boundary->reportGhostMap("after first reset").pass;
        result.no_open_ghost_gate = true; // Not applicable: this is an open-boundary branch.
    }
    else
    {
        result.no_open_ghost_gate = mask_result.pass && mask_result.kept_indicator == 0;
        result.ghost_map_gate = result.no_open_ghost_gate;
        result.face_ghost_gate = true;
        std::cout << "[Cylinder3DCompressible][NoOpenGhostGate] active_owners="
                  << mask_result.kept_indicator << " ghost_creation=not-constructed"
                  << " pass=" << (result.no_open_ghost_gate ? "yes" : "NO") << std::endl;
    }

    kernel_correction_matrix.exec();
    kernel_gradient_update.exec();
    if (farfield_boundary)
    {
        ghost_kernel_gradient_update.exec();
    }
    if (x_face_ghost_correction)
    {
        result.face_ghost_gate = x_face_ghost_correction->exec(true).pass && result.face_ghost_gate;
    }

    // Wall-flux contract, run only now that the kernel correction matrix and the
    // ghost kernel gradients are in place -- otherwise the gate would measure a
    // different reconstruction than the production path. The rate accumulators
    // are zeroed, the wall contact half step runs alone (no update(), so no state
    // is advanced), and the imbalance is judged relative to the gross flux.
    {
        Real *dmass_dt = fluid_particles.getVariableDataByName<Real>("MassChangeRate");
        Real *dE_dt = fluid_particles.getVariableDataByName<Real>("TotalEnergyChangeRate");
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            dmass_dt[i] = 0.0;
            dE_dt[i] = 0.0;
        }
        density_gradient.exec();
        velocity_gradient.exec();
        pressure_gradient.exec();
        if (has_first_order_reconstruction)
        {
            exec_first_order_reconstruction();
            result.first_order_mask_gate =
                reportFirstOrderReconstructionMask(
                    fluid_particles, cfg, fluid_wall_contact, "after first reconstruction")
                    .pass;
        }
        wall_flux_only.exec(0.0);
        // The floor is a fraction of the reference mass flow through the cylinder
        // frontal area, so it scales with the case rather than being a bare epsilon.
        const Real reference_mass_flow = cfg.rho_inf * cfg.u_inf * cfg.D * cfg.DW;
        result.wall_flux_gate =
            reportWallFluxBalance(fluid_particles, fluid_wall_contact, "wall-contact only",
                                  cfg.wall_flux_relative_tolerance,
                                  cfg.wall_flux_absolute_floor_factor * reference_mass_flow)
                .pass;
        for (size_t i = 0; i != total_real_particles; ++i)
        {
            dmass_dt[i] = 0.0;
            dE_dt[i] = 0.0;
        }
    }

    if (geometry_only)
    {
        // Time-loop gates are not exercised here; they stay false and
        // --geometry-only reports only the gates it actually ran.
        std::cout << "[Cylinder3DCompressible][Summary] mode=geometry-only"
                  << " geometry=" << result.geometry_gate
                  << " open_mask=" << result.open_mask_gate
                  << " ghost_map=" << result.ghost_map_gate
                  << " face_ghost=" << result.face_ghost_gate
                  << " wall_pair=" << result.wall_pair_gate
                  << " y_wall=" << result.y_wall_gate
                  << " periodic_topology=" << result.periodic_topology_gate
                  << " no_open_ghost=" << result.no_open_ghost_gate
                  << " wall_flux=" << result.wall_flux_gate
                  << " first_order_mask=" << result.first_order_mask_gate
                  << " init_state=" << result.init_state_gate
                  << " (time-loop gates not evaluated)" << std::endl;
        return result;
    }

    //----------------------------------------------------------------------
    //  Output.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<int>(fluid_block, "Indicator");
    write_states.addToWrite<Real>(fluid_block, "Density");
    write_states.addToWrite<Real>(fluid_block, "Pressure");
    write_states.addToWrite<Real>(fluid_block, "TotalEnergy");
    write_states.addToWrite<Vecd>(fluid_block, "Velocity");
    write_states.addToWrite<Vecd>(fluid_block, "Momentum");
    write_states.addToWrite<Vecd>(fluid_block, "NormalDirection");
    if (has_first_order_reconstruction)
    {
        write_states.addToWrite<int>(fluid_block, "FirstOrderReconstructionMask");
    }
    write_states.addToWrite<Vecd>(cylinder_wall, "NormalDirection");
    if (use_y_walls)
    {
        write_states.addToWrite<Vecd>(y_walls, "NormalDirection");
    }
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_viscous_force(cylinder_wall, "ViscousForceFromFluid");
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_pressure_force(cylinder_wall, "PressureForceFromFluid");
    ReducedQuantityRecording<MaximumSpeed> write_maximum_speed(fluid_block);

    // Surface-closure check on the discrete cylinder, before any step is taken.
    //
    // At t = 0 the pressure field is exactly p_inf everywhere (the initial
    // condition perturbs only the streamwise velocity), so the continuous integral
    // of the pressure traction over a closed surface vanishes identically and the
    // discrete pressure load is pure closure residual of the lattice-generated
    // surface. This matters at Ma = 0.3 because the constant-pressure part enters
    // Cd amplified by p_inf/q_inf = 2/(gamma*Ma^2) = 15.9, so even an O(dp/R)
    // closure error would surface as an O(1) spurious drag that never decays.
    // Measured at velocity_noise_amplitude = 0: |Fp|/(p_inf*D*DW) = 1.8e-15, i.e.
    // the closure is exact to machine precision and the amplification is not
    // triggered.
    //
    // The measurement survives velocity_noise_amplitude > 0. The only other term
    // in PressureForceFromFluid is riemann_solvers_k.DissipativePJump(u_jump)
    // (fluid_structure_interaction.hpp:55), and this case's solver defines
    //     HLLCFSIRiemannSolver::DissipativePJump(Real) { return 0.0; }
    // (eulerian_riemann_solver.h:103) -- identically zero, independent of u_jump.
    // So the per-particle velocity scatter cannot enter the pressure load, and the
    // baseline stays a pure closure residual: measured 1.71e-14 at amplitude 0 and
    // bit-identical at 0.1. Re-check this if the solver type is ever changed to one
    // with a live pressure jump.
    //
    // The viscous part is NOT a closure residual: a freestream against a no-slip
    // wall is an impulsively started shear layer, whose t -> 0 traction is
    // physically singular. So Fv > 0 here is expected, and it does shift slightly
    // with the perturbation (0.591 -> 0.586) because the shear itself changes.
    viscous_force_from_fluid.exec();
    pressure_force_from_fluid.exec();
    const CylinderLoadReport uniform_flow_load =
        reportCylinderLoad(wall_particles, cfg, "freestream baseline (t=0, closure residual)");
    const Real closure_scale = cfg.p_inf * cfg.D * cfg.DW;
    std::cout << "[Cylinder3DCompressible][LoadBaseline]"
              << " p_inf_over_q_inf=" << 2.0 / (cfg.gamma * cfg.mach_inf * cfg.mach_inf)
              << " |Fp|/(p_inf*D*DW)=" << uniform_flow_load.pressure_force.norm() / closure_scale
              << " Cd_from_pressure=" << uniform_flow_load.pressure_force[0] / forceScale(cfg)
              << " Cd_from_viscous=" << uniform_flow_load.viscous_force[0] / forceScale(cfg)
              << " noise=" << cfg.velocity_noise_amplitude
              << " -- pressure part is the closure residual (perturbation-independent:"
              << " DissipativePJump is identically 0 for this solver);"
              << " viscous part is the impulsive-start shear and is physical."
              << std::endl;

    const Real output_interval = cfg.end_time / static_cast<Real>(cfg.output_interval);
    // In smoke mode the production interval (screen_output_interval *
    // restart_output_factor = 10000 steps by default) would never be reached
    // within the ~51-step budget, so the restart round-trip could not be
    // exercised at all. Write one checkpoint mid-run instead.
    const size_t restart_output_interval =
        smoke_only
            ? SMAX(size_t(1), static_cast<size_t>(cfg.smoke_min_steps) / 2)
            : static_cast<size_t>(cfg.screen_output_interval) * static_cast<size_t>(cfg.restart_output_factor);
    const ConservationBudget reference_budget =
        accumulateConservationBudget(state, 0, total_real_particles);
    // Running total of the applied change rates, for the report only. Default
    // constructed, not copied from reference_budget: only the accumulated_* fields
    // are ever read from it, and carrying stale totals here is a trap.
    ConservationBudget applied_rate_total;
    Real *dmass_dt_field = fluid_particles.getVariableDataByName<Real>("MassChangeRate");
    Real *dE_dt_field = fluid_particles.getVariableDataByName<Real>("TotalEnergyChangeRate");
    if (dmass_dt_field == nullptr || dE_dt_field == nullptr)
    {
        throw std::runtime_error("Conservation report needs MassChangeRate / TotalEnergyChangeRate.");
    }
    std::vector<Real> y_wall_probe_saved_mass_rate;
    std::vector<Real> y_wall_probe_saved_energy_rate;
    if (y_wall_flux_only)
    {
        y_wall_probe_saved_mass_rate.resize(total_real_particles);
        y_wall_probe_saved_energy_rate.resize(total_real_particles);
    }
    // Sampled throughout the run rather than once at step 50, so the totals'
    // trajectory is visible instead of a single early snapshot.
    const size_t conservation_report_interval =
        SMAX(size_t(1), static_cast<size_t>(cfg.smoke_min_steps));
    bool conservation_gate_ever_failed = false;
    // In smoke mode the run is bounded by the step count, not by physical time.
    const size_t smoke_step_limit = restart_step + static_cast<size_t>(cfg.smoke_min_steps) + 1;
    // Roughly 20 load samples over the smoke budget, whatever the budget is.
    const size_t load_trace_interval =
        SMAX(size_t(1), static_cast<size_t>(cfg.smoke_min_steps) / 20);
    size_t iteration = restart_step;
    bool dt_collapsed = false;
    bool wrote_force = false;
    // "Never evaluated" and "evaluated and failed" must stay distinguishable:
    // a bare !fifty_step_finite guard is true in both cases, so the late
    // fallback would re-evaluate a failed step-50 gate on the final state and
    // overwrite the failure with a pass.
    bool fifty_step_gate_evaluated = false;
    // Latched so a defect seen at one output step cannot be replaced by a clean
    // later sample. Same reason the conservation gate latches.
    bool state_gate_ever_failed = false;
    bool force_gate_ever_failed = false;
    bool smoke_limit_reached = false;
    TickCount t1 = TickCount::now();
    TimeInterval interval;

    write_states.writeToFile(0);
    while (physical_time < cfg.end_time && !dt_collapsed && !smoke_limit_reached)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval && physical_time < cfg.end_time &&
               !dt_collapsed && !smoke_limit_reached)
        {
            const Real dt = get_acoustic_dt.exec();
            if (dt < Real(1.0e-12) || !isFiniteRealValue(dt))
            {
                std::cout << "[Cylinder3DCompressible][DtCollapse] dt=" << dt
                          << " iteration=" << iteration << " time=" << physical_time << std::endl;
                reportFluidFiniteState(fluid_particles, "dt-collapse");
                // Evolved state: the write-path Momentum == Mass*Velocity
                // contract does not apply mid-scheme, see checkStateConsistency.
                checkStateConsistency(state, 0, total_real_particles, "dt-collapse", cfg.u_inf,
                                      1.0e-12, false);
                dt_collapsed = true;
                break;
            }

            // Time loop order fixed by Design Notes 1.7: the gradients are
            // refreshed after every state update, so the second half step
            // reconstructs from the post-first-half velocity rather than the
            // stale pre-first-half one.
            if (farfield_boundary)
            {
                farfield_boundary->resetBoundaryConditions();
            }
            density_gradient.exec();
            velocity_gradient.exec();
            pressure_gradient.exec();
            if (has_first_order_reconstruction)
            {
                exec_first_order_reconstruction();
            }
            if (farfield_boundary)
            {
                ghost_kernel_gradient_update.exec();
            }
            if (x_face_ghost_correction)
            {
                result.face_ghost_gate = x_face_ghost_correction->exec().pass && result.face_ghost_gate;
            }
            viscous_force.exec();
            momentum_relaxation.exec(dt);

            if (farfield_boundary)
            {
                farfield_boundary->resetBoundaryConditions();
            }
            density_gradient.exec();
            velocity_gradient.exec();
            pressure_gradient.exec();
            if (has_first_order_reconstruction)
            {
                exec_first_order_reconstruction();
            }
            if (farfield_boundary)
            {
                ghost_kernel_gradient_update.exec();
            }
            if (x_face_ghost_correction)
            {
                result.face_ghost_gate = x_face_ghost_correction->exec().pass && result.face_ghost_gate;
            }
            if (y_wall_flux_only && iteration % static_cast<size_t>(cfg.wall_flux_probe_interval) == 0)
            {
                // This contact-only execution shares rate arrays with production.
                // Snapshot, zero, sample and restore them so it cannot affect the
                // following combined second-half update.
                std::copy(dmass_dt_field, dmass_dt_field + total_real_particles,
                          y_wall_probe_saved_mass_rate.begin());
                std::copy(dE_dt_field, dE_dt_field + total_real_particles,
                          y_wall_probe_saved_energy_rate.begin());
                std::fill(dmass_dt_field, dmass_dt_field + total_real_particles, Real(0));
                std::fill(dE_dt_field, dE_dt_field + total_real_particles, Real(0));
                y_wall_flux_only->exec(0.0);
                const Real reference_mass_flow = cfg.rho_inf * cfg.u_inf * cfg.D * cfg.DW;
                reportWallFluxBalance(fluid_particles, *y_wall_probe_contact,
                                      "y-wall probe step " + std::to_string(iteration),
                                      cfg.wall_flux_relative_tolerance,
                                      cfg.wall_flux_absolute_floor_factor * reference_mass_flow);
                std::copy(y_wall_probe_saved_mass_rate.begin(), y_wall_probe_saved_mass_rate.end(),
                          dmass_dt_field);
                std::copy(y_wall_probe_saved_energy_rate.begin(), y_wall_probe_saved_energy_rate.end(),
                          dE_dt_field);
            }
            density_and_energy_relaxation.exec(dt);
            // Accumulate with the same dt that advanced the state, right after the
            // second half step wrote the rates and before the next step
            // overwrites them.
            accumulateAppliedRates(applied_rate_total, dmass_dt_field, dE_dt_field,
                                   0, total_real_particles, dt);
            // Floor after the rate accounting: the correction is deliberately
            // outside the solver's applied rates (see Cylinder3DPositivityFloor).
            if (positivity_floor)
            {
                positivity_floor->exec();
            }

            if (farfield_boundary)
            {
                farfield_boundary->resetBoundaryConditions();
            }

            integration_time += dt;
            physical_time += dt;
            // Counted from the run's own start so a restart still hits the gate.
            if (iteration - restart_step == static_cast<size_t>(cfg.smoke_min_steps) &&
                !fifty_step_gate_evaluated)
            {
                result.fifty_step_finite =
                    reportFluidFiniteState(fluid_particles, "50-step gate") &&
                    checkStateConsistency(state, 0, total_real_particles, "50-step gate", cfg.u_inf,
                                          1.0e-10, false)
                        .pass;
                fifty_step_gate_evaluated = true;
            }
            // Conservation report. Only its finiteness/positivity verdict
            // is gated, and that verdict is latched.
            if (iteration % conservation_report_interval == 0)
            {
                ConservationBudget sample = accumulateConservationBudget(state, 0, total_real_particles);
                sample.accumulated_mass_rate = applied_rate_total.accumulated_mass_rate;
                sample.accumulated_energy_rate = applied_rate_total.accumulated_energy_rate;
                if (!reportConservationDrift(reference_budget, sample, cfg,
                                             "step " + std::to_string(iteration))
                         .pass)
                {
                    conservation_gate_ever_failed = true;
                }
                reportBulkFlow(state, fluid_particles, "step " + std::to_string(iteration));
                result.conservation_gate = !conservation_gate_ever_failed;
            }
            // Load history during a bounded smoke run. Without it the load is
            // sampled once at the end, which cannot distinguish a decaying
            // impulsive-start transient from a converged force -- the only
            // signature that separates them is the trend in time.
            if (smoke_only && iteration % load_trace_interval == 0)
            {
                viscous_force_from_fluid.exec();
                pressure_force_from_fluid.exec();
                const CylinderLoadReport trace = reportCylinderLoad(wall_particles, cfg, "load trace");
                std::cout << "[Cylinder3DCompressible][LoadTrace] N=" << iteration
                          << " t=" << physical_time
                          << " t/(D/u_inf)=" << physical_time * cfg.u_inf / cfg.D
                          << " Cd=" << trace.cd
                          << " Cd_p=" << trace.pressure_force[0] / forceScale(cfg)
                          << " Cd_v=" << trace.viscous_force[0] / forceScale(cfg)
                          << " Cl=" << trace.cl << std::endl;
            }
            if (iteration % static_cast<size_t>(cfg.screen_output_interval) == 0)
            {
                write_maximum_speed.writeToFile(iteration);
                std::cout << std::fixed << std::setprecision(9)
                          << "N=" << iteration << " t=" << physical_time << " dt=" << dt << std::endl;
            }
            ++iteration;
            if (cfg.enable_restart && restart_output_interval > 0 &&
                iteration % restart_output_interval == 0 && iteration != restart_step)
            {
                restart_io->writeToFile(iteration);
                cleanupOldRestartCheckpoints("restart", cfg.restart_keep_last_n);
            }
            if (smoke_only && iteration >= smoke_step_limit)
            {
                std::cout << "[Cylinder3DCompressible][Smoke] reached step limit " << smoke_step_limit
                          << "; this is a numerical interface check only, no physical Cd/Cl or St claim."
                          << std::endl;
                smoke_limit_reached = true;
            }
        }
        // Compute the loads before writing, so the VTP frame carries this step's
        // cylinder forces rather than the previous step's.
        viscous_force_from_fluid.exec();
        pressure_force_from_fluid.exec();
        write_states.writeToFile();
        write_total_viscous_force.writeToFile(iteration);
        write_total_pressure_force.writeToFile(iteration);
        const CylinderLoadReport load = reportCylinderLoad(wall_particles, cfg, "output step");
        // The far-field problem must retain a positive drag direction. In this
        // unforced closed periodic isolation experiment, however, U_bulk decays
        // under wall/cylinder dissipation; a late force sign is therefore a
        // physical-interpretation warning, not evidence of numerical collapse.
        // Finiteness remains a hard gate in every boundary mode.
        if (!load.finite || (!use_x_periodic && !load.drag_direction_ok))
        {
            force_gate_ever_failed = true;
        }
        else if (use_x_periodic && !load.drag_direction_ok)
        {
            std::cout << "[Cylinder3DCompressible][DragDirectionWarning] t=" << physical_time
                      << " Cd=" << load.cd
                      << " -- unforced periodic isolation: inspect U_bulk before interpreting drag."
                      << std::endl;
        }
        result.force_gate = !force_gate_ever_failed;
        // The result was previously discarded, so a non-finite or non-positive
        // state found at an output step never reached the exit code.
        if (!checkStateConsistency(state, 0, total_real_particles, "output step", cfg.u_inf,
                                   1.0e-10, false)
                 .pass)
        {
            state_gate_ever_failed = true;
        }
        wrote_force = true;
    }

    TickCount t4 = TickCount::now();
    TimeInterval tt = t4 - t1 - interval;
    std::cout << "Total wall time: " << tt.seconds() << " s" << std::endl;

    result.iterations = iteration;
    // Late fallback for the case where the in-loop gate never fired (the loop
    // exited on physical time between two checks). It must apply the SAME two
    // criteria as the in-loop version: replacing the state-consistency term with
    // a finiteness-only check would turn a failed gate into a reported pass.
    if (!fifty_step_gate_evaluated &&
        iteration - restart_step >= static_cast<size_t>(cfg.smoke_min_steps))
    {
        result.fifty_step_finite =
            reportFluidFiniteState(fluid_particles, "late 50-step gate") &&
            checkStateConsistency(state, 0, total_real_particles, "late 50-step gate", cfg.u_inf,
                                  1.0e-10, false)
                .pass;
        fifty_step_gate_evaluated = true;
    }
    // Final conservation budget over the whole run. Combined with the latch so a
    // failure at any earlier sample still fails the gate.
    if (iteration > restart_step)
    {
        ConservationBudget final_sample = accumulateConservationBudget(state, 0, total_real_particles);
        final_sample.accumulated_mass_rate = applied_rate_total.accumulated_mass_rate;
        final_sample.accumulated_energy_rate = applied_rate_total.accumulated_energy_rate;
        if (!reportConservationDrift(reference_budget, final_sample, cfg, "final").pass)
        {
            conservation_gate_ever_failed = true;
        }
        reportBulkFlow(state, fluid_particles, "final");
        result.conservation_gate = !conservation_gate_ever_failed;
    }
    // Carries the latched per-output-step consistency verdict, so a defect seen
    // mid-run fails the run even if the final state happens to be clean.
    result.final_finite = reportFluidFiniteState(fluid_particles, "final") && !state_gate_ever_failed;
    if (farfield_boundary)
    {
        result.ghost_map_gate =
            farfield_boundary->reportGhostMap("final").pass && result.ghost_map_gate;
    }
    // force_gate is set inside the loop from the load report; require that the
    // loads were actually written at least once as well.
    result.force_gate = result.force_gate && wrote_force;
    result.dt_not_collapsed = !dt_collapsed;
    std::cout << "[Cylinder3DCompressible][Summary]"
              << " geometry=" << result.geometry_gate
              << " open_mask=" << result.open_mask_gate
              << " ghost_map=" << result.ghost_map_gate
              << " face_ghost=" << result.face_ghost_gate
              << " wall_pair=" << result.wall_pair_gate
              << " y_wall=" << result.y_wall_gate
              << " periodic_topology=" << result.periodic_topology_gate
              << " no_open_ghost=" << result.no_open_ghost_gate
              << " wall_flux=" << result.wall_flux_gate
              << " first_order_mask=" << result.first_order_mask_gate
              << " init_state=" << result.init_state_gate
              << " fifty_step_finite=" << result.fifty_step_finite
              << " final_finite=" << result.final_finite
              << " conservation=" << result.conservation_gate
              << " force=" << result.force_gate
              << " dt_not_collapsed=" << result.dt_not_collapsed
              << " iterations=" << result.iterations << std::endl;
    return result;
}

int main(int argc, char **argv)
{
    try
    {
        bool contract_tests = false;
        bool geometry_only = false;
        bool smoke = false;
        int smoke_steps_override = 0;
        for (int i = 1; i < argc; ++i)
        {
            const std::string arg(argv[i]);
            if (arg == "--contract-tests")
            {
                contract_tests = true;
            }
            else if (arg == "--geometry-only")
            {
                geometry_only = true;
            }
            else if (arg == "--smoke")
            {
                smoke = true;
            }
            // Step-count override for the bounded smoke run. The 50-step default
            // covers only ~0.1 D/u_inf, which is far too short to tell an
            // impulsive-start acoustic transient from a steady-state load; this
            // lets the transient be measured without touching config.ini or
            // starting an unbounded production run.
            else if (arg.rfind("--smoke-steps=", 0) == 0)
            {
                smoke_steps_override = std::stoi(arg.substr(std::string("--smoke-steps=").size()));
                if (smoke_steps_override <= 0)
                {
                    throw std::runtime_error("--smoke-steps must be positive.");
                }
            }
        }
        if (contract_tests)
        {
            return runContractTests() ? 0 : 1;
        }

        CompressibleCylinderConfig cfg = loadConfig(resolveDefaultConfigPath().string());
        if (smoke_steps_override > 0)
        {
            if (!smoke)
            {
                throw std::runtime_error("--smoke-steps only applies together with --smoke.");
            }
            std::cout << "[Cylinder3DCompressible][Smoke] step budget overridden: "
                      << cfg.smoke_min_steps << " -> " << smoke_steps_override << std::endl;
            cfg.smoke_min_steps = smoke_steps_override;
        }
        printConfigSummary(cfg);

        SmokeRunResult result = runCylinder3D(cfg, geometry_only, smoke);
        // Static gates are evaluated in every mode; the time-loop gates only
        // when the loop actually ran, so --geometry-only does not report
        // unevaluated gates as passing.
        const bool static_gates = result.geometry_gate && result.open_mask_gate &&
                                  result.ghost_map_gate && result.wall_pair_gate &&
                                  result.periodic_topology_gate && result.no_open_ghost_gate &&
                                  result.wall_flux_gate && result.first_order_mask_gate &&
                                  result.init_state_gate;
        const bool time_loop_gates = result.fifty_step_finite && result.final_finite &&
                                     result.conservation_gate && result.force_gate &&
                                     result.dt_not_collapsed;
        const bool pass = geometry_only ? static_gates : (static_gates && time_loop_gates);
        return pass ? 0 : 1;
    }
    catch (const std::exception &e)
    {
        std::cerr << "[Cylinder3DCompressible][Error] " << e.what() << std::endl;
        return 1;
    }
}
