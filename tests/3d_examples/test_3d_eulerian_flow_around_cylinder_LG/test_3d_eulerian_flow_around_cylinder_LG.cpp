#include "sphinxsys.h"
#include "eulerian_open_boundary.h"
#include "cylinder_3d_data.hpp"
#include "cylinder_3d_geometry.hpp"
#include "../shared/ck_time_step.hpp"

#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

using namespace SPH;
using namespace SPH::cylinder_3d;

struct SmokeRunResult
{
    bool geometry_gate = false;
    bool fifty_step_finite = false;
    bool final_finite = false;
    bool boundary_gate = false;
    bool periodic_gate = false;
    bool force_gate = false;
    bool dt_not_collapsed = false;
    size_t iterations = 0;
};

bool runContractTests()
{
    Cylinder3DConfig cfg = loadConfig(resolveDefaultConfigPath().string());
    bool pass = true;
    const auto expect_true = [&](bool condition, const std::string &name)
    {
        std::cout << "[Contract] " << name << " : " << (condition ? "PASS" : "FAIL") << std::endl;
        pass = pass && condition;
    };

    expect_true(std::fabs(cfg.c_f - cfg.sound_speed_factor * cfg.u_f) < 1.0e-12,
                "derived sound speed");
    expect_true(std::fabs(cfg.mu_f - cfg.rho0_f * cfg.u_f * cfg.D / cfg.re) < 1.0e-12,
                "derived dynamic viscosity");

    Cylinder3DConfig invalid_radius = cfg;
    invalid_radius.cylinder_radius = 0.0;
    expect_true(([&]()
                 {
                     try
                     {
                         deriveAndValidate(invalid_radius);
                         return false;
                     }
                     catch (const std::runtime_error &)
                     {
                         return true;
                     }
                 })(),
                "reject non-positive radius");

    Cylinder3DConfig invalid_pressure = cfg;
    invalid_pressure.outlet_pressure = 1.0;
    expect_true(([&]()
                 {
                     try
                     {
                         deriveAndValidate(invalid_pressure);
                         return false;
                     }
                     catch (const std::runtime_error &)
                     {
                         return true;
                     }
                 })(),
                "reject nonzero outlet pressure");

    Cylinder3DConfig invalid_mach = cfg;
    invalid_mach.sound_speed_factor = 5.0;
    expect_true(([&]()
                 {
                     try
                     {
                         deriveAndValidate(invalid_mach);
                         return false;
                     }
                     catch (const std::runtime_error &)
                     {
                         return true;
                     }
                 })(),
                "reject high Mach setup");

    Cylinder3DConfig invalid_position = cfg;
    invalid_position.cylinder_center_x = 0.05;
    expect_true(([&]()
                 {
                     try
                     {
                         deriveAndValidate(invalid_position);
                         return false;
                     }
                     catch (const std::runtime_error &)
                     {
                         return true;
                     }
                 })(),
                "reject cylinder touching boundary");

    return pass;
}

SmokeRunResult runCylinder3D(const Cylinder3DConfig &cfg, bool geometry_only = false)
{
    SmokeRunResult result;
    const Real BW = 4.0 * cfg.dp;
    BoundingBoxd system_domain_bounds(
        Vecd(-cfg.sponge_width - BW, -cfg.sponge_width - BW, -BW),
        Vecd(cfg.DL + cfg.sponge_width + BW, cfg.DH + cfg.sponge_width + BW, cfg.DW + BW));
    SPHSystem sph_system(system_domain_bounds, cfg.dp);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(false);

    FluidBody fluid_block(sph_system, makeShared<Cylinder3DFluidBlock>("Cylinder3DFluid", cfg));
    fluid_block.defineComponentLevelSetShape("OuterBoundary");
    fluid_block.defineMatterMaterial<WeaklyCompressibleFluid>(cfg.rho0_f, cfg.c_f);
    fluid_block.addMaterialProperty<Viscosity>(cfg.mu_f);
    fluid_block.generateParticles<BaseParticles, Lattice>();

    SolidBody cylinder_wall(sph_system, makeShared<Cylinder3DWallBlock>("Cylinder", cfg));
    cylinder_wall.defineAdaptationRatios(1.3, 1.0);
    cylinder_wall.defineBodyLevelSetShape();
    cylinder_wall.defineMatterMaterial<Solid>();
    cylinder_wall.generateParticles<BaseParticles, Lattice>();

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
        std::cout << "[Cylinder3D][Restart] enable=true step=" << restart_step
                  << " output_interval_steps="
                  << static_cast<size_t>(cfg.screen_output_interval) * static_cast<size_t>(cfg.restart_output_factor)
                  << " keep_last_n=" << cfg.restart_keep_last_n << std::endl;
    }

    InnerRelation fluid_inner(fluid_block);
    ContactRelation fluid_wall_contact(fluid_block, RealBodyVector{&cylinder_wall});
    ComplexRelation fluid_wall_complex(fluid_inner, fluid_wall_contact);
    InnerRelation cylinder_inner(cylinder_wall);
    ContactRelation cylinder_contact(cylinder_wall, RealBodyVector{&fluid_block});

    BoundingBoxd periodic_z_bounds(
        Vecd(-cfg.sponge_width, -cfg.sponge_width, 0.0),
        Vecd(cfg.DL + cfg.sponge_width, cfg.DH + cfg.sponge_width, cfg.DW));
    PeriodicAlongAxis periodic_along_z(periodic_z_bounds, zAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition_z(fluid_block, periodic_along_z);
    PeriodicConditionUsingCellLinkedList periodic_condition_wall_z(cylinder_wall, periodic_along_z);

    SimpleDynamics<Cylinder3DPeriodicWallNormal> cylinder_normal_direction(cylinder_wall, cfg);
    InteractionWithUpdate<FreeSurfaceIndicationComplex> surface_indicator(fluid_inner, fluid_wall_contact);
    InteractionDynamics<SmearedSurfaceIndication> smeared_surface(fluid_inner);
    SimpleDynamics<NormalDirectionFromBodyShape> fluid_normal_direction(fluid_block);

    Dynamics1Level<EulerianPressureRelaxationWithWallNormalOnly>
        pressure_relaxation(fluid_inner, fluid_wall_contact);
    InteractionWithUpdate<fluid_dynamics::EulerianIntegration2ndHalfWithWallRiemann>
        density_relaxation(fluid_inner, fluid_wall_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(fluid_inner, fluid_wall_contact);
    SimpleDynamics<Cylinder3DInitialCondition> initial_condition(fluid_block, cfg);
    channel_ck::AcousticTimeStep<> get_acoustic_dt(fluid_block, cfg.acoustic_cfl);
    InteractionWithUpdate<Cylinder3DFarFieldBoundary> farfield_boundary(fluid_inner, cfg);

    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(cylinder_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>>
        pressure_force_from_fluid(cylinder_contact);
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");

    sph_system.initializeSystemCellLinkedLists();
    periodic_condition_z.update_cell_linked_list_.exec();
    periodic_condition_wall_z.update_cell_linked_list_.exec();
    sph_system.initializeSystemConfigurations();
    cylinder_normal_direction.exec();
    if (cfg.enable_restart && restart_step > 0)
    {
        physical_time = restart_io->readRestartFiles(restart_step);
        periodic_condition_z.bounding_.exec();
        fluid_block.updateCellLinkedList();
        periodic_condition_z.update_cell_linked_list_.exec();
        cylinder_wall.updateCellLinkedList();
        periodic_condition_wall_z.update_cell_linked_list_.exec();
        fluid_wall_complex.updateConfiguration();
        cylinder_contact.updateConfiguration();
        cylinder_normal_direction.exec();
        std::cout << "[Cylinder3D][Restart] loaded step=" << restart_step
                  << " physical_time=" << physical_time << std::endl;
    }
    else
    {
        initial_condition.exec();
        farfield_boundary.exec();
    }
    surface_indicator.exec();
    smeared_surface.exec();
    fluid_normal_direction.exec();

    result.geometry_gate = reportGeometryGate(fluid_block, cylinder_wall, fluid_inner, fluid_wall_contact, cfg).pass;
    BaseParticles &fluid_particles = fluid_block.getBaseParticles();
    reportFluidFiniteState(fluid_particles, "after init");
    result.boundary_gate = reportBoundaryGate(fluid_particles, cfg, "after init").pass;
    if (geometry_only)
    {
        result.fifty_step_finite = true;
        result.final_finite = true;
        result.periodic_gate = reportPeriodicGate(fluid_particles, cfg, "after init").pass;
        result.force_gate = true;
        result.dt_not_collapsed = true;
        std::cout << "[Cylinder3D][Summary]"
                  << " geometry=" << result.geometry_gate
                  << " fifty_step_finite=" << result.fifty_step_finite
                  << " final_finite=" << result.final_finite
                  << " boundary=" << result.boundary_gate
                  << " periodic=" << result.periodic_gate
                  << " force=" << result.force_gate
                  << " dt_not_collapsed=" << result.dt_not_collapsed
                  << " iterations=" << result.iterations << std::endl;
        return result;
    }

    BodyStatesRecordingToVtp write_states(sph_system);
    write_states.addToWrite<int>(fluid_block, "Indicator");
    write_states.addToWrite<int>(fluid_block, "SmearedSurface");
    write_states.addToWrite<Real>(fluid_block, "Density");
    write_states.addToWrite<Real>(fluid_block, "Pressure");
    write_states.addToWrite<Vecd>(fluid_block, "Velocity");
    write_states.addToWrite<Vecd>(fluid_block, "NormalDirection");
    write_states.addToWrite<Vecd>(cylinder_wall, "NormalDirection");
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_viscous_force(cylinder_wall, "ViscousForceFromFluid");
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_pressure_force(cylinder_wall, "PressureForceFromFluid");
    ReducedQuantityRecording<MaximumSpeed> write_maximum_speed(fluid_block);

    const Real output_interval = cfg.end_time / static_cast<Real>(cfg.output_interval);
    const size_t restart_output_interval =
        static_cast<size_t>(cfg.screen_output_interval) * static_cast<size_t>(cfg.restart_output_factor);
    size_t iteration = restart_step;
    bool dt_collapsed = false;
    bool wrote_force = false;
    TickCount t1 = TickCount::now();
    TimeInterval interval;

    write_states.writeToFile(0);
    while (physical_time < cfg.end_time && !dt_collapsed)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval && physical_time < cfg.end_time && !dt_collapsed)
        {
            const Real dt = get_acoustic_dt.exec();
            if (dt < Real(1.0e-12) || !isFiniteReal(dt))
            {
                std::cout << "[Cylinder3D][DtCollapse] dt=" << dt
                          << " iteration=" << iteration << " time=" << physical_time << std::endl;
                reportFluidFiniteState(fluid_particles, "dt-collapse");
                dt_collapsed = true;
                break;
            }
            viscous_force.exec();
            pressure_relaxation.exec(dt);
            density_relaxation.exec(dt);
            farfield_boundary.exec();

            integration_time += dt;
            physical_time += dt;
            if (iteration == static_cast<size_t>(cfg.smoke_min_steps) && !result.fifty_step_finite)
            {
                result.fifty_step_finite = reportFluidFiniteState(fluid_particles, "50-step gate");
            }
            if (iteration % static_cast<size_t>(cfg.screen_output_interval) == 0)
            {
                write_maximum_speed.writeToFile(iteration);
                std::cout << std::fixed << std::setprecision(9)
                          << "N=" << iteration << " t=" << physical_time << " dt=" << dt << std::endl;
            }
            ++iteration;
            if (cfg.enable_restart && restart_output_interval > 0 &&
                iteration % restart_output_interval == 0 &&
                iteration != restart_step)
            {
                restart_io->writeToFile(iteration);
                cleanupOldRestartCheckpoints("restart", cfg.restart_keep_last_n);
            }
        }
        write_states.writeToFile();
        viscous_force_from_fluid.exec();
        pressure_force_from_fluid.exec();
        write_total_viscous_force.writeToFile(iteration);
        write_total_pressure_force.writeToFile(iteration);
        wrote_force = true;
    }

    TickCount t4 = TickCount::now();
    TimeInterval tt = t4 - t1 - interval;
    std::cout << "Total wall time: " << tt.seconds() << " s" << std::endl;

    result.iterations = iteration;
    if (!result.fifty_step_finite && iteration >= static_cast<size_t>(cfg.smoke_min_steps))
    {
        result.fifty_step_finite = reportFluidFiniteState(fluid_particles, "late 50-step gate");
    }
    result.final_finite = reportFluidFiniteState(fluid_particles, "final");
    result.boundary_gate = reportBoundaryGate(fluid_particles, cfg, "final").pass && result.boundary_gate;
    result.periodic_gate = reportPeriodicGate(fluid_particles, cfg, "final").pass;
    result.force_gate = wrote_force;
    result.dt_not_collapsed = !dt_collapsed;
    std::cout << "[Cylinder3D][Summary]"
              << " geometry=" << result.geometry_gate
              << " fifty_step_finite=" << result.fifty_step_finite
              << " final_finite=" << result.final_finite
              << " boundary=" << result.boundary_gate
              << " periodic=" << result.periodic_gate
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
        }
        if (contract_tests)
        {
            return runContractTests() ? 0 : 1;
        }

        Cylinder3DConfig cfg = loadConfig(resolveDefaultConfigPath().string());
        std::cout << "[Cylinder3D][Config] dp=" << cfg.dp
                  << " c_f=" << cfg.c_f
                  << " mu_f=" << cfg.mu_f
                  << " D=" << cfg.D
                  << " DW=" << cfg.DW << std::endl;
        SmokeRunResult result = runCylinder3D(cfg, geometry_only);
        const bool pass = result.geometry_gate && result.fifty_step_finite && result.final_finite &&
                          result.boundary_gate && result.periodic_gate && result.force_gate &&
                          result.dt_not_collapsed;
        return pass ? 0 : 1;
    }
    catch (const std::exception &e)
    {
        std::cerr << "[Cylinder3D][Error] " << e.what() << std::endl;
        return 1;
    }
}
