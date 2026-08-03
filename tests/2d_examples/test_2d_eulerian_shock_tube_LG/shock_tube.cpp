/**
 * @file 	shock_tube.cpp
 * @brief 	2D Eulerian SPH shock tube (Lax problem), solid-wall paradigm.
 *          x ends are real SolidBody reflective walls (ContactRelation +
 *          compressible MUSCL-WithWall integrators); y direction is periodic.
 * @author 	KIYOYOZU
 */
#include "shock_tube.h"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace SPH;
//----------------------------------------------------------------------
//	CSV probe: dump all real particles (x, rho, p, u) for post-processing.
//----------------------------------------------------------------------
class ShockTubeCsvProbe
{
  public:
    ShockTubeCsvProbe(BaseParticles &particles)
        : particles_(particles),
          // Source dir is baked in by CMake as CASE_SOURCE_DIR so the probe
          // always writes to <case>/output regardless of the executable's cwd.
          output_dir_((std::filesystem::path(CASE_SOURCE_DIR) / "output").string()),
          pos_(particles.getVariableDataByName<Vecd>("Position")),
          rho_(particles.getVariableDataByName<Real>("Density")),
          p_(particles.getVariableDataByName<Real>("Pressure")),
          vel_(particles.getVariableDataByName<Vecd>("Velocity"))
    {
        std::filesystem::create_directories(output_dir_);
    }

    void write(size_t step, Real physical_time) const
    {
        std::ofstream ofs(output_dir_ + "/particles_" + std::to_string(step) + ".csv");
        ofs << "x,y,rho,p,u\n";
        ofs << std::setprecision(12);
        for (size_t i = 0; i != particles_.TotalRealParticles(); ++i)
        {
            ofs << pos_[i][0] << "," << pos_[i][1] << ","
                << rho_[i] << "," << p_[i] << "," << vel_[i][0] << "\n";
        }
        std::ofstream meta(output_dir_ + "/time_" + std::to_string(step) + ".txt");
        meta << physical_time << "\n";
        std::cout << "[Probe] wrote particles_" << step << ".csv (N=" << particles_.TotalRealParticles()
                  << ", t=" << physical_time << ")\n";
    }

  protected:
    BaseParticles &particles_;
    std::string output_dir_;
    Vecd *pos_;
    Real *rho_, *p_;
    Vecd *vel_;
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(false);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody fluid_block(sph_system, makeShared<FluidBlock>("FluidBlock"));
    fluid_block.defineComponentLevelSetShape("OuterBoundary");
    fluid_block.defineMatterMaterial<CompressibleFluid>(heat_capacity_ratio);
    fluid_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_left(sph_system, makeShared<WallBlock>("WallLeft", -wall_thickness, 0.0));
    wall_left.defineAdaptationRatios(1.3, 1.0);
    wall_left.defineBodyLevelSetShape();
    wall_left.defineMatterMaterial<Solid>();
    wall_left.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_right(sph_system, makeShared<WallBlock>("WallRight", L, L + wall_thickness));
    wall_right.defineAdaptationRatios(1.3, 1.0);
    wall_right.defineBodyLevelSetShape();
    wall_right.defineMatterMaterial<Solid>();
    wall_right.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //----------------------------------------------------------------------
    InnerRelation fluid_inner(fluid_block);
    ContactRelation fluid_wall_contact(fluid_block, RealBodyVector{&wall_left, &wall_right});
    //----------------------------------------------------------------------
    //	y-direction periodicity (x ends are solid walls).
    //----------------------------------------------------------------------
    BoundingBoxd periodic_bounds(Vec2d(0.0, 0.0), Vec2d(L, H));
    PeriodicAlongAxis periodic_along_y(periodic_bounds, yAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition_y(fluid_block, periodic_along_y);
    PeriodicConditionUsingCellLinkedList periodic_condition_wall_left_y(wall_left, periodic_along_y);
    PeriodicConditionUsingCellLinkedList periodic_condition_wall_right_y(wall_right, periodic_along_y);
    //----------------------------------------------------------------------
    //	Define the numerical methods. Construction order matters:
    // normals -> indicator -> gradient buffer registration
    // -> the two MUSCL half steps.
    //----------------------------------------------------------------------
    SimpleDynamics<WallNormal> wall_left_normal(wall_left, Vecd(1.0, 0.0));
    std::cout << "[dbg] wall_left_normal ok\n";
    SimpleDynamics<WallNormal> wall_right_normal(wall_right, Vecd(-1.0, 0.0));
    std::cout << "[dbg] wall_right_normal ok\n";
    SimpleDynamics<NormalDirectionFromBodyShape> fluid_normal_direction(fluid_block);
    std::cout << "[dbg] fluid_normal_direction ok\n";

    // Register the primitive and gradient buffers before constructing the
    // MUSCL integrators, which acquire these variables by name.
    BaseParticles &fluid_particles_for_vars = fluid_block.getBaseParticles();
    fluid_particles_for_vars.registerStateVariableData<Vecd>("Velocity");
    fluid_particles_for_vars.registerStateVariableData<Vecd>("Momentum");
    fluid_particles_for_vars.registerStateVariableData<Real>("Pressure");
    fluid_particles_for_vars.registerStateVariableData<Real>("TotalEnergy");
    fluid_particles_for_vars.registerStateVariableData<Real>("Density");
    fluid_particles_for_vars.registerStateVariableData<Real>("Mass");
    std::cout << "[dbg] pre-registered compressible vars ok\n";

    SimpleDynamics<ShockTubeInitialCondition> initial_condition(fluid_block);
    std::cout << "[dbg] initial_condition ok\n";

    InteractionWithUpdate<FreeSurfaceIndicationComplex> surface_indicator(fluid_inner, fluid_wall_contact);
    std::cout << "[dbg] surface_indicator ok\n";

    // Keep inner and wall kernel gradients on the same raw stencil so a
    // uniform state remains balanced at the reflective boundaries.
    InteractionWithUpdate<fluid_dynamics::DensityGradient<Inner<NoKernelCorrection>>>
        density_gradient(fluid_inner);
    std::cout << "[dbg] density_gradient ok\n";
    InteractionWithUpdate<fluid_dynamics::VelocityGradient<Inner<NoKernelCorrection>>>
        velocity_gradient(fluid_inner);
    std::cout << "[dbg] velocity_gradient ok\n";
    InteractionWithUpdate<fluid_dynamics::PressureGradient<Inner<NoKernelCorrection>>>
        pressure_gradient(fluid_inner);
    std::cout << "[dbg] pressure_gradient ok\n";

    // First-order reconstruction: keep the MUSCL-WithWall interaction and its
    // conservative HLLC flux, but disable slope reconstruction in this case.
    fluid_dynamics::MUSCLHLLCBridgeConfig muscl_bridge_cfg;
    // Use standard HLLC at both inner and wall interfaces; the optional
    // limiter is not used for this zero-slope shock-tube discretization.
    muscl_bridge_cfg.use_hllc_dissipation_limiter = false;
    muscl_bridge_cfg.gamma = heat_capacity_ratio;
    fluid_dynamics::MUSCLHLLCBridgeConfig wall_bridge_cfg;
    wall_bridge_cfg.use_hllc_dissipation_limiter = false;
    wall_bridge_cfg.gamma = heat_capacity_ratio;
    InteractionWithUpdate<fluid_dynamics::EulerianCompressibleIntegration1stHalfMUSCLWithWall>
        momentum_relaxation(DynamicsArgs(fluid_inner, muscl_bridge_cfg),
                            DynamicsArgs(fluid_wall_contact, wall_bridge_cfg));
    InteractionWithUpdate<fluid_dynamics::EulerianCompressibleIntegration2ndHalfMUSCLWithWall>
        density_and_energy_relaxation(DynamicsArgs(fluid_inner, muscl_bridge_cfg),
                                      DynamicsArgs(fluid_wall_contact, wall_bridge_cfg));

    ReduceDynamics<fluid_dynamics::EulerianCompressibleAcousticTimeStepSize>
        get_fluid_time_step_size(fluid_block, 0.1);
    //----------------------------------------------------------------------
    //	Prepare cell linked list, configuration and initial condition.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    periodic_condition_y.update_cell_linked_list_.exec();
    periodic_condition_wall_left_y.update_cell_linked_list_.exec();
    periodic_condition_wall_right_y.update_cell_linked_list_.exec();
    sph_system.initializeSystemConfigurations();
    wall_left_normal.exec();
    wall_right_normal.exec();
    fluid_normal_direction.exec();
    surface_indicator.exec();
    initial_condition.exec();

    // The shock-tube reference solution is one-dimensional. Keep the
    // conservative MUSCL-HLLC bridge at first order to damp the lattice
    // checkerboard mode in smooth regions.
    std::fill_n(fluid_particles_for_vars.getVariableDataByName<Vecd>("DensityGradient"),
                fluid_particles_for_vars.TotalRealParticles(), Vecd::Zero());
    std::fill_n(fluid_particles_for_vars.getVariableDataByName<Matd>("VelocityGradient"),
                fluid_particles_for_vars.TotalRealParticles(), Matd::Zero());
    std::fill_n(fluid_particles_for_vars.getVariableDataByName<Vecd>("PressureGradient"),
                fluid_particles_for_vars.TotalRealParticles(), Vecd::Zero());
    //----------------------------------------------------------------------
    //	Define I/O: VTP states + CSV probe.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    write_real_body_states.addToWrite<int>(fluid_block, "Indicator");
    write_real_body_states.addToWrite<Real>(fluid_block, "Density");
    write_real_body_states.addToWrite<Real>(fluid_block, "Pressure");
    write_real_body_states.addToWrite<Vecd>(fluid_block, "Velocity");
    ReducedQuantityRecording<MaximumSpeed> write_maximum_speed(fluid_block);

    BaseParticles &fluid_particles = fluid_block.getBaseParticles();
    ShockTubeCsvProbe csv_probe(fluid_particles);
    std::cout << "[Debug] fluid particles = " << fluid_particles.TotalRealParticles()
              << " (first-order HLLC reconstruction)\n";
    //----------------------------------------------------------------------
    //	Setup for time-stepping control.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real output_interval = 0.02; /**< VTP/CSV output cadence. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time.
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile(0);
    csv_probe.write(0, physical_time);
    //----------------------------------------------------------------------
    // Main loop uses fixed zero slopes for first-order interface states.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval && physical_time < end_time)
        {
            Real dt = get_fluid_time_step_size.exec();

            // First-order HLLC reconstruction: do not refresh MUSCL slopes
            // between half steps, so the zero slopes remain invariant.
            momentum_relaxation.exec(dt);

            // The second half uses the same piecewise-constant interface states.
            density_and_energy_relaxation.exec(dt);

            integration_time += dt;
            physical_time += dt;
            if (number_of_iterations % screen_output_interval == 0)
            {
                write_maximum_speed.writeToFile(number_of_iterations);
                std::cout << std::fixed << std::setprecision(9)
                          << "N=" << number_of_iterations << "  t=" << physical_time
                          << "  dt=" << dt << "\n";
            }
            number_of_iterations++;
        }

        TickCount t2 = TickCount::now();
        write_real_body_states.writeToFile();
        csv_probe.write(number_of_iterations, physical_time);

        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    csv_probe.write(number_of_iterations, physical_time);
    return 0;
}
