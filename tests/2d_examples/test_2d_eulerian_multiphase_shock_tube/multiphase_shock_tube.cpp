/**
 * @file 	multiphase_shock_tube.cpp
 * @brief 	2D Eulerian SPH multiphase shock tube (gas-water interface),
 *          Kapila five-equation model with stiffened gas EOS.
 *          x ends are real SolidBody reflective walls (ContactRelation +
 *          five-equation WithWall integrators); y direction is periodic.
 *          First-order (piecewise constant) Godunov scheme: interface states
 *          come from the five-equation HLLC Riemann solver with Wood sound
 *          speed, no slope reconstruction.
 * @author 	KIYOYOZU
 */
#include "multiphase_shock_tube.h"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace SPH;
//----------------------------------------------------------------------
//	CSV probe: dump all real particles (x, y, rho, p, u, alpha) for
//	post-processing. Written to <case>/output regardless of the cwd.
//----------------------------------------------------------------------
class MultiPhaseShockTubeCsvProbe
{
  public:
    MultiPhaseShockTubeCsvProbe(BaseParticles &particles)
        : particles_(particles),
          // CSV dumps go to output/ next to the resolved config.ini, so each
          // case folder (cases/<name>/) keeps its own data when running from it
          output_dir_((std::filesystem::path(mp_cfg_detail::resolvePath())
                           .parent_path() /
                       "output")
                          .string()),
          pos_(particles.getVariableDataByName<Vecd>("Position")),
          rho_(particles.getVariableDataByName<Real>("Density")),
          p_(particles.getVariableDataByName<Real>("Pressure")),
          vel_(particles.getVariableDataByName<Vecd>("Velocity")),
          alpha_(particles.getVariableDataByName<Real>("VolumeFraction"))
    {
        std::filesystem::create_directories(output_dir_);
    }

    void write(size_t step, Real physical_time) const
    {
        std::ofstream ofs(output_dir_ + "/particles_" + std::to_string(step) + ".csv");
        ofs << "x,y,rho,p,u,alpha\n";
        ofs << std::setprecision(12);
        for (size_t i = 0; i != particles_.TotalRealParticles(); ++i)
        {
            ofs << pos_[i][0] << "," << pos_[i][1] << ","
                << rho_[i] << "," << p_[i] << ","
                << vel_[i][0] << "," << alpha_[i] << "\n";
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
    Real *alpha_;
};
//----------------------------------------------------------------------
//	Conservation diagnostics: total mass, x-momentum and total energy.
//----------------------------------------------------------------------
class ConservationProbe
{
  public:
    ConservationProbe(BaseParticles &particles)
        : particles_(particles),
          mass_(particles.getVariableDataByName<Real>("Mass")),
          mom_(particles.getVariableDataByName<Vecd>("Momentum")),
          E_(particles.getVariableDataByName<Real>("TotalEnergy")) {};

    void write(Real physical_time) const
    {
        Real total_mass = 0.0;
        Real total_mom_x = 0.0;
        Real total_energy = 0.0;
        for (size_t i = 0; i != particles_.TotalRealParticles(); ++i)
        {
            total_mass += mass_[i];
            total_mom_x += mom_[i][0];
            total_energy += E_[i];
        }
        std::cout << "[Conservation] t=" << physical_time
                  << "  mass=" << total_mass
                  << "  mom_x=" << total_mom_x
                  << "  energy=" << total_energy << "\n";
    }

  protected:
    BaseParticles &particles_;
    Real *mass_, *E_;
    Vecd *mom_;
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
    // run from a case folder (cases/<name>/): vtp output lands in its own
    // output/ next to the config.ini that resolvePath() picked up
    std::string out_root = std::filesystem::path(mp_cfg_detail::resolvePath())
                               .parent_path()
                               .generic_string();
    IO::getEnvironment().resetOutputFolder(out_root + "/output");
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody fluid_block(sph_system, makeShared<FluidBlock>("FluidBlock"));
    fluid_block.defineComponentLevelSetShape("OuterBoundary");
    // The body material only satisfies the framework; the physical EOS is
    // carried by the MultiphaseMixture below.
    fluid_block.defineMatterMaterial<StiffenedGas>(gamma_gas, p_inf_gas);
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
    //	Two-phase mixture: phase 1 = gas, phase 2 = water. The mixture holds
    //	references to the two standalone stiffened-gas materials and is shared
    //	by the initial condition, the integrators and the time-step controller.
    //----------------------------------------------------------------------
    StiffenedGas gas_material(gamma_gas, p_inf_gas);
    StiffenedGas water_material(gamma_water, p_inf_water);
    MultiphaseMixture mixture(gas_material, water_material);
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
    //	Define the numerical methods. Construction order matters: normals ->
    //	variable registration -> initial condition -> indicator -> the two
    //	five-equation half steps.
    //----------------------------------------------------------------------
    SimpleDynamics<WallNormal> wall_left_normal(wall_left, Vecd(1.0, 0.0));
    SimpleDynamics<WallNormal> wall_right_normal(wall_right, Vecd(-1.0, 0.0));
    SimpleDynamics<NormalDirectionFromBodyShape> fluid_normal_direction(fluid_block);

    // Register the conservative and primitive buffers before constructing the
    // initial condition and integrators, which acquire these variables by name.
    // Density and Mass are not created by default for an Eulerian fluid body.
    BaseParticles &fluid_particles_for_vars = fluid_block.getBaseParticles();
    fluid_particles_for_vars.registerStateVariableData<Real>("Density");
    fluid_particles_for_vars.registerStateVariableData<Real>("Mass");
    fluid_particles_for_vars.registerStateVariableData<Vecd>("Velocity");
    fluid_particles_for_vars.registerStateVariableData<Vecd>("Momentum");
    fluid_particles_for_vars.registerStateVariableData<Real>("Pressure");
    fluid_particles_for_vars.registerStateVariableData<Real>("TotalEnergy");
    fluid_particles_for_vars.registerStateVariableData<Real>("VolumeFraction");

    SimpleDynamics<MultiPhaseShockTubeInitialCondition> initial_condition(fluid_block, mixture);

    InteractionWithUpdate<FreeSurfaceIndicationComplex> surface_indicator(fluid_inner, fluid_wall_contact);

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

    // Acoustic time step (Wood sound speed), common to both orders.
    ReduceDynamics<fluid_dynamics::EulerianMultiphaseAcousticTimeStepSize>
        get_fluid_time_step_size(fluid_block, mixture, MP_CFG.acoustic_cfl);
    //----------------------------------------------------------------------
    //	Define I/O: VTP states + CSV probe + conservation diagnostics.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    write_real_body_states.addToWrite<int>(fluid_block, "Indicator");
    write_real_body_states.addToWrite<Real>(fluid_block, "Density");
    write_real_body_states.addToWrite<Real>(fluid_block, "Pressure");
    write_real_body_states.addToWrite<Vecd>(fluid_block, "Velocity");
    write_real_body_states.addToWrite<Real>(fluid_block, "VolumeFraction");
    ReducedQuantityRecording<MaximumSpeed> write_maximum_speed(fluid_block);

    BaseParticles &fluid_particles = fluid_block.getBaseParticles();
    MultiPhaseShockTubeCsvProbe csv_probe(fluid_particles);
    ConservationProbe conservation_probe(fluid_particles);
    std::cout << "[Debug] fluid particles = " << fluid_particles.TotalRealParticles()
              << "  riemann_order = " << MP_CFG.riemann_order << "\n";
    //----------------------------------------------------------------------
    //	Setup for time-stepping control.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real output_interval = MP_CFG.output_interval; /**< VTP/CSV cadence. */
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
    conservation_probe.write(physical_time);
    //----------------------------------------------------------------------
    //	Main loop. riemann_order == 2 uses MUSCL-reconstructed interface states
    //	(gradients recomputed each step); riemann_order == 1 is the piecewise-
    //	constant baseline with reflective walls.
    //----------------------------------------------------------------------
    if (MP_CFG.riemann_order == 2)
    {
        // MUSCL gradients on the raw (uncorrected) kernel stencil: keeps the
        // inner and wall stencils identical so a uniform state stays balanced
        // at the reflective boundaries (validated recipe of
        // test_2d_eulerian_shock_tube_LG). The fluxes also use the raw dW.
        InteractionWithUpdate<fluid_dynamics::DensityGradient<Inner<NoKernelCorrection>>>
            density_gradient(fluid_inner);
        InteractionWithUpdate<fluid_dynamics::VelocityGradient<Inner<NoKernelCorrection>>>
            velocity_gradient(fluid_inner);
        InteractionWithUpdate<fluid_dynamics::PressureGradient<Inner<NoKernelCorrection>>>
            pressure_gradient(fluid_inner);
        InteractionWithUpdate<fluid_dynamics::VolumeFractionGradient>
            alpha_gradient(fluid_inner);

        fluid_dynamics::SecondOrderConfig soc;
        // NOTE: soc.gamma (ideal-gas energy reconstruction inside
        // reconstruct_primitives_muscl) is intentionally NOT set here: the
        // five-equation bridge overwrites that energy with the mixture EOS.
        std::string limiter = MP_CFG.limiter;
        std::transform(limiter.begin(), limiter.end(), limiter.begin(), ::tolower);
        if (limiter == "mc")
            soc.limiter = fluid_dynamics::SlopeLimiter::MC;
        else if (limiter == "vanleer")
            soc.limiter = fluid_dynamics::SlopeLimiter::VanLeer;
        else if (limiter == "none")
            soc.limiter = fluid_dynamics::SlopeLimiter::None;
        else if (limiter == "minmod")
            soc.limiter = fluid_dynamics::SlopeLimiter::Minmod;
        else
        {
            std::cerr << "[MUSCL] WARNING: unknown muscl_limiter '" << MP_CFG.limiter
                      << "', falling back to minmod.\n";
            soc.limiter = fluid_dynamics::SlopeLimiter::Minmod;
        }

        // Reconstruction scope: full = rho/vel/p/alpha (good for moderate
        // impedance ratios, e.g. two-gas Sod); vel_p = vel/p only with rho and
        // alpha piecewise constant (needed for stiff interfaces like gas-water,
        // where reconstructing rho/alpha across the gamma/p_inf jump drives
        // velocity spikes at the contact).
        std::string recon = MP_CFG.reconstruct;
        std::transform(recon.begin(), recon.end(), recon.begin(), ::tolower);
        if (recon != "full" && recon != "vel_p")
        {
            std::cerr << "[MUSCL] WARNING: unknown muscl_reconstruct '" << MP_CFG.reconstruct
                      << "', falling back to full.\n";
            recon = "full";
        }
        soc.piecewise_rho_alpha = (recon == "vel_p");

        // MUSCL with reflective walls (inner + wall contact), same composition
        // as the first-order WithWall integrators.
        InteractionWithUpdate<fluid_dynamics::EulerianMultiphaseIntegration1stHalfMUSCLWithWall>
            momentum_relaxation(DynamicsArgs(fluid_inner, mixture, soc),
                                DynamicsArgs(fluid_wall_contact, mixture, soc));
        InteractionWithUpdate<fluid_dynamics::EulerianMultiphaseIntegration2ndHalfMUSCLWithWall>
            density_energy_alpha_relaxation(DynamicsArgs(fluid_inner, mixture, soc),
                                            DynamicsArgs(fluid_wall_contact, mixture, soc));

        while (physical_time < end_time)
        {
            Real integration_time = 0.0;
            while (integration_time < output_interval && physical_time < end_time)
            {
                Real dt = get_fluid_time_step_size.exec();
                if (!std::isfinite(dt) || dt <= TinyReal)
                {
                    std::cerr << "[FATAL] time step collapsed (dt=" << dt
                              << ") at t=" << physical_time
                              << " -- configuration unstable; stopping.\n";
                    physical_time = end_time;
                    break;
                }

                // In vel_p mode rho/alpha stay piecewise constant, so their
                // gradients are not needed.
                if (!soc.piecewise_rho_alpha)
                {
                    density_gradient.exec();
                    alpha_gradient.exec();
                }
                velocity_gradient.exec();
                pressure_gradient.exec();
                momentum_relaxation.exec(dt);
                density_energy_alpha_relaxation.exec(dt);

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
            conservation_probe.write(physical_time);

            TickCount t3 = TickCount::now();
            interval += t3 - t2;
        }
    }
    else
    {
        // First-order five-equation integration (piecewise constant, HLLC fluxes).
        InteractionWithUpdate<fluid_dynamics::EulerianMultiphaseIntegration1stHalfWithWall>
            momentum_relaxation(DynamicsArgs(fluid_inner, mixture),
                                DynamicsArgs(fluid_wall_contact, mixture));
        InteractionWithUpdate<fluid_dynamics::EulerianMultiphaseIntegration2ndHalfWithWall>
            density_energy_alpha_relaxation(DynamicsArgs(fluid_inner, mixture),
                                            DynamicsArgs(fluid_wall_contact, mixture));

        while (physical_time < end_time)
        {
            Real integration_time = 0.0;
            while (integration_time < output_interval && physical_time < end_time)
            {
                Real dt = get_fluid_time_step_size.exec();
                if (!std::isfinite(dt) || dt <= TinyReal)
                {
                    std::cerr << "[FATAL] time step collapsed (dt=" << dt
                              << ") at t=" << physical_time
                              << " -- configuration unstable; stopping.\n";
                    physical_time = end_time;
                    break;
                }

                momentum_relaxation.exec(dt);
                density_energy_alpha_relaxation.exec(dt);

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
            conservation_probe.write(physical_time);

            TickCount t3 = TickCount::now();
            interval += t3 - t2;
        }
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    csv_probe.write(number_of_iterations, physical_time);
    conservation_probe.write(physical_time);
    return 0;
}
