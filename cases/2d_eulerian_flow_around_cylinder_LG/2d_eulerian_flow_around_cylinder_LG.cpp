/**
 * @file 	2d_eulerian_flow_around_cylinder_LG.cpp
 * @brief 	使用默认 kernel 的弱可压黏性圆柱绕流测试算例。
 * @details 本算例考虑二维 Eulerian 流动经过圆柱的过程。
 * @details 本文件只保留算例主流程编排：读取配置、构造物体/关系、选择 relaxation/restart/remap 路径并推进时间循环。
 * @author 	Zhentong Wang and Xiangyu Hu
 */
#include "cylinder_lg_calculation.hpp"

#include <exception>
#include <iomanip>
#include <memory>

using namespace SPH;
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
/**
 * @brief 2D Eulerian cylinder LG 算例入口，负责配置加载、SPH 系统装配、初始化和主时间推进。
 */
int main(int ac, char *av[])
{
    try
    {
    cylinder_lg::SimulationConfig simulation_config = cylinder_lg::loadCaseConfig();
    if (simulation_config.enable_staged_refinement)
    {
        return cylinder_lg::runStagedRefinementWorkflow(simulation_config, ac, av);
    }
    applyCaseConfig(simulation_config);
    std::cout << "Loaded case config: " << simulation_config.config_file << std::endl;
    const bool configured_restart_mode =
        simulation_config.enable_restart && simulation_config.restart_step != 0;
    const bool remapped_continuation_mode =
        simulation_config.enable_local_refinement && simulation_config.load_remapped_state;
    const bool preserve_output_mode = configured_restart_mode || remapped_continuation_mode;

    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vec2d(-DL_sponge, -DH_sponge), Vec2d(DL, DH + DH_sponge));
    cylinder_lg::preserveOutputFolderForRestart(preserve_output_mode);
    SPHSystem sph_system(system_domain_bounds, global_resolution);
    cylinder_lg::restoreOutputFolderForRestart(preserve_output_mode);
    // Tag for run particle relaxation for the initial body fitted distribution.
    sph_system.setRunParticleRelaxation(simulation_config.run_particle_relaxation);
    // Tag for computation start with relaxed body fitted particles distribution.
    sph_system.setReloadParticles(configured_restart_mode ? false : simulation_config.reload_particles);
    // Handle command line arguments and override the tags for particle relaxation and reload.
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBlock"));
    MultiPolygonShape local_refinement_region = createLocalRefinementRegion();
    if (enable_local_refinement)
    {
        std::cout << "[LocalRefinement] Enabled adaptive particle spacing: level="
                  << local_refinement_level
                  << ", spacing_factor=" << local_refinement_spacing_factor
                  << ", region=[" << local_refinement_region_x_min << ", "
                  << local_refinement_region_x_max << "] x ["
                  << local_refinement_region_y_min << ", "
                  << local_refinement_region_y_max << "]" << std::endl;
        water_block.defineAdaptation<AdaptiveWithinShapeBySpacingFactor>(
            1.3, 1.0, local_refinement_level, local_refinement_spacing_factor);
    }
    water_block.defineComponentLevelSetShape("OuterBoundary");
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    if (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
    {
        water_block.generateParticles<BaseParticles, Reload>(water_block.getName());
        if (load_remapped_state)
        {
            water_block.getBaseParticles()
                .reloadExtraVariable<Real>("Density")
                .reloadExtraVariable<Real>("Pressure")
                .reloadExtraVariable<Vecd>("Velocity");
        }
    }
    else if (enable_local_refinement)
    {
        water_block.generateParticles<BaseParticles, DeterministicLattice>(local_refinement_region);
    }
    else
    {
        water_block.generateParticles<BaseParticles, Lattice>();
    }

    SolidBody cylinder(sph_system, makeShared<Cylinder>("Cylinder"));
    cylinder.defineAdaptationRatios(1.3, 2.0);
    cylinder.defineBodyLevelSetShape();
    cylinder.defineMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? cylinder.generateParticles<BaseParticles, Reload>(cylinder.getName())
        : cylinder.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //	Note that the same relation should be defined only once.
    //----------------------------------------------------------------------
    std::unique_ptr<InnerRelation> water_block_inner_single;
    std::unique_ptr<AdaptiveInnerRelation> water_block_inner_adaptive;
    if (enable_local_refinement)
    {
        water_block_inner_adaptive = std::make_unique<AdaptiveInnerRelation>(water_block);
    }
    else
    {
        water_block_inner_single = std::make_unique<InnerRelation>(water_block);
    }
    BaseInnerRelation &water_block_inner =
        enable_local_refinement ? static_cast<BaseInnerRelation &>(*water_block_inner_adaptive)
                                : static_cast<BaseInnerRelation &>(*water_block_inner_single);
    InnerRelation cylinder_inner(cylinder);
    std::unique_ptr<ContactRelation> water_block_contact_single;
    std::unique_ptr<AdaptiveContactRelation> water_block_contact_adaptive;
    if (enable_local_refinement)
    {
        water_block_contact_adaptive = std::make_unique<AdaptiveContactRelation>(water_block, RealBodyVector{&cylinder});
    }
    else
    {
        water_block_contact_single = std::make_unique<ContactRelation>(water_block, RealBodyVector{&cylinder});
    }
    BaseContactRelation &water_block_contact =
        enable_local_refinement ? static_cast<BaseContactRelation &>(*water_block_contact_adaptive)
                                : static_cast<BaseContactRelation &>(*water_block_contact_single);
    std::unique_ptr<ContactRelation> cylinder_contact_single;
    std::unique_ptr<AdaptiveContactRelation> cylinder_contact_adaptive;
    if (enable_local_refinement)
    {
        cylinder_contact_adaptive = std::make_unique<AdaptiveContactRelation>(cylinder, RealBodyVector{&water_block});
    }
    else
    {
        cylinder_contact_single = std::make_unique<ContactRelation>(cylinder, RealBodyVector{&water_block});
    }
    BaseContactRelation &cylinder_contact =
        enable_local_refinement ? static_cast<BaseContactRelation &>(*cylinder_contact_adaptive)
                                : static_cast<BaseContactRelation &>(*cylinder_contact_single);
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_wall_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<DeterministicRandomizeParticlePosition> random_inserted_body_particles(cylinder);
        SimpleDynamics<DeterministicRandomizeParticlePosition> random_water_body_particles(water_block);
        BodyStatesRecordingToVtp write_real_body_states(sph_system);
        ReloadParticleIO write_real_body_particle_reload_files({&cylinder, &water_block});
        RelaxationStepLevelSetCorrectionInner relaxation_step_inner(cylinder_inner);
        RelaxationStepLevelSetCorrectionComplex relaxation_step_complex(
            DynamicsArgs(water_block_inner, std::string("OuterBoundary")), water_block_contact);
        std::unique_ptr<SimpleDynamics<UpdateSmoothingLengthRatioByShape>> update_smoothing_length_ratio;
        if (enable_local_refinement)
        {
            update_smoothing_length_ratio =
                std::make_unique<SimpleDynamics<UpdateSmoothingLengthRatioByShape>>(water_block, local_refinement_region);
        }
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        random_water_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        relaxation_step_complex.SurfaceBounding().exec();
        write_real_body_states.writeToFile(0);

        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            if (enable_local_refinement)
            {
                update_smoothing_length_ratio->exec();
            }
            relaxation_step_complex.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
                write_real_body_states.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process finish !" << std::endl;

        write_real_body_particle_reload_files.writeToFile(0);

        return 0;
    }
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    InteractionWithUpdate<FreeSurfaceIndicationComplex> surface_indicator(water_block_inner, water_block_contact);
    InteractionDynamics<SmearedSurfaceIndication> smeared_surface(water_block_inner);
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> cylinder_kernel_correction_matrix(cylinder_inner, cylinder_contact);
    std::unique_ptr<InteractionWithUpdate<LinearGradientCorrectionMatrixComplex>> water_block_kernel_correction_matrix_single;
    std::unique_ptr<InteractionWithUpdate<LinearGradientCorrectionMatrixComplex>> water_block_kernel_correction_matrix_adaptive;
    if (enable_local_refinement)
    {
        water_block_kernel_correction_matrix_adaptive =
            std::make_unique<InteractionWithUpdate<LinearGradientCorrectionMatrixComplex>>(
                DynamicsArgs(water_block_inner, Real(0.1)), water_block_contact);
    }
    else
    {
        water_block_kernel_correction_matrix_single =
            std::make_unique<InteractionWithUpdate<LinearGradientCorrectionMatrixComplex>>(
                water_block_inner, water_block_contact);
    }
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> &water_block_kernel_correction_matrix =
        enable_local_refinement ? *water_block_kernel_correction_matrix_adaptive
                                : *water_block_kernel_correction_matrix_single;
    InteractionDynamics<KernelGradientCorrectionComplex> kernel_gradient_update(water_block_inner, water_block_contact);
    SimpleDynamics<NormalDirectionFromBodyShape> cylinder_normal_direction(cylinder);

    Dynamics1Level<EulerianPressureRelaxationWithWallNormalOnly> pressure_relaxation(
        water_block_inner, water_block_contact, simulation_config.pressure_contact_only_normal,
        simulation_config.riemann_limiter_parameter);
    InteractionWithUpdate<fluid_dynamics::EulerianIntegration2ndHalfWithWallRiemann> density_relaxation(
        water_block_inner, water_block_contact);

    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_block_contact);
    InteractionDynamics<FluidWallViscousForceRecorder> fluid_wall_viscous_force_recorder(water_block_contact);
    SimpleDynamics<FluidWallForceReset> reset_fluid_wall_force(water_block);
    SimpleDynamics<NormalDirectionFromBodyShape> water_block_normal_direction(water_block);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block, acoustic_cfl);
    InteractionWithUpdate<FarFieldBoundary> variable_reset_in_boundary_condition(water_block_inner);
    //----------------------------------------------------------------------
    //	Compute the force exerted on solid body due to fluid pressure and viscosity
    //----------------------------------------------------------------------
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(cylinder_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(cylinder_contact);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    write_real_body_states.addToWrite<int>(water_block, "Indicator");
    if (enable_local_refinement)
    {
        write_real_body_states.addToWrite<Real>(water_block, "VolumetricMeasure");
        write_real_body_states.addToWrite<Real>(water_block, "SmoothingLengthRatio");
        write_real_body_states.addToWrite<Real>(water_block, "Density");
    }
    write_real_body_states.addToWrite<Vecd>(water_block, "Velocity");
    write_real_body_states.addToWrite<Real>(water_block, "Pressure");
    BaseParticles &water_block_particles = water_block.getBaseParticles();
    water_block_particles.addEvolvingVariable<Real>("Density");
    water_block_particles.addEvolvingVariable<Real>("Mass");
    water_block_particles.addEvolvingVariable<Real>("Pressure");
    water_block_particles.addEvolvingVariable<Real>("DensityChangeRate");
    water_block_particles.addEvolvingVariable<Real>("MassChangeRate");
    water_block_particles.addEvolvingVariable<Vecd>("Velocity");
    water_block_particles.addEvolvingVariable<Vecd>("Momentum");
    water_block_particles.addEvolvingVariable<Vecd>("MomentumChangeRate");
    water_block_particles.addEvolvingVariable<Vecd>("Force");
    water_block_particles.addEvolvingVariable<Vecd>("ForcePrior");
    if (enable_local_refinement && !load_remapped_state)
    {
        initializeEulerianStateFromDensityAndVelocity(water_block);
    }
    else if (enable_local_refinement)
    {
        initializeEulerianStateFromRemappedReload(water_block);
    }
    // Regression test is disabled by default for this case on deployment environments.
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    if (enable_local_refinement && !reportRelationFiniteState(water_block, water_block_inner, water_block_contact, "after initial configuration"))
    {
        throw std::runtime_error("Local refinement detected non-finite relation after initial configuration.");
    }
    cylinder_normal_direction.exec();
    surface_indicator.exec();
    smeared_surface.exec();
    water_block_normal_direction.exec();
    variable_reset_in_boundary_condition.exec();
    if (enable_local_refinement && !reportRelationFiniteState(water_block, water_block_inner, water_block_contact, "after far field initialization"))
    {
        throw std::runtime_error("Local refinement detected non-finite relation after far field initialization.");
    }
    cylinder_kernel_correction_matrix.exec();
    water_block_kernel_correction_matrix.exec();
    if (enable_local_refinement && !reportRelationFiniteState(water_block, water_block_inner, water_block_contact, "after correction matrix"))
    {
        throw std::runtime_error("Local refinement detected non-finite relation after correction matrix.");
    }
    kernel_gradient_update.exec();
    if (enable_local_refinement)
    {
        if (!reportRelationFiniteState(water_block, water_block_inner, water_block_contact, "after configuration"))
        {
            throw std::runtime_error("Local refinement detected non-finite relation after configuration.");
        }
    }
    LocalRefinementConservationStats local_refinement_initial_conservation;
    if (enable_local_refinement)
    {
        local_refinement_initial_conservation =
            computeLocalRefinementConservationStats(water_block);
        if (local_refinement_runtime_diagnostics)
        {
            std::cout << "[LocalRefinement] Runtime diagnostics enabled." << std::endl;
        }
    }
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    if (enable_local_refinement && load_remapped_state)
    {
        physical_time = remap_physical_time;
        std::cout << "[LocalRefinement][Remap] Loaded remapped state from reload at t="
                  << physical_time << "s." << std::endl;
    }
    const int screen_output_interval = simulation_config.screen_output_interval;
    const int restart_output_interval =
        screen_output_interval * simulation_config.restart_output_factor;
    const Real end_time = simulation_config.end_time;
    const Real output_interval = simulation_config.output_interval; /**< time stamps for output. */
    bool enable_restart = simulation_config.enable_restart;
    int restart_step_config = simulation_config.restart_step;
    size_t restart_step = 0;
    if (enable_restart && restart_step_config == -1)
    {
        const int detected_step = cylinder_lg::detectLatestRestartStep();
        if (detected_step > 0)
        {
            restart_step = static_cast<size_t>(detected_step);
        }
        else
        {
            throw std::runtime_error("Cannot continue restart: no valid restart file found.");
        }
    }
    else if (enable_restart && restart_step_config > 0)
    {
        restart_step = static_cast<size_t>(restart_step_config);
    }

    sph_system.setRestartStep(restart_step);
    std::unique_ptr<RestartIO> restart_io;
    if (enable_restart)
    {
        restart_io = std::make_unique<RestartIO>(sph_system);
    }
    const bool is_restart = enable_restart && restart_step > 0;
    const bool is_output_continuation = is_restart || remapped_continuation_mode;
    if (is_restart)
    {
        std::cout << "Restart mode: reading from step " << restart_step << std::endl;
        cylinder_lg::requireRestartAttributes(
            restart_step,
            {"Density", "Mass", "Pressure", "Velocity", "Momentum",
             "MomentumChangeRate", "MassChangeRate"});
        physical_time = restart_io->readRestartFiles(restart_step);
        water_block.updateCellLinkedList();
        cylinder.updateCellLinkedList();
        water_wall_complex.updateConfiguration();
        cylinder_inner.updateConfiguration();
        cylinder_contact.updateConfiguration();
        cylinder_normal_direction.exec();
        surface_indicator.exec();
        smeared_surface.exec();
        water_block_normal_direction.exec();
        cylinder_kernel_correction_matrix.exec();
        water_block_kernel_correction_matrix.exec();
        kernel_gradient_update.exec();

        const std::string output_folder = sph_system.getIOEnvironment().OutputFolder();
        for (const auto &filename : cylinder_lg::restartDatFilenamesToTruncate(output_folder))
        {
            cylinder_lg::truncateDatFileToTime(output_folder + "/" + filename, physical_time);
        }
        cylinder_lg::cleanupFutureVtpFiles(output_folder, physical_time);
    }
    else if (remapped_continuation_mode)
    {
        const std::string output_folder = sph_system.getIOEnvironment().OutputFolder();
        for (const auto &filename : cylinder_lg::restartDatFilenamesToTruncate(output_folder))
        {
            cylinder_lg::truncateDatFileToTime(output_folder + "/" + filename, physical_time);
        }
        cylinder_lg::cleanupFutureVtpFiles(output_folder, physical_time);
        std::cout << "[RESTART] Remapped continuation output will continue from t = "
                  << physical_time << std::endl;
    }
    else
    {
        std::cout << (enable_restart
                          ? "Fresh start: checkpoint output is enabled; no restart file is loaded."
                          : "Fresh start: restart I/O is disabled.")
                  << std::endl;
    }
    size_t number_of_iterations = is_restart ? restart_step : 0;
    const int restart_keep_last_n = simulation_config.restart_keep_last_n;
    const bool write_final_restart = enable_restart;
    size_t next_restart_output_iteration =
        restart_step + static_cast<size_t>(restart_output_interval);
    if (enable_restart)
    {
        std::cout << "Restart output interval: every " << restart_output_interval
                  << " iteration(s)." << std::endl;
        std::cout << "restart_keep_last_n = " << restart_keep_last_n
                  << (restart_keep_last_n < 0 ? " (keep all)" : "") << std::endl;
        if (write_final_restart)
        {
            std::cout << "Final restart output is enabled." << std::endl;
        }
    }
    size_t last_restart_written_iteration = restart_step;
    cylinder_lg::RestartSafeReducedQuantityRecording<QuantitySummation<Vecd>>
        write_total_viscous_force_from_fluid(is_output_continuation, cylinder, "ViscousForceFromFluid");
    cylinder_lg::RestartSafeReducedQuantityRecording<QuantitySummation<Vecd>>
        write_total_pressure_force_from_fluid_body(is_output_continuation, cylinder, "PressureForceFromFluid");
    FluidWallForceRecording write_fluid_wall_force(water_block, is_output_continuation);
    cylinder_lg::RestartSafeReducedQuantityRecording<MaximumSpeed> write_maximum_speed(
        is_output_continuation, water_block);
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    if (!is_output_continuation)
    {
        write_real_body_states.writeToFile(0);
    }
    else
    {
        std::cout << "[RESTART] Initial output is skipped at t = " << physical_time
                  << " to avoid duplicate frames." << std::endl;
    }
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        reset_fluid_wall_force.exec();
        const bool run_local_refinement_runtime_diagnostics =
            enable_local_refinement && local_refinement_runtime_diagnostics;
        auto report_local_refinement_conservation =
            [&](const std::string &stage, bool force_report = false)
        {
            if (enable_local_refinement && force_report)
            {
                reportLocalRefinementConservationStats(
                    water_block, stage, local_refinement_initial_conservation);
            }
        };
        auto run_acoustic_step = [&](Real dt, bool run_viscous_stage)
        {
            if (enable_local_refinement)
            {
                if (!isFiniteReal(dt) || dt <= Real(0.0))
                {
                    reportFluidFiniteState(water_block, "before invalid dt");
                    throw std::runtime_error("Local refinement detected invalid dt before dynamics.");
                }
            }
            if (run_local_refinement_runtime_diagnostics &&
                !reportRelationFiniteState(water_block, water_block_inner, water_block_contact, "before pressure chain"))
            {
                throw std::runtime_error("Local refinement detected non-finite relation before pressure chain.");
            }
            if (run_viscous_stage)
            {
                viscous_force.exec();
                size_t first_bad_particle = water_block.getBaseParticles().TotalRealParticles();
                if (run_local_refinement_runtime_diagnostics &&
                    !reportFluidFiniteState(water_block, "after viscous_force", &first_bad_particle))
                {
                    reportViscousForceNeighborhoodDiagnostics(water_block, water_block_inner, water_block_contact, first_bad_particle, "after viscous_force");
                    throw std::runtime_error("Local refinement detected non-finite state after viscous_force.");
                }
                fluid_wall_viscous_force_recorder.exec();
                if (run_local_refinement_runtime_diagnostics &&
                    !reportFluidFiniteState(water_block, "after fluid_wall_viscous_force_recorder"))
                {
                    throw std::runtime_error("Local refinement detected non-finite state after fluid wall force recorder.");
                }
            }
            pressure_relaxation.exec(dt);
            if (run_local_refinement_runtime_diagnostics &&
                !reportFluidFiniteState(water_block, "after pressure_relaxation"))
            {
                throw std::runtime_error("Local refinement detected non-finite state after pressure_relaxation.");
            }
            report_local_refinement_conservation("after pressure_relaxation");
            density_relaxation.exec(dt);
            if (run_local_refinement_runtime_diagnostics &&
                !reportFluidFiniteState(water_block, "after density_relaxation"))
            {
                throw std::runtime_error("Local refinement detected non-finite state after density_relaxation.");
            }
            report_local_refinement_conservation("after density_relaxation");
            size_t density_non_positive_particle = water_block.getBaseParticles().TotalRealParticles();
            if (run_local_refinement_runtime_diagnostics &&
                !reportFluidPositiveState(
                    water_block, "after density_relaxation",
                    &density_non_positive_particle, true))
            {
                report_local_refinement_conservation("after density_relaxation failure", true);
                reportDensityRelaxationMassFluxDiagnostic(
                    water_block, water_block_inner, water_block_contact,
                    density_non_positive_particle, dt, "after density_relaxation");
                throw std::runtime_error("Local refinement detected non-positive density or mass after density_relaxation.");
            }
            if (run_local_refinement_runtime_diagnostics)
            {
                reportFarFieldBoundaryRiskScan(water_block, water_block_inner, "before far_field_boundary");
            }

            integration_time += dt;
            physical_time += dt;
            variable_reset_in_boundary_condition.exec();
            if (run_local_refinement_runtime_diagnostics &&
                !reportFluidFiniteState(water_block, "after far_field_boundary"))
            {
                throw std::runtime_error("Local refinement detected non-finite state after far field boundary.");
            }
            report_local_refinement_conservation("after far_field_boundary");
            size_t first_non_positive_particle = water_block.getBaseParticles().TotalRealParticles();
            if (run_local_refinement_runtime_diagnostics &&
                !reportFluidPositiveState(
                    water_block, "after far_field_boundary",
                    &first_non_positive_particle))
            {
                report_local_refinement_conservation("after far_field_boundary failure", true);
                writeFarFieldBoundaryDiagnostic(
                    water_block, water_block_inner, first_non_positive_particle,
                    "after far_field_boundary");
                throw std::runtime_error("Local refinement detected non-positive density or mass after far field boundary.");
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	dt = " << dt << "\n";
            }
            number_of_iterations++;
        };

        while (integration_time < output_interval)
        {
            Real dt = get_fluid_time_step_size.exec();
            run_acoustic_step(dt, true);
        }
        if (enable_local_refinement && !local_refinement_runtime_diagnostics)
        {
            // 正式运行只在输出窗口末尾做一次重检查，避免每个 acoustic step 全量扫描拖慢并行计算。
            if (!reportFluidFiniteState(water_block, "after output interval"))
            {
                throw std::runtime_error("Local refinement detected non-finite state after output interval.");
            }
            size_t first_non_positive_particle = water_block.getBaseParticles().TotalRealParticles();
            if (!reportFluidPositiveState(
                    water_block, "after output interval", &first_non_positive_particle))
            {
                report_local_refinement_conservation("after output interval failure", true);
                throw std::runtime_error("Local refinement detected non-positive density or mass after output interval.");
            }
        }

        TickCount t2 = TickCount::now();
        write_real_body_states.writeToFile();
        if (enable_restart &&
            number_of_iterations >= next_restart_output_iteration &&
            number_of_iterations != restart_step)
        {
            restart_io->writeToFile(number_of_iterations);
            last_restart_written_iteration = number_of_iterations;
            cylinder_lg::cleanupOldCheckpoints("restart", restart_keep_last_n);
            while (next_restart_output_iteration <= number_of_iterations)
            {
                next_restart_output_iteration += static_cast<size_t>(restart_output_interval);
            }
        }
        viscous_force_from_fluid.exec();
        pressure_force_from_fluid.exec();
        write_total_viscous_force_from_fluid.writeToFile(number_of_iterations);
        write_total_pressure_force_from_fluid_body.writeToFile(number_of_iterations);
        write_fluid_wall_force.writeToFile(number_of_iterations, integration_time);

        write_maximum_speed.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    if (enable_restart &&
        write_final_restart &&
        number_of_iterations != restart_step &&
        number_of_iterations != last_restart_written_iteration)
    {
        restart_io->writeToFile(number_of_iterations);
        cylinder_lg::cleanupOldCheckpoints("restart", restart_keep_last_n);
    }

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "[FATAL] " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    catch (...)
    {
        std::cerr << "[FATAL] Unknown exception." << std::endl;
        return EXIT_FAILURE;
    }
}
