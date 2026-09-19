/**
 * @file 	multiphase_shock_tube_asr.cpp
 * @brief 	2D Eulerian SPH multiphase shock tube with particle-band ASR.
 *          First-order five-equation Godunov scheme, symmetric kernel
 *          gradient (Eq. 3), local-h acoustic time step, reflective walls
 *          at the x ends and y periodicity.
 *
 *          Configs select static graded lattices, split/merge events, band
 *          tracking and the optional Shepard filter.
 * @author 	KIYOYOZU
 */
#include "multiphase_shock_tube_asr.h"

#include "kernel_hyperbolic.h"
#include "particle_reserve.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
using namespace SPH;
//----------------------------------------------------------------------
//	CSV probe: particle dump with adaptation variables for post-processing.
//----------------------------------------------------------------------
class AsrCsvProbe
{
  public:
    AsrCsvProbe(BaseParticles &particles, const std::string &tag,
                const std::string &case_dir)
        : particles_(particles), tag_(tag),
          // CSV dumps live next to their config file (per-case folder)
          output_dir_((std::filesystem::path(case_dir) / "output").string()),
          pos_(particles.getVariableDataByName<Vecd>("Position")),
          rho_(particles.getVariableDataByName<Real>("Density")),
          p_(particles.getVariableDataByName<Real>("Pressure")),
          vel_(particles.getVariableDataByName<Vecd>("Velocity")),
          alpha_(particles.getVariableDataByName<Real>("VolumeFraction")),
          Vol_(particles.getVariableDataByName<Real>("VolumetricMeasure")),
          band_(particles.getVariableDataByName<int>("ParticleBand")),
          h_ratio_(particles.getVariableDataByName<Real>("SmoothingLengthRatio")),
          ref_spacing_(particles.getVariableDataByName<Real>("ReferenceSpacing")),
          rho_raw_(particles.getVariableDataByName<Real>("DensityRaw")),
          p_raw_(particles.getVariableDataByName<Real>("PressureRaw"))
    {
        std::filesystem::create_directories(output_dir_);
    }

    void write(size_t step, Real physical_time) const
    {
        std::string prefix = tag_.empty() ? "" : tag_ + "_";
        std::ofstream ofs(output_dir_ + "/" + prefix + "particles_" + std::to_string(step) + ".csv");
        ofs << "x,y,rho,p,u,v,alpha,Vol,band,h_ratio,ref_spacing,rho_raw,p_raw\n";
        ofs << std::setprecision(12);
        for (size_t i = 0; i != particles_.TotalRealParticles(); ++i)
        {
            ofs << pos_[i][0] << "," << pos_[i][1] << ","
                << rho_[i] << "," << p_[i] << ","
                << vel_[i][0] << "," << vel_[i][1] << "," << alpha_[i] << ","
                << Vol_[i] << "," << band_[i] << ","
                << h_ratio_[i] << "," << ref_spacing_[i] << ","
                << rho_raw_[i] << "," << p_raw_[i] << "\n";
        }
        std::ofstream meta(output_dir_ + "/" + prefix + "time_" + std::to_string(step) + ".txt");
        meta << physical_time << "\n";
        std::cout << "[Probe] wrote " << prefix << "particles_" << step
                  << ".csv (N=" << particles_.TotalRealParticles()
                  << ", t=" << physical_time << ")\n";
    }

  protected:
    BaseParticles &particles_;
    std::string tag_, output_dir_;
    Vecd *pos_, *vel_;
    Real *rho_, *p_, *alpha_, *Vol_, *h_ratio_, *ref_spacing_, *rho_raw_, *p_raw_;
    int *band_;
};
//----------------------------------------------------------------------
//	Conservation diagnostics: mass, momentum, energy and mixture volume.
//----------------------------------------------------------------------
class AsrConservationProbe
{
  public:
    AsrConservationProbe(BaseParticles &particles)
        : particles_(particles),
          mass_(particles.getVariableDataByName<Real>("Mass")),
          mom_(particles.getVariableDataByName<Vecd>("Momentum")),
          E_(particles.getVariableDataByName<Real>("TotalEnergy")),
          Vol_(particles.getVariableDataByName<Real>("VolumetricMeasure")),
          alpha_(particles.getVariableDataByName<Real>("VolumeFraction")) {};

    void write(Real physical_time, const std::string &label = "") const
    {
        Real total_mass = 0.0, total_mom_x = 0.0, total_energy = 0.0;
        Real total_vol = 0.0, total_alpha_vol = 0.0;
        for (size_t i = 0; i != particles_.TotalRealParticles(); ++i)
        {
            total_mass += mass_[i];
            total_mom_x += mom_[i][0];
            total_energy += E_[i];
            total_vol += Vol_[i];
            total_alpha_vol += alpha_[i] * Vol_[i];
        }
        std::cout << "[Conservation]" << label << " t=" << physical_time
                  << "  N=" << particles_.TotalRealParticles()
                  << std::setprecision(15)
                  << "  mass=" << total_mass
                  << "  mom_x=" << total_mom_x
                  << "  energy=" << total_energy
                  << "  vol=" << total_vol
                  << "  alpha_vol=" << total_alpha_vol << "\n";
    }

  protected:
    BaseParticles &particles_;
    Real *mass_, *E_, *Vol_, *alpha_;
    Vecd *mom_;
};
//----------------------------------------------------------------------
//	Wall impulse accumulator (acceptance probe). Records the impulse the fluid
//	delivers to each wall, i.e. minus the wall-part of the momentum
//	update that the integrator applies (the probe replays the identical
//	reflective-ghost HLLC wall flux of the 1st-half wall pass, term by
//	term, from the same pre-half state the integrator is about to
//	consume). Walls are static here (no FSI coupling), so the mirrored
//	wall velocity is exactly -vel_i; adding FSI would require reading
//	the wall AverageVelocity like the integrator does. The correction's
//	wall gauge term is NOT included: with the correction enabled the
//	integrator + correction wall force differs by the gauge, so probe
//	values are the raw-scheme impulse (shock-tube runs use correction off).
//----------------------------------------------------------------------
class AsrWallImpulseProbe
{
  public:
    AsrWallImpulseProbe(BaseParticles &particles, BaseContactRelation &wall_contact,
                        MultiphaseMixture &mixture)
        : particles_(particles), wall_contact_(wall_contact), solver_(mixture),
          Vol_(particles.getVariableDataByName<Real>("VolumetricMeasure")),
          rho_(particles.getVariableDataByName<Real>("Density")),
          p_(particles.getVariableDataByName<Real>("Pressure")),
          E_(particles.getVariableDataByName<Real>("TotalEnergy")),
          alpha_(particles.getVariableDataByName<Real>("VolumeFraction")),
          vel_(particles.getVariableDataByName<Vecd>("Velocity")),
          impulse_(wall_contact.getContactParticles().size(), Vecd::Zero()) {};

    void accumulate(Real dt)
    {
        StdVec<BaseParticles *> contact_particles = wall_contact_.getContactParticles();
        StdVec<ParticleConfiguration> &contact_configuration = wall_contact_.contact_configuration_;
        StdVec<Real *> wall_Vol;
        for (size_t k = 0; k != contact_particles.size(); ++k)
            wall_Vol.push_back(contact_particles[k]->getVariableDataByName<Real>("VolumetricMeasure"));

        for (size_t index_i = 0; index_i != particles_.TotalRealParticles(); ++index_i)
        {
            Real rho_i = rho_[index_i], p_i = p_[index_i];
            Real E_i = E_[index_i] / Vol_[index_i], alpha_i = alpha_[index_i];
            Vecd vel_i = vel_[index_i];
            MultiphaseFluidState state_i(rho_i, vel_i, p_i, E_i, alpha_i);
            Vecd vel_reflect = -vel_i;
            MultiphaseFluidState state_g(rho_i, vel_reflect, p_i, E_i, alpha_i);

            for (size_t k = 0; k != contact_configuration.size(); ++k)
            {
                Neighborhood &contact_neighborhood = contact_configuration[k][index_i];
                for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
                {
                    Real dW_ijV_j = contact_neighborhood.dW_ij_[n] *
                                    wall_Vol[k][contact_neighborhood.j_[n]];
                    Vecd e_ij = contact_neighborhood.e_ij_[n];
                    MultiphaseFluidStarState interface_state =
                        solver_.getInterfaceState(state_i, state_g, e_ij);
                    Matd convect_flux = interface_state.rho_ * interface_state.vel_ *
                                        interface_state.vel_.transpose();
                    Vecd force_on_fluid = -2.0 * Vol_[index_i] * dW_ijV_j *
                                          (convect_flux + interface_state.p_ * Matd::Identity()) * e_ij;
                    impulse_[k] -= force_on_fluid * dt;
                }
            }
        }
    }

    void write(Real physical_time) const
    {
        std::cout << "[WallImpulse] t=" << physical_time << std::setprecision(12);
        for (size_t k = 0; k != impulse_.size(); ++k)
            std::cout << "  wall[" << k << "]=(" << impulse_[k][0]
                      << "," << impulse_[k][1] << ")";
        std::cout << "\n";
    }

  protected:
    BaseParticles &particles_;
    BaseContactRelation &wall_contact_;
    MultiphaseHLLCRiemannSolver solver_;
    Real *Vol_, *rho_, *p_, *E_, *alpha_;
    Vecd *vel_;
    StdVec<Vecd> impulse_;
};
//----------------------------------------------------------------------
//	First-moment gate of the cached symmetric gradient,
//	M_i = |sum_j Vol_j dW_ij e_ij|. A uniform quiescent state stays at rest
//	only up to this residual, so it must be on the level of the uniform
//	reference lattice.
//----------------------------------------------------------------------
class AsrFirstMomentProbe
{
  public:
    AsrFirstMomentProbe(BaseParticles &particles, ParticleConfiguration &inner_configuration)
        : particles_(particles), inner_configuration_(inner_configuration),
          Vol_(particles.getVariableDataByName<Real>("VolumetricMeasure")),
          grad_corr_(particles.getVariableDataByName<Vecd>("GradientCorrection")) {};

    void write(const std::string &label) const
    {
        Real max_moment = 0.0, sum_sq = 0.0, max_corrected = 0.0;
        size_t count = particles_.TotalRealParticles();
        for (size_t i = 0; i != count; ++i)
        {
            Vecd moment = Vecd::Zero();
            Real vol_sum = 0.0;
            Neighborhood &neighborhood = inner_configuration_[i];
            for (size_t n = 0; n != neighborhood.current_size_; ++n)
            {
                moment += Vol_[neighborhood.j_[n]] * neighborhood.dW_ij_[n] * neighborhood.e_ij_[n];
                vol_sum += Vol_[neighborhood.j_[n]];
            }
            Real norm = moment.norm();
            max_moment = SMAX(max_moment, norm);
            sum_sq += norm * norm;
            // inner-only residual after correction: -> 0 away from the wall;
            // wall-adjacent particles keep the wall part of c_i as a residual
            max_corrected = SMAX(max_corrected, (moment - grad_corr_[i] * vol_sum).norm());
        }
        std::cout << "[FirstMoment]" << label
                  << "  max=" << std::setprecision(12) << max_moment
                  << "  rms=" << std::sqrt(sum_sq / Real(count))
                  << "  corrected_max_inner=" << max_corrected << "\n";
    }

  protected:
    BaseParticles &particles_;
    ParticleConfiguration &inner_configuration_;
    Real *Vol_;
    Vecd *grad_corr_;
};
//----------------------------------------------------------------------
//	Pair-list audit of the inner configuration.
//	1) reciprocity: every cached pair (i,j) appears in j's list with the
//	   same dW (machine precision) -- precondition for exact conservation;
//	2) completeness: brute-force enumeration with the y minimum-image
//	   distance finds no missing or extra pair relative to the support
//	   criterion r < rc_ref / min(h_ratio_i, h_ratio_j).
//----------------------------------------------------------------------
class AsrSymmetryAudit
{
  public:
    AsrSymmetryAudit(BaseParticles &particles, ParticleConfiguration &inner_configuration,
                     Real cutoff_ref, Real H)
        : particles_(particles), inner_configuration_(inner_configuration),
          cutoff_ref_(cutoff_ref), H_(H),
          pos_(particles.getVariableDataByName<Vecd>("Position")),
          h_ratio_(particles.getVariableDataByName<Real>("SmoothingLengthRatio")) {};

    bool run() const
    {
        size_t count = particles_.TotalRealParticles();
        size_t asymmetric_pairs = 0, dW_mismatch = 0, ambiguous_asymmetric = 0;
        size_t missing_pairs = 0, extra_pairs = 0, ambiguous_pairs = 0;

        // dW comparison tolerance: the two search paths can round the pair
        // distance differently at the 1e-14 relative level; through the kernel
        // Lipschitz bound this allows an absolute cached-dW difference of
        // ~1e-13 * dw_max, which matters only near the cutoff where dW itself
        // vanishes quadratically and amplifies the noise to a large relative
        // (but physically negligible) difference
        Real dw_max = 0.0;
        for (size_t i = 0; i != count; ++i)
        {
            Neighborhood &ni = inner_configuration_[i];
            for (size_t n = 0; n != ni.current_size_; ++n)
                dw_max = SMAX(dw_max, std::abs(ni.dW_ij_[n]));
        }

        // reciprocity, image-aware: every entry (i -> j at distance r) must
        // have a counterpart entry (j -> i at the same distance) with an
        // identical cached dW (precondition for pairwise flux antisymmetry)
        for (size_t i = 0; i != count; ++i)
        {
            Neighborhood &ni = inner_configuration_[i];
            for (size_t n = 0; n != ni.current_size_; ++n)
            {
                size_t j = ni.j_[n];
                Neighborhood &nj = inner_configuration_[j];
                bool found = false;
                for (size_t m = 0; m != nj.current_size_; ++m)
                {
                    if (nj.j_[m] != i)
                        continue;
                    Real r_scale = SMAX(ni.r_ij_[n], TinyReal);
                    if (std::abs(nj.r_ij_[m] - ni.r_ij_[n]) <= 1e-12 * r_scale)
                    {
                        Real dw_tol = 1e-13 * dw_max + 1e-12 * std::abs(ni.dW_ij_[n]);
                        if (std::abs(nj.dW_ij_[m] - ni.dW_ij_[n]) <= dw_tol)
                        {
                            found = true;
                            break;
                        }
                        if (dW_mismatch < 5)
                            std::cout << "  [audit] dW mismatch: i=" << i << " j=" << j
                                      << " r=" << ni.r_ij_[n] << "/" << nj.r_ij_[m]
                                      << " dW=" << ni.dW_ij_[n] << "/" << nj.dW_ij_[m]
                                      << " hr=" << h_ratio_[i] << "/" << h_ratio_[j]
                                      << " pos_i=(" << pos_[i][0] << "," << pos_[i][1] << ")"
                                      << " pos_j=(" << pos_[j][0] << "," << pos_[j][1] << ")\n";
                        dW_mismatch++;
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    // pairs at the cutoff within FP noise are exempt: the
                    // mirrored displacement arithmetic can round the distance
                    // to either side of the cutoff on the two search paths,
                    // and the cached dW vanishes at the cutoff anyway
                    Real pair_cutoff = cutoff_ref_ / SMIN(h_ratio_[i], h_ratio_[j]);
                    if (std::abs(ni.r_ij_[n] - pair_cutoff) < 1e-9 * pair_cutoff)
                    {
                        ambiguous_asymmetric++;
                        continue;
                    }
                    asymmetric_pairs++;
                    if (asymmetric_pairs <= 5)
                        std::cout << "  [audit] no counterpart: i=" << i << " j=" << j
                                  << " r=" << ni.r_ij_[n]
                                  << " pos_i=(" << pos_[i][0] << "," << pos_[i][1] << ")"
                                  << " pos_j=(" << pos_[j][0] << "," << pos_[j][1] << ")"
                                  << " hr_i=" << h_ratio_[i] << " hr_j=" << h_ratio_[j] << "\n";
                }
            }
        }

        // completeness: enumerate all periodic images (k = -1, 0, +1) of every
        // pair and check presence/absence against the support criterion
        // r < rc_ref / min(h_ratio_i, h_ratio_j)
        for (size_t i = 0; i != count; ++i)
            for (size_t j = i + 1; j != count; ++j)
            {
                Real pair_cutoff = cutoff_ref_ / SMIN(h_ratio_[i], h_ratio_[j]);
                for (int k = -1; k <= 1; ++k)
                {
                    Vecd disp = pos_[i] - pos_[j] - Vecd(0.0, Real(k) * H_);
                    Real distance = disp.norm();
                    bool expected = distance < pair_cutoff;
                    // images closer than 1e-9 relative to the cutoff are not judged
                    if (std::abs(distance - pair_cutoff) < 1e-9 * pair_cutoff)
                    {
                        ambiguous_pairs++;
                        continue;
                    }
                    bool found_i = listContainsImage(i, j, distance);
                    bool found_j = listContainsImage(j, i, distance);
                    if (expected && (!found_i || !found_j))
                    {
                        missing_pairs++;
                        if (missing_pairs <= 5)
                            std::cout << "  [audit] missing: i=" << i << " j=" << j
                                      << " k=" << k << " r=" << distance
                                      << " cutoff=" << pair_cutoff
                                      << " in_i=" << found_i << " in_j=" << found_j << "\n";
                    }
                    if (!expected && (found_i || found_j))
                    {
                        extra_pairs++;
                        if (extra_pairs <= 5)
                            std::cout << "  [audit] extra: i=" << i << " j=" << j
                                      << " k=" << k << " r=" << distance
                                      << " cutoff=" << pair_cutoff << "\n";
                    }
                }
            }

        std::cout << "[SymmetryAudit] N=" << count
                  << "  asymmetric_pairs=" << asymmetric_pairs
                  << "  dW_mismatch=" << dW_mismatch
                  << "  missing_pairs=" << missing_pairs
                  << "  extra_pairs=" << extra_pairs
                  << "  ambiguous(cutoff-band)=" << ambiguous_pairs
                  << "  ambiguous_asymmetric=" << ambiguous_asymmetric << "\n";
        bool pass = asymmetric_pairs == 0 && dW_mismatch == 0 &&
                    missing_pairs == 0 && extra_pairs == 0;
        std::cout << "[SymmetryAudit] " << (pass ? "PASS" : "FAIL") << "\n";
        return pass;
    }

  protected:
    BaseParticles &particles_;
    ParticleConfiguration &inner_configuration_;
    Real cutoff_ref_, H_;
    Vecd *pos_;
    Real *h_ratio_;

    bool listContainsImage(size_t i, size_t j, Real distance) const
    {
        const Neighborhood &ni = inner_configuration_[i];
        for (size_t n = 0; n != ni.current_size_; ++n)
            if (ni.j_[n] == j &&
                std::abs(ni.r_ij_[n] - distance) <= 1e-9 * SMAX(distance, TinyReal))
                return true;
        return false;
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    // config file selection: --config=path/config.ini (default: config.ini
    // in the working directory; CSV dumps go to output/ next to the config)
    std::string config_file = "config.ini";
    for (int i = 1; i < ac; ++i)
    {
        std::string arg(av[i]);
        if (arg.rfind("--config=", 0) == 0)
            config_file = arg.substr(9);
    }
    std::string case_dir = std::filesystem::absolute(config_file).parent_path().string();
    AsrConfig cfg = loadAsrConfig(config_file);
    std::cout << "[ASR] config = " << config_file
              << "  preset = " << cfg.preset
              << "  banding = " << (cfg.initial_banding ? "graded" : "uniform")
              << "  kernel = " << cfg.kernel << "\n";

    const Real dp = cfg.dp;
    const Real L = cfg.L;
    const Real H = cfg.H;
    const Real wall_thickness = cfg.wall_thickness;
    const Real BW = 4.0 * dp;
    BoundingBoxd system_domain_bounds(Vec2d(-wall_thickness - BW, -BW),
                                      Vec2d(L + wall_thickness + BW, H + BW));
    //----------------------------------------------------------------------
    //	System, bodies and the particle-band adaptation.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, dp);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(false);

    FluidBody fluid_block(sph_system, makeShared<AsrFluidBlock>("FluidBlock", L, H));
    // the band structure follows ds_max_factor in graded AND uniform mode
    // (USR baselines set ds_max_factor = 1 explicitly; split/merge event tests need
    // the full band set on a uniform lattice to drive split/merge)
    fluid_block.defineAdaptation<ParticleBandAdaptation>(
        cfg.ds_max_factor, cfg.band_coef, cfg.band_width_factor, cfg.h_spacing_ratio);
    ParticleBandAdaptation &adaptation =
        DynamicCast<ParticleBandAdaptation>(&fluid_block, fluid_block.getSPHAdaptation());
    if (cfg.kernel == "hyperbolic")
    adaptation.resetKernel<KernelHyperbolic>();
    fluid_block.defineComponentLevelSetShape("OuterBoundary");
    fluid_block.defineMatterMaterial<StiffenedGas>(cfg.gamma1, cfg.p_inf1);

    ParticleBuffer<ReserveSizeFactor> particle_buffer(ReserveSizeFactor(cfg.buffer_growth));
    fluid_block.generateParticlesWithReserve<BaseParticles, BandedLattice>(particle_buffer, cfg);

    // Walls are generated AFTER the fluid: their particles are exact mirror
    // images of the fluid lattice across each wall face (see MirrorWallLattice).
    BaseParticles &fluid_particles_for_walls = fluid_block.getBaseParticles();
    SolidBody wall_left(sph_system, makeShared<AsrWallBlock>("WallLeft", -wall_thickness, 0.0, H));
    // The wall periodic bounding band equals the wall kernel cutoff; it must
    // reach at least the fluid cutoff, else coarse wall rows adjacent to the
    // y seam miss their periodic images, M_wall gains a transverse component
    // and any (p*_wall - p_i) mismatch jets fluid along y at the corner rows
    wall_left.defineAdaptationRatios(1.6 * cfg.ds_max_factor, 1.0);
    wall_left.defineBodyLevelSetShape();
    wall_left.defineMatterMaterial<Solid>();
    wall_left.generateParticles<BaseParticles, MirrorWallLattice>(
        fluid_particles_for_walls, 0.0, 1.0, wall_thickness);

    SolidBody wall_right(sph_system, makeShared<AsrWallBlock>("WallRight", L, L + wall_thickness, H));
    wall_right.defineAdaptationRatios(1.6 * cfg.ds_max_factor, 1.0);
    wall_right.defineBodyLevelSetShape();
    wall_right.defineMatterMaterial<Solid>();
    wall_right.generateParticles<BaseParticles, MirrorWallLattice>(
        fluid_particles_for_walls, L, -1.0, wall_thickness);
    //----------------------------------------------------------------------
    //	Two-phase mixture and body relations. The inner relation caches the
    //	symmetric kernel gradient of Eq. (3); the wall contact uses the
    //	adaptive builder so that wall pairs see the fluid's local h
    //	(h_ratio_min = fluid ratio, wall ratio mapped through relative_h_ref
    //	is always larger), keeping the wall stencil the mirror image of the
    //	inner stencil at h_i.
    //----------------------------------------------------------------------
    StiffenedGas gas_material(cfg.gamma1, cfg.p_inf1);
    StiffenedGas water_material(cfg.gamma2, cfg.p_inf2);
    MultiphaseMixture mixture(gas_material, water_material);

    ParticleBandInnerRelation fluid_inner(fluid_block);
    AdaptiveContactRelation fluid_wall_contact(fluid_block, RealBodyVector{&wall_left, &wall_right});

    BoundingBoxd periodic_bounds(Vec2d(0.0, 0.0), Vec2d(L, H));
    PeriodicAlongAxis periodic_along_y(periodic_bounds, yAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition_y(fluid_block, periodic_along_y);
    PeriodicConditionUsingCellLinkedList periodic_condition_wall_left_y(wall_left, periodic_along_y);
    PeriodicConditionUsingCellLinkedList periodic_condition_wall_right_y(wall_right, periodic_along_y);
    //----------------------------------------------------------------------
    //	Numerics: wall normals, variable registration, initial condition.
    //----------------------------------------------------------------------
    SimpleDynamics<AsrWallNormal> wall_left_normal(wall_left, Vecd(1.0, 0.0));
    SimpleDynamics<AsrWallNormal> wall_right_normal(wall_right, Vecd(-1.0, 0.0));

    BaseParticles &fluid_particles = fluid_block.getBaseParticles();
    fluid_particles.registerStateVariableData<Real>("Density");
    fluid_particles.registerStateVariableData<Real>("Mass");
    fluid_particles.registerStateVariableData<Vecd>("Velocity");
    fluid_particles.registerStateVariableData<Vecd>("Momentum");
    fluid_particles.registerStateVariableData<Real>("Pressure");
    fluid_particles.registerStateVariableData<Real>("TotalEnergy");
    fluid_particles.registerStateVariableData<Real>("VolumeFraction");
    // Spawn/Remove copy the evolving set only; the conserved fluid state must
    // survive split/merge events (rho/p/vel are recovered from it afterwards)
    fluid_particles.addEvolvingVariable<Real>("Mass");
    fluid_particles.addEvolvingVariable<Vecd>("Momentum");
    fluid_particles.addEvolvingVariable<Real>("TotalEnergy");
    fluid_particles.addEvolvingVariable<Real>("VolumeFraction");

    SimpleDynamics<ShockTubeAsrInitialCondition> initial_condition(fluid_block, mixture, cfg);
    //----------------------------------------------------------------------
    //	Build cell linked lists, periodic ghosts and configurations.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    periodic_condition_y.update_cell_linked_list_.exec();
    periodic_condition_wall_left_y.update_cell_linked_list_.exec();
    periodic_condition_wall_right_y.update_cell_linked_list_.exec();
    sph_system.initializeSystemConfigurations();
    wall_left_normal.exec();
    wall_right_normal.exec();
    initial_condition.exec();
    // test hook: force every particle into one band so that the split
    // (gamma > gamma_s) or merge (gamma < gamma_m) criterion fires on the
    // whole uniform lattice deterministically
    if (cfg.force_band >= 0)
    {
        int *band = fluid_particles.getVariableDataByName<int>("ParticleBand");
        Real *ref_spacing = fluid_particles.getVariableDataByName<Real>("ReferenceSpacing");
        for (size_t i = 0; i != fluid_particles.TotalRealParticles(); ++i)
        {
            band[i] = cfg.force_band;
            ref_spacing[i] = adaptation.BandSpacing(cfg.force_band);
        }
        std::cout << "[ASR] force_band = " << cfg.force_band
                  << " applied (ds_band = " << adaptation.BandSpacing(cfg.force_band) << ")\n";
    }
    //----------------------------------------------------------------------
    //	Time step (local h) and first-order five-equation integrators.
    //----------------------------------------------------------------------
    ReduceDynamics<fluid_dynamics::EulerianMultiphaseAcousticTimeStepSizeLocalH>
        get_fluid_time_step_size(fluid_block, mixture, cfg.acoustic_cfl);

    InteractionWithUpdate<fluid_dynamics::EulerianMultiphaseIntegration1stHalfWithWall>
        momentum_relaxation(DynamicsArgs(fluid_inner, mixture),
                            DynamicsArgs(fluid_wall_contact, mixture));
    InteractionWithUpdate<fluid_dynamics::EulerianMultiphaseIntegration2ndHalfWithWall>
        density_energy_alpha_relaxation(DynamicsArgs(fluid_inner, mixture),
                                        DynamicsArgs(fluid_wall_contact, mixture));

    // riemann_order == 2: MUSCL-reconstructed interface states, gradients on
    // the raw (uncorrected) kernel stencil over the inner neighborhood only
    // (validated recipe of the sibling uniform case)
    fluid_dynamics::SecondOrderConfig soc;
    if (cfg.limiter == "mc")
        soc.limiter = fluid_dynamics::SlopeLimiter::MC;
    else if (cfg.limiter == "vanleer")
        soc.limiter = fluid_dynamics::SlopeLimiter::VanLeer;
    else if (cfg.limiter == "none")
        soc.limiter = fluid_dynamics::SlopeLimiter::None;
    else
        soc.limiter = fluid_dynamics::SlopeLimiter::Minmod;
    soc.piecewise_rho_alpha = (cfg.reconstruct == "vel_p");
    InteractionWithUpdate<fluid_dynamics::DensityGradient<Inner<NoKernelCorrection>>>
        density_gradient(fluid_inner);
    InteractionWithUpdate<fluid_dynamics::VelocityGradient<Inner<NoKernelCorrection>>>
        velocity_gradient(fluid_inner);
    InteractionWithUpdate<fluid_dynamics::PressureGradient<Inner<NoKernelCorrection>>>
        pressure_gradient(fluid_inner);
    InteractionWithUpdate<fluid_dynamics::VolumeFractionGradient>
        alpha_gradient(fluid_inner);
    InteractionWithUpdate<fluid_dynamics::EulerianMultiphaseIntegration1stHalfMUSCLWithWall>
        momentum_relaxation_muscl(DynamicsArgs(fluid_inner, mixture, soc),
                                  DynamicsArgs(fluid_wall_contact, mixture, soc));
    InteractionWithUpdate<fluid_dynamics::EulerianMultiphaseIntegration2ndHalfMUSCLWithWall>
        density_energy_alpha_relaxation_muscl(DynamicsArgs(fluid_inner, mixture, soc),
                                              DynamicsArgs(fluid_wall_contact, mixture, soc));
    const bool second_order = (cfg.riemann_order == 2);

    // well-balanced consistent-flux correction: cancels the spurious
    // first-moment x frozen-flux residual at band boundaries for any locally
    // uniform state (identity on a uniform lattice, see module doc)
    MultiphaseBackgroundPressureCorrection steady_correction(fluid_inner, fluid_wall_contact,
                                                             mixture);

    UpdateSmoothingLengthByBand update_smoothing_length(fluid_block, fluid_inner);
    // zeroth-order gradient correction: c_i closes the kernel first moment over
    // the inner + wall stencil, must be (re)computed after every configuration
    // rebuild and before the fluxes consume it
    ComputeGradientCorrection compute_gradient_correction(fluid_inner, &fluid_wall_contact);
    compute_gradient_correction.exec(); // initial c_i from the t=0 configuration
    UpdateParticleBands update_particle_bands(fluid_block, fluid_inner,
                                              cfg.interface_alpha_tol, cfg.band_hysteresis,
                                              cfg.band_tracking == "dijkstra",
                                              cfg.shock_band, cfg.shock_rel_jump);
    ParticleSplittingByBand particle_splitting(fluid_block, fluid_inner,
                                               cfg.gamma_split, cfg.split_lambda, H);
    ParticleMergingByBand particle_merging(fluid_block, fluid_inner, mixture,
                                           cfg.gamma_merge);
    ShepardDensityFilter shepard_filter(fluid_block, fluid_inner, mixture);
    //----------------------------------------------------------------------
    //	I/O and diagnostics.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    write_real_body_states.addToWrite<Real>(fluid_block, "Density");
    write_real_body_states.addToWrite<Real>(fluid_block, "Pressure");
    write_real_body_states.addToWrite<Vecd>(fluid_block, "Velocity");
    write_real_body_states.addToWrite<Real>(fluid_block, "VolumeFraction");
    write_real_body_states.addToWrite<int>(fluid_block, "ParticleBand");
    write_real_body_states.addToWrite<Real>(fluid_block, "SmoothingLengthRatio");

    // tag outputs so that graded/uniform runs do not overwrite
    std::string tag = std::filesystem::path(config_file).stem().string();
    if (!cfg.initial_banding)
        tag += "_uniform";
    AsrCsvProbe csv_probe(fluid_particles, tag, case_dir);
    AsrConservationProbe conservation_probe(fluid_particles);
    AsrWallImpulseProbe wall_impulse_probe(fluid_particles, fluid_wall_contact, mixture);
    AsrFirstMomentProbe moment_probe(fluid_particles, fluid_inner.inner_configuration_);
    AsrSymmetryAudit symmetry_audit(fluid_particles, fluid_inner.inner_configuration_,
                                    adaptation.getKernel()->CutOffRadius(), H);

    std::cout << "[ASR] fluid particles = " << fluid_particles.TotalRealParticles()
              << "  bands = " << adaptation.BandCount() + 1
              << "  h_ref = " << adaptation.ReferenceSmoothingLength()
              << "  N_r = " << adaptation.ReferenceNeighborNumber() << "\n";
    std::cout << "[ASR] riemann_order = " << cfg.riemann_order
              << (cfg.riemann_order == 2 ? " (MUSCL: limiter=" + cfg.limiter +
                                           ", reconstruct=" + cfg.reconstruct + ")"
                                         : " (first-order)")
              << "\n";
    if (cfg.check_symmetry && !symmetry_audit.run())
    {
        std::cerr << "[FATAL] pair-list symmetry audit failed; aborting." << std::endl;
        return 1;
    }
    if (cfg.check_symmetry)
    {
        // band-interaction restriction (paper: neighbors span at most ±1 band)
        int *band = fluid_particles.getVariableDataByName<int>("ParticleBand");
        size_t cross_pairs = 0;
        size_t total = fluid_particles.TotalRealParticles();
        for (size_t i = 0; i != total; ++i)
        {
            Neighborhood &nb = fluid_inner.inner_configuration_[i];
            for (size_t n = 0; n != nb.current_size_; ++n)
                if (std::abs(band[i] - band[nb.j_[n]]) > 1)
                    cross_pairs++;
        }
        std::cout << "[BandCheck] cross-band(>1) neighbor pairs = " << cross_pairs << "\n";
        if (cross_pairs != 0)
        {
            std::cerr << "[FATAL] band-interaction restriction violated." << std::endl;
            return 1;
        }
    }
    if (cfg.print_moment)
        moment_probe.write(" t=0");

    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    Real output_interval = cfg.output_interval;
    // event_window mode: count steps since the last split/merge event
    size_t steps_since_event = std::numeric_limits<size_t>::max() / 4;

    TickCount t1 = TickCount::now();
    TimeInterval interval;

    write_real_body_states.writeToFile(0);
    csv_probe.write(0, physical_time);
    conservation_probe.write(physical_time);
    //----------------------------------------------------------------------
    //	Main loop. Adaptation events (when enabled) run every adapt_interval
    //	steps and are always followed by a full cell-linked-list, periodic
    //	ghost and configuration rebuild before any neighbor data is consumed.
    //----------------------------------------------------------------------
    const bool adaptation_enabled = cfg.adapt_h || cfg.adapt_split_merge || cfg.adapt_bands;

    while (physical_time < cfg.end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval && physical_time < cfg.end_time)
        {
            if (adaptation_enabled && number_of_iterations % cfg.adapt_interval == 0)
            {
                // band tracking first: it sets the spacing targets that the
                // smoothing-length update and the split/merge criteria consume
                size_t band_changes = 0;
                if (cfg.adapt_bands)
                    band_changes = update_particle_bands.exec();
                if (cfg.adapt_h)
                    update_smoothing_length.exec();
                size_t splits = 0, merges = 0;
                if (cfg.adapt_split_merge)
                {
                    conservation_probe.write(physical_time, " pre-event");
                    splits = particle_splitting.exec();
                    merges = particle_merging.exec();
                    conservation_probe.write(physical_time, " post-event");
                }
                if (band_changes > 0 && splits == 0 && merges == 0)
                    std::cout << "[BandTrack] iter=" << number_of_iterations
                              << "  band_changes=" << band_changes
                              << "  x_if=" << update_particle_bands.DetectedInterface()
                              << "  x_shock=" << update_particle_bands.DetectedShock()
                              << "  N=" << fluid_particles.TotalRealParticles() << "\n";
                if (splits > 0 || merges > 0)
                {
                    steps_since_event = 0;
                    std::cout << "[AdaptEvent] iter=" << number_of_iterations
                              << "  band_changes=" << band_changes
                              << "  x_if=" << update_particle_bands.DetectedInterface()
                              << "  x_shock=" << update_particle_bands.DetectedShock()
                              << "  splits=" << splits << "  merges=" << merges
                              << "  N=" << fluid_particles.TotalRealParticles() << "\n";
                }
                sph_system.initializeSystemCellLinkedLists();
                periodic_condition_y.update_cell_linked_list_.exec();
                periodic_condition_wall_left_y.update_cell_linked_list_.exec();
                periodic_condition_wall_right_y.update_cell_linked_list_.exec();
                sph_system.initializeSystemConfigurations();
                // recompute c_i against the rebuilt dW/positions before any
                // flux consumes the new configuration
                compute_gradient_correction.exec();
                if (cfg.check_symmetry && (splits > 0 || merges > 0))
                {
                    // no stale indices: after the compression the rebuilt pair
                    // lists must stay symmetric and complete, and the band
                    // restriction must still hold
                    if (!symmetry_audit.run())
                    {
                        std::cerr << "[FATAL] post-event pair-list audit failed." << std::endl;
                        return 1;
                    }
                }
            }

            Real dt = get_fluid_time_step_size.exec();
            if (!std::isfinite(dt) || dt <= TinyReal)
            {
                std::cerr << "[FATAL] time step collapsed (dt=" << dt
                          << ") at t=" << physical_time
                          << " -- configuration unstable; stopping.\n";
                physical_time = cfg.end_time;
                break;
            }

            // accumulate the wall impulse from the same pre-half state the
            // integrator is about to consume (see probe doc)
            wall_impulse_probe.accumulate(dt);
            if (second_order)
            {
                // in vel_p mode rho/alpha stay piecewise constant, so their
                // gradients are not needed
                if (!soc.piecewise_rho_alpha)
                {
                    density_gradient.exec();
                    alpha_gradient.exec();
                }
                velocity_gradient.exec();
                pressure_gradient.exec();
            }
            // the frozen fluxes must be built from the state the halves
            // consume, not from the mid-step state polluted by the raw
            // band-boundary kick -- hence the snapshot before the 1st half
            if (cfg.steady_correction)
                steady_correction.snapshotState();
            if (second_order)
                momentum_relaxation_muscl.exec(dt);
            else
                momentum_relaxation.exec(dt);
            // the frozen momentum flux must be removed before the 2nd half:
            // its mass/E/alpha fluxes inherit the uncorrected band-boundary
            // pseudo-velocity, which feeds the alpha upwind term at material
            // interfaces and corrupts the EOS pressure of the next step
            if (cfg.steady_correction)
                steady_correction.execMomentumBetweenHalves(
                    dt, cfg.corr_upwind_select, cfg.corr_interface_density_ratio);
            if (second_order)
                density_energy_alpha_relaxation_muscl.exec(dt);
            else
                density_energy_alpha_relaxation.exec(dt);
            // the frozen mass/energy divergence is injected by the 2nd half
            // itself, so its removal (with the full state recovery) comes after
            if (cfg.steady_correction)
                steady_correction.execMassEnergyAfterHalves(
                    dt, cfg.corr_upwind_select, cfg.corr_interface_density_ratio);

            // Shepard override after the 2nd-half state recovery, so the next
            // step's fluxes see the filtered density/pressure
            bool shepard_active =
                cfg.shepard_filter == "on" ||
                (cfg.shepard_filter == "event_window" &&
                 steps_since_event < size_t(cfg.shepard_window_steps));
            if (shepard_active)
                shepard_filter.exec();
            if (steps_since_event < std::numeric_limits<size_t>::max() / 4)
                steps_since_event++;

            integration_time += dt;
            physical_time += dt;
            if (number_of_iterations % cfg.screen_output_interval == 0)
            {
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
        wall_impulse_probe.write(physical_time);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }

    TickCount t4 = TickCount::now();
    TimeInterval tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    csv_probe.write(number_of_iterations, physical_time);
    conservation_probe.write(physical_time);
    wall_impulse_probe.write(physical_time);
    if (cfg.print_moment)
        moment_probe.write(" t=end");

    // auto post-process: plot the latest profile next to output/ (the script
    // sits in the case family root, i.e. the parent of cases/)
    {
        // the postprocess script lives in the case-family root (the directory
        // that contains cases/), two levels above the case directory
        std::filesystem::path script = std::filesystem::path(case_dir).parent_path()
                                           .parent_path() / "postprocess_asr.py";
        if (std::filesystem::exists(script))
        {
            // pass the case directory itself so cases living outside the
            // default cases/ tree (e.g. quiescent gates) are also handled
            std::string cmd = "python \"" + script.string() + "\" \"" + case_dir + "\"";
            std::cout << "[ASR] postprocess: " << cmd << std::endl;
            int rc = std::system(cmd.c_str());
            if (rc != 0)
                std::cerr << "[ASR] WARNING: postprocess exited with code " << rc << std::endl;
        }
        else
        {
            std::cout << "[ASR] postprocess_asr.py not found next to cases/, skipping." << std::endl;
        }
    }
    return 0;
}
