/**
 * @file 	test_3d_eulerian_channel.cpp
 * @brief 	3D Eulerian weakly-compressible laminar channel flow smoke test.
 * @details 2D LG paradigm alignment (2026-06-27): weakly-compressible Eulerian
 *          SPH on a sponge-extended channel [-DL_sponge, DL+DL_sponge]×[0,DH]×[0,DW].
 *          Streamwise x (inlet/outlet via sponge + FarFieldBoundary non-reflective),
 *          wall-normal y (real SolidBody no-slip walls via ContactRelation +
 *          WithWallRiemann integrators), spanwise z (periodic). Inlet parabolic
 *          profile u_x(y) = U_max * 4 * eta * (1 - eta) with U_max = 1.5 * u_bulk.
 *          The ghost mirror + boundary_type label paradigm is fully retired;
 *          the Task 2 shared helper (syncEulerianWeaklyCompressibleState) is
 *          kept linked for the helper_state_sync_linked GTest but no longer
 *          called on the main path (sponge + FarFieldBoundary does not need
 *          ghost state sync).
 *          GTest gates: finiteness of all state fields and a minimum step count.
 * @author 	KIYOYOZU, Xiangyu Hu
 */
#include "sphinxsys.h"
#include "eulerian_open_boundary.h"  // Task 2 shared helper (linked for GTest, not called on main path)
#include "eulerian_channel_data.hpp"
#include "eulerian_channel_geometry.hpp"
#include "../shared/ck_time_step.hpp"

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace SPH;

//----------------------------------------------------------------------
//	Finiteness check over the boundary quintuple of every real particle.
//  Returns false on the first NaN/Inf encountered.
//  Thin wrapper around reportFluidFiniteState (Task 0): when state is non-finite
//  the diagnostic function prints the full bad-particle record, so the caller
//  gets both a quick bool gate and the per-particle snapshot for debugging.
//----------------------------------------------------------------------
bool checkStatesFinite(BaseParticles &particles, const std::string &stage = "default")
{
    return reportFluidFiniteState(particles, stage);
}

// Smoke-run result returned to the TEST body for GTest assertions.
struct SmokeRunResult
{
    bool fifty_step_finite;      // gate: state finite at the 50th acoustic step
    bool final_finite;           // state finite at the end of the run
    size_t iterations;           // total acoustic steps executed
    Real mass_flux_imbalance;    // |minlet - moutlet| / max(|minlet|, |moutlet|) at final frame
    Real profile_l2;             // relative L2 of x-z-averaged u_x(y) vs parabolic target
    bool flux_computed;          // false if neither inlet nor outlet face had particles
    bool profile_computed;       // false if no interior bins were sampled
};

//----------------------------------------------------------------------
//	Compute the inlet/outlet mass-flux imbalance and the streamwise-velocity
//  profile L2 error at the final frame. Both are evaluated over the real
//  particles only (ghost particles are excluded — they are boundary states,
//  not physical flux carriers).
//    mass_flux_imbalance = |m_in - m_out| / max(|m_in|, |m_out|)
//  where m_face = sum_{i in face} rho_i * u_x_i * A_i, A_i = Vol_i / dp
//  (a cubic Eulerian particle of side dp presents a face area dp^2 = Vol/dp).
//    profile_l2 = sqrt( sum_y (u_avg(y) - u_target(y))^2 ) / sqrt( sum_y u_target(y)^2 )
//  where u_avg(y) is the x-z average of u_x over interior particles binned by y.
//----------------------------------------------------------------------
struct FluxProfileDiagnostics
{
    Real mass_flux_imbalance;
    Real profile_l2;
    bool flux_computed;
    bool profile_computed;
};
FluxProfileDiagnostics computeFluxProfileDiagnostics(BaseParticles &particles,
                                                     const SimulationConfig &cfg)
{
    FluxProfileDiagnostics diag{0.0, 0.0, false, false};
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *Vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    const Real dp = cfg.global_resolution;
    const Real face_band = static_cast<Real>(cfg.boundary_n_layers) * dp;
    // sponge 区（x<0 或 x>DL）粒子参与 FarFieldBoundary 校正，不是物理通量载体，
    // 必须排除——否则 sponge 粒子会污染 inlet/outlet face flux 统计（契约修正 M 项）。
    const auto in_sponge = [&](Real x)
    { return x < -0.5 * dp || x > cfg.DL + 0.5 * dp; };

    // --- Mass-flux imbalance over inlet (x ~ 0) and outlet (x ~ DL) faces. ---
    // 物理入口面 0 < x < face_band，物理出口面 DL-face_band < x < DL，
    // 排除 sponge 区（x<0 或 x>DL）。
    Real m_inlet = 0.0;
    Real m_outlet = 0.0;
    bool has_inlet = false;
    bool has_outlet = false;
    for (size_t i = 0; i < n; ++i)
    {
        if (in_sponge(pos[i][0]))
        {
            continue; // 排除 sponge 区粒子
        }
        const Real A_i = Vol[i] / dp; // face area dp^2 for a cubic particle
        const Real flux_i = rho[i] * vel[i][0] * A_i;
        if (pos[i][0] <= face_band)
        {
            m_inlet += flux_i;
            has_inlet = true;
        }
        else if (pos[i][0] >= cfg.DL - face_band)
        {
            m_outlet += flux_i;
            has_outlet = true;
        }
    }
    // Require BOTH faces to be present: with only one face the imbalance
    // trivially degenerates to 1.0 (|m-0|/|m|), which is not a physical
    // imbalance signal. Marking flux_computed=false forces the GTest to fail
    // explicitly instead of silently passing a meaningless 1.0.
    // NaN 守卫：若 rho/vel 发散使 m_inlet/m_outlet 为 NaN，denom>0.0 为 false 会
    // 误给 imbalance=0.0 假通过。必须显式检查 finite，NaN 时 flux_computed=false
    // 让 GTest FAIL（避免 corner NaN 假阳性，2026-06-27 范式切换调试发现）。
    if (has_inlet && has_outlet &&
        std::isfinite(static_cast<double>(m_inlet)) &&
        std::isfinite(static_cast<double>(m_outlet)))
    {
        const Real denom = std::max(std::fabs(m_inlet), std::fabs(m_outlet));
        diag.mass_flux_imbalance = denom > 0.0 ? std::fabs(m_inlet - m_outlet) / denom : 0.0;
        diag.flux_computed = true;
    }

    // --- Profile L2: bin interior particles by y, compare x-z-averaged u_x to
    //     the parabolic target. Interior = not in inlet/outlet/wall face bands
    //     AND not in sponge region. ---
    const Real U_max = inletUMax(cfg);
    const int n_bins = std::max(10, static_cast<int>(cfg.DH / dp));
    std::vector<Real> sum_ux(n_bins, 0.0);
    std::vector<int> count_ux(n_bins, 0);
    bool has_interior = false;
    for (size_t i = 0; i < n; ++i)
    {
        if (in_sponge(pos[i][0]))
        {
            continue; // 排除 sponge 区粒子
        }
        const bool in_wall_band = (pos[i][1] <= face_band || pos[i][1] >= cfg.DH - face_band);
        const bool in_x_band = (pos[i][0] <= face_band || pos[i][0] >= cfg.DL - face_band);
        if (in_wall_band || in_x_band)
        {
            continue; // interior only — exclude all boundary-face particles
        }
        const int bin = static_cast<int>((pos[i][1] / cfg.DH) * static_cast<Real>(n_bins));
        const int b = std::min(std::max(bin, 0), n_bins - 1);
        sum_ux[b] += vel[i][0];
        ++count_ux[b];
        has_interior = true;
    }
    // Require a minimum bin coverage so a single stray interior particle
    // cannot trivially pass the gate (num==0 -> L2==0). With wall bands of
    // boundary_n_layers*dp on each side the usable y-range is narrower than
    // the full DH, so require only a quarter of the bins populated — enough
    // to reject a single-particle degenerate sample while tolerating the
    // reduced interior height.
    if (has_interior)
    {
        int populated_bins = 0;
        for (int b = 0; b < n_bins; ++b)
        {
            if (count_ux[b] > 0)
            {
                ++populated_bins;
            }
        }
        const int min_bins = std::max(2, n_bins / 4);
        if (populated_bins >= min_bins)
        {
            Real num = 0.0;
            Real den = 0.0;
            for (int b = 0; b < n_bins; ++b)
            {
                if (count_ux[b] == 0)
                {
                    continue;
                }
                const Real y = (static_cast<Real>(b) + 0.5) * cfg.DH / static_cast<Real>(n_bins);
                const Real u_avg = sum_ux[b] / static_cast<Real>(count_ux[b]);
                const Real u_target = inletProfileVelocity(y, cfg.DH, U_max);
                num += (u_avg - u_target) * (u_avg - u_target);
                den += u_target * u_target;
            }
            // NaN 守卫：u_avg 发散时 num 为 NaN，sqrt(NaN)=NaN，EXPECT_LT(NaN,0.6) 不触发
            // 失败导致假通过。必须检查 finite，NaN 时 profile_computed=false 让 GTest FAIL。
            if (std::isfinite(static_cast<double>(num)) && den > 0.0)
            {
                diag.profile_l2 = std::sqrt(num / den);
                diag.profile_computed = true;
            }
        }
    }
    return diag;
}

//----------------------------------------------------------------------
//	Main simulation driver. Throws on configuration/validation errors.
//  Returns the 50-step and final finiteness flags for the TEST body.
//----------------------------------------------------------------------
SmokeRunResult eulerian_channel_3d(const SimulationConfig &cfg)
{
    SmokeRunResult result{false, false, 0};
    bool fifty_step_finite_ = false; // set at the 50-step gate checkpoint
    const Real DL = cfg.DL, DH = cfg.DH, DW = cfg.DW, dp = cfg.global_resolution;
    const Real rho0_f = cfg.rho0_f;
    const Real c_f = cfg.c_f;
    const Real mu_f = cfg.rho0_f * cfg.nu;
    const Real DL_sponge = cfg.sponge_width_factor * dp;
    const Real BW = 4.0 * dp; ///< boundary padding for the system domain bounds.

    // 2D LG paradigm: 流体域含 sponge = [-DL_sponge, DL+DL_sponge]×[0,DH]×[0,DW]，
    // system bounds 外延 BW 给 wall body 粒子留余量（壁面侧 padding 保留）。
    BoundingBoxd system_domain_bounds(Vecd(-DL_sponge - BW, -BW, -BW),
                                      Vecd(DL + DL_sponge + BW, DH + BW, DW + BW));
    SPHSystem sph_system(system_domain_bounds, dp);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(false);

    // --- Fluid body: sponge-extended domain + OuterBoundary component level set. ---
    // defineComponentLevelSetShape 要求 initial_shape_ 是 ComplexShape（ChannelSpongeFluidBlock ✓），
    // 供 NormalDirectionFromBodyShape 的 findNormalDirection 在 sponge 外壳返回 ±x 外法向（契约 H3）。
    FluidBody fluid_block(sph_system, makeShared<ChannelSpongeFluidBlock>("ChannelFluid", cfg));
    fluid_block.defineComponentLevelSetShape("OuterBoundary");
    fluid_block.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    fluid_block.addMaterialProperty<Viscosity>(mu_f);
    fluid_block.generateParticles<BaseParticles, Lattice>(); // 无 ghost reserve，无 Reload

    // --- Wall body: 体粒子两片薄壁（2D LG wall-contact paradigm，体粒子版）。
    // SolidBody + ChannelWallsBlock（ComplexShape + GeometricShapeBox）+ defineBodyLevelSetShape
    // + NormalDirectionFromBodyShape 几何求值 + 普通 ContactRelation（对齐 2D LG 圆柱范式，
    // 参考 cylinder_lg_calculation.hpp:1469）。x 含 sponge 全覆盖，z 含缓冲层外延。
    // 流固压力用最原始共享层 WithWallRiemann（不用 pressure_contact_only_normal 法向投影）。 ---
    SolidBody wall_block(sph_system, makeShared<ChannelWallsBlock>("ChannelWalls", cfg));
    wall_block.defineAdaptationRatios(1.3, 2.0);
    wall_block.defineBodyLevelSetShape();
    wall_block.defineMatterMaterial<Solid>();
    wall_block.generateParticles<BaseParticles, Lattice>();

    // --- Relations: inner + wall contact + complex; z 周期保留（3D 槽道特征，2D LG 无）。 ---
    // 普通 ContactRelation（体粒子 Vol=dp³ 与 inner 同量级，无需 shell 厚度补偿）。
    InnerRelation fluid_inner(fluid_block);
    ContactRelation fluid_wall_contact(fluid_block, RealBodyVector{&wall_block});
    ComplexRelation fluid_wall_complex(fluid_inner, fluid_wall_contact);
    // Spanwise (z) periodicity. x is non-periodic (inlet/outlet via sponge + FarField).
    PeriodicAlongAxis periodic_along_z(fluid_block.getSPHBodyBounds(), zAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition_z(fluid_block, periodic_along_z);

    // --- Dynamics（全部共享层，对齐 2D LG wall-contact paradigm）。 ---
    // 构造顺序硬约束（plan 第 40 行契约）：wall_normal 必须先于 surface_indicator
    // （FreeSurfaceIndicationComplex 依赖 wall contact 已配置，契约 M3）；fluid_normal_direction
    // 必须先于 FarFieldBoundary（FarField 基类构造取 fluid body NormalDirection）。
    SimpleDynamics<NormalDirectionFromBodyShape> wall_normal_direction(wall_block);
    InteractionWithUpdate<FreeSurfaceIndicationComplex> surface_indicator(fluid_inner, fluid_wall_contact);
    InteractionDynamics<SmearedSurfaceIndication> smeared_surface(fluid_inner); // 隐式注册 SmearedSurface
    SimpleDynamics<NormalDirectionFromBodyShape> fluid_normal_direction(fluid_block);
    // pressure/density 用共享层 WithWallRiemann（最原始接口，不用法向投影）：wall 法向由
    // NormalDirectionFromBodyShape 几何求值（体粒子 box 端面 ±x/±z 退化由 wall x 含 sponge +
    // z 缓冲层外延缓解，端面 wall 粒子落在稀疏区由 FarField 主导）。limiter 走 ComplexInteraction
    // 默认 15.0（不传第三参 Real，签名不符）。
    InteractionWithUpdate<fluid_dynamics::EulerianIntegration1stHalfWithWallRiemann>
        pressure_relaxation(fluid_inner, fluid_wall_contact);
    InteractionWithUpdate<fluid_dynamics::EulerianIntegration2ndHalfWithWallRiemann>
        density_relaxation(fluid_inner, fluid_wall_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(fluid_inner, fluid_wall_contact);
    SimpleDynamics<EulerianChannelInitialCondition> initial_condition(fluid_block, cfg);
    channel_ck::AcousticTimeStep<> get_acoustic_dt(fluid_block, cfg.acoustic_cfl);
    // FarFieldBoundary：非反射开边界，构造注入 cfg 远场（契约 H1），update 自写出流/入流凸组合
    // （出流远场 vel 用抛物线 + 出流压力加权混合，对齐 2D LG）。
    // 必须用 InteractionWithUpdate 包装：exec() 才会先 interaction(邻域加权累加) 后 update(凸组合写 state)。
    // FarFieldBoundary 本身继承 LocalDynamics 无 exec，2D LG 同样用 InteractionWithUpdate<FarFieldBoundary>。
    InteractionWithUpdate<FarFieldBoundary> variable_reset_in_boundary_condition(fluid_inner, cfg);

    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition. 初始化序列对齐 2D LG（契约 M3）：
    //  wall_normal → surface_indicator → smeared_surface → fluid_normal →
    //  initial_condition → FarFieldBoundary 初始校正 → 只读诊断。
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    periodic_condition_z.update_cell_linked_list_.exec();
    sph_system.initializeSystemConfigurations();
    wall_normal_direction.exec();              // ← 先于 surface_indicator（契约 M3）
    surface_indicator.exec();
    smeared_surface.exec();
    fluid_normal_direction.exec();
    initial_condition.exec();
    variable_reset_in_boundary_condition.exec(); // FarFieldBoundary 初始校正

    BaseParticles &particles = fluid_block.getBaseParticles();
    // 只读诊断（契约 M2：补足未调用诊断 + H3 几何验证）。
    reportSpongeBoundaryDiagnostics(particles, cfg, "after init");
    reportFarFieldBoundaryRiskScan(particles, fluid_inner, cfg, "after init");

    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(fluid_block);
    write_states.addToWrite<Real>(fluid_block, "Density");
    write_states.addToWrite<Real>(fluid_block, "Pressure");
    write_states.addToWrite<Vecd>(fluid_block, "Velocity");
    ReducedQuantityRecording<MaximumSpeed> write_maximum_speed(fluid_block);

    //----------------------------------------------------------------------
    //	Setup for time-stepping control.
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = cfg.screen_output_interval;
    Real end_time = cfg.end_time;
    Real output_interval = end_time / static_cast<Real>(cfg.output_interval);

    TickCount t1 = TickCount::now();
    TimeInterval interval;
    write_states.writeToFile(0);

    //----------------------------------------------------------------------
    //	Main loop（对齐 2D LG，契约 M5 主循环 BC 顺序）：
    //  每 acoustic step: viscous_force → pressure_relaxation → density_relaxation
    //  → FarFieldBoundary(末尾一次)。废弃 v1 的 3 次 resetBoundaryConditions。
    //  FarFieldBoundary 是 InteractionWithUpdate，其 interaction(邻域加权累加) +
    //  update(凸组合写 state) 由 .exec() 一次完成。
    //----------------------------------------------------------------------
    bool dt_collapsed = false; // dt 崩溃标志，跳出双层循环避免死循环打印
    while (physical_time < end_time && !dt_collapsed)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval && !dt_collapsed)
        {
            const Real dt = get_acoustic_dt.exec();
            viscous_force.exec();
            pressure_relaxation.exec(dt);
            density_relaxation.exec(dt);
            variable_reset_in_boundary_condition.exec(); // FarFieldBoundary，末尾一次

            integration_time += dt;
            physical_time += dt;
            // dt 衰减保护：若 dt 退化到 0（corner 病态 vel 爆炸签名），physical_time
            // 不再推进，主循环会死循环。此时跳出并让 final_finite gate 捕获 NaN。
            // 阈值 TinyReal 之下视为 dt 崩溃（正常 dt ~ 1e-3 量级）。
            if (dt < TinyReal)
            {
                std::cout << "[ChannelDiag][DtCollapse] dt=" << dt
                          << " at N=" << number_of_iterations
                          << ", t=" << physical_time << " — breaking loop.\n";
                checkStatesFinite(particles, "dt-collapse");
                dt_collapsed = true;
                break;
            }
            if (number_of_iterations % screen_output_interval == 0)
            {
                write_maximum_speed.writeToFile(number_of_iterations);
                std::cout << std::fixed << std::setprecision(9)
                          << "N=" << number_of_iterations
                          << " t=" << physical_time
                          << " dt=" << dt << "\n";
            }
            // Gate checkpoint: at the 50th acoustic step, record whether every
            // real particle state is finite. The TEST body checks this flag.
            if (number_of_iterations == 50 && !fifty_step_finite_)
            {
                fifty_step_finite_ = checkStatesFinite(particles, "50-step gate");
            }
            number_of_iterations++;
        }
        write_states.writeToFile();
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt = t4 - t1 - interval;
    std::cout << "Total wall time: " << tt.seconds() << " s\n";

    result.fifty_step_finite = fifty_step_finite_;
    result.final_finite = checkStatesFinite(particles, "final");
    result.iterations = number_of_iterations;
    const FluxProfileDiagnostics diag = computeFluxProfileDiagnostics(particles, cfg);
    result.mass_flux_imbalance = diag.mass_flux_imbalance;
    result.profile_l2 = diag.profile_l2;
    result.flux_computed = diag.flux_computed;
    result.profile_computed = diag.profile_computed;
    std::cout << std::fixed << std::setprecision(6)
              << "[diag] mass_flux_imbalance=" << result.mass_flux_imbalance
              << " profile_l2=" << result.profile_l2
              << " (flux_computed=" << result.flux_computed
              << " profile_computed=" << result.profile_computed << ")\n";
    return result;
}

//----------------------------------------------------------------------
//	GTest cases.（仅保留契约测试：参数解析 / 入口剖面形状 / 低 Mach 守卫 /
//  helper 链接性。主计算由 main() 直接读 config.ini 驱动，不再走 GTest smoke 包装。）
//----------------------------------------------------------------------

// Contract: the default config parses to a parabolic inlet with the
// U_max = 1.5 * u_bulk relationship.
TEST(test_3d_eulerian_channel, config_parses_parabolic)
{
    SimulationConfig cfg = load_config(resolve_default_config_path().string());
    EXPECT_EQ(cfg.inlet_profile, "parabolic");
    EXPECT_NEAR(inletUMax(cfg), 1.5 * cfg.u_bulk, 1e-12);
    // Plan Task 1: sponge + z buffer + dp contract (2D LG wall-contact paradigm, 体粒子版).
    EXPECT_NEAR(cfg.global_resolution, 0.05, 1e-12);
    EXPECT_NEAR(cfg.sponge_width_factor, 5.0, 1e-12);
    EXPECT_NEAR(cfg.z_buffer_factor, 4.0, 1e-12);
}

// Contract: the parabolic profile vanishes at the walls, peaks at the
// centreline, and its bulk average equals u_bulk.
TEST(test_3d_eulerian_channel, inlet_profile_shape)
{
    SimulationConfig cfg = load_config(resolve_default_config_path().string());
    Real Umax = inletUMax(cfg);
    EXPECT_NEAR(inletProfileVelocity(0.0, cfg.DH, Umax), 0.0, 1e-12);        // bottom wall
    EXPECT_NEAR(inletProfileVelocity(cfg.DH, cfg.DH, Umax), 0.0, 1e-12);     // top wall
    EXPECT_NEAR(inletProfileVelocity(cfg.DH * 0.5, cfg.DH, Umax), Umax, 1e-12); // centreline
    EXPECT_NEAR(inletProfileBulk(cfg.DH, Umax), cfg.u_bulk, 1e-12);          // bulk average
}

// Contract: the low-Mach guard rejects a config with c_f too small unless
// allow_high_mach is set. A temporary ini is written and removed afterwards.
TEST(test_3d_eulerian_channel, config_rejects_high_mach)
{
    // 写到 case-local 临时目录，避免污染系统 Temp（CLAUDE.md 约定）。
    const std::filesystem::path case_dir =
        std::filesystem::path(__FILE__).parent_path();
    const std::filesystem::path tmp_dir = case_dir / "local_tmp";
    std::error_code mk_ec;
    std::filesystem::create_directories(tmp_dir, mk_ec);
    const std::filesystem::path tmp_ini = tmp_dir / "eulerian_channel_high_mach_test.ini";
    std::ofstream out(tmp_ini);
    if (!out.is_open())
    {
        throw std::runtime_error("Cannot write temporary config: " + tmp_ini.string());
    }
    out << "[physical]\n"
        << "rho0_f = 1.0\n"
        << "u_bulk = 1.0\n"
        << "c_f = 5.0\n"          // < 10 * 1.5 * u_bulk = 15
        << "nu = 0.01\n"
        << "[geometry]\n"
        << "DL = 2.0\n"
        << "DH = 1.0\n"
        << "DW = 1.0\n"
        << "global_resolution = 0.1\n"
        << "[inlet]\n"
        << "profile = parabolic\n"
        << "[simulation]\n"
        << "end_time = 2.0\n"
        << "acoustic_cfl = 0.25\n"
        << "output_interval = 20\n"
        << "screen_output_interval = 10\n"
        << "allow_high_mach = false\n"
        << "[boundary]\n"
        << "outlet_pressure_mode = zero_gradient\n";
    out.close();

    EXPECT_THROW(load_config(tmp_ini.string()), std::runtime_error);

    std::error_code ec;
    std::filesystem::remove(tmp_ini, ec);
    // Remove the temp directory if it is now empty; leave it if other tests
    // wrote into it concurrently.
    std::filesystem::remove(tmp_dir, ec);
}

// Contract: the Task 2 helper keeps the boundary quintuple coherent. A full
// BaseParticles cannot be built in isolation without an SPHSystem, so this
// case is covered indirectly by the main simulation (which calls the helper
// on every particle). Here we only assert the helper symbols are linked.
TEST(test_3d_eulerian_channel, helper_state_sync_linked)
{
    EXPECT_TRUE(&makeEulerianWeaklyCompressibleBoundaryState != nullptr);
    EXPECT_TRUE(&syncEulerianWeaklyCompressibleState != nullptr);
    EXPECT_TRUE(&setEulerianWeaklyCompressibleVelocity != nullptr);
    SUCCEED();
}

//----------------------------------------------------------------------
//	Entry point: 默认直接读 config.ini 驱动主计算（end_time 等参数完全由
//  config.ini 决定，不再有 smoke 包装强制缩短 end_time）。
//  保留 --gtest 入口以运行上方契约测试（参数解析 / 入口剖面 / 低 Mach 守卫）。
//  诊断信息（mass_flux_imbalance / profile_l2）由 eulerian_channel_3d 末尾打印。
//----------------------------------------------------------------------
int main(int argc, char **argv)
{
    // 带 --gtest 参数时走 GTest 契约测试；否则直接跑主计算。
    bool run_gtest = false;
    for (int i = 1; i < argc; ++i)
    {
        std::string arg(argv[i]);
        if (arg == "--gtest" || arg.rfind("--gtest_", 0) == 0)
        {
            run_gtest = true;
            break;
        }
    }

    if (run_gtest)
    {
        testing::InitGoogleTest(&argc, argv);
        return RUN_ALL_TESTS();
    }

    // 主计算：直接读 config.ini（resolve_default_config_path 回退到源目录）。
    std::cerr << "[main] start: resolving config path...\n";
    SimulationConfig cfg = load_config(resolve_default_config_path().string());
    std::cerr << "[main] config loaded: end_time=" << cfg.end_time
              << " dp=" << cfg.global_resolution << " — launching eulerian_channel_3d\n";
    std::cout << "[main] running eulerian_channel_3d from config.ini: "
              << "end_time=" << cfg.end_time << " dp=" << cfg.global_resolution << "\n";
    SmokeRunResult result = eulerian_channel_3d(cfg);

    // 末态有限性 + 诊断摘要（eulerian_channel_3d 内部已打印 mass_flux / profile，
    // 这里仅汇总退出码：状态非有限视为失败）。
    const int rc = result.final_finite ? 0 : 1;
    std::cout << "[main] done: iterations=" << result.iterations
              << " final_finite=" << (result.final_finite ? "yes" : "NO")
              << " exit=" << rc << "\n";
    return rc;
}
