/**
 * @file 	eulerian_channel_geometry.hpp
 * @brief 	Geometry and boundary-condition layer for the 3D Eulerian channel flow.
 * @details Fluid domain sponge-extended to [-DL_sponge, DL+DL_sponge] x [0,DH] x [0,DW].
 *          Top/bottom (y) are no-slip volumetric walls (ChannelWallsBlock +
 *          NormalDirectionFromBodyShape + ContactRelation + WithWallRiemann).
 *          Inlet/outlet (x) are handled by the sponge layers plus a non-reflective
 *          FarFieldBoundary (shared-layer NonReflectiveBoundaryCorrection, cfg-driven).
 *          Spanwise (z) periodicity is handled by the SPHSystem periodic BC.
 *          The earlier ghost-mirror + boundary_type-label (3/10/5) paradigm and the
 *          v1 self-written correction class have been retired (2026-06-27) and removed.
 */
#ifndef EULERIAN_CHANNEL_GEOMETRY_H
#define EULERIAN_CHANNEL_GEOMETRY_H

#include "eulerian_channel_data.hpp"
#include "eulerian_open_boundary.h" // shared-layer state-sync helper (syncEulerianWeaklyCompressibleState)
#include "sphinxsys.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace SPH
{

//----------------------------------------------------------------------
//  Sponge-extended fluid domain aligned with the 2D LG paradigm.
//  流体物理核心 [0,DL]×[0,DH]×[0,DW]，sponge 仅在 x 方向（入口/出口）外延，
//  厚度 DL_sponge = sponge_width_factor * dp（默认 5×dp，对齐 2D LG）。
//  完整流体域 = [-DL_sponge, DL+DL_sponge]×[0,DH]×[0,DW]，sponge 区有真实粒子
//  （Lattice 生成），使入口/出口面粒子拥有完整核支持，corner 粒子不再三面核缺失。
//  y 方向（上下壁）不加 sponge，由真实 wall body 接管；z 方向保持周期。
//  关键：子 shape 命名 "OuterBoundary"，供 defineComponentLevelSetShape +
//  NormalDirectionFromBodyShape 的 findNormalDirection 使用，使 sponge 外壳
//  粒子获得严格沿 ±x 的外法向（入口 n_x<0、出口 n_x>0），避免基类入流判定
//  fabs(n_y)>fabs(n_x) 在角落误判（契约 H3）。
//----------------------------------------------------------------------
class ChannelSpongeFluidBlock : public ComplexShape
{
  public:
    explicit ChannelSpongeFluidBlock(const std::string &shape_name, const SimulationConfig &cfg)
        : ComplexShape(shape_name)
    {
        const Real DL_sponge = cfg.sponge_width_factor * cfg.global_resolution;
        // 域 [-DL_sponge, DL+DL_sponge]×[0,DH]×[0,DW]：box 中心 x = DL/2（物理核心中点）。
        const Vecd halfsize(0.5 * (cfg.DL + 2.0 * DL_sponge), 0.5 * cfg.DH, 0.5 * cfg.DW);
        const Vecd translation(0.5 * cfg.DL, 0.5 * cfg.DH, 0.5 * cfg.DW);
        add<GeometricShapeBox>(Transform(translation), halfsize, "OuterBoundary");
    }
};

//----------------------------------------------------------------------
//  Wall body geometry（2D LG wall-contact paradigm，体粒子版，2026-06-27）。
//  两片薄壁：下壁 y∈[-dp/2, dp/2]、上壁 y∈[DH-dp/2, DH+dp/2]，厚度 dp（体粒子，Vol=dp³）。
//  SolidBody + ComplexShape(GeometricShapeBox) + defineBodyLevelSetShape +
//  NormalDirectionFromBodyShape 几何求值 + 普通 ContactRelation（对齐 2D LG 圆柱范式，
//  参考 cylinder_lg_calculation.hpp:1469 Cylinder : MultiPolygonShape）。
//  3D 无 MultiPolygonShape（仅 for_2D_build），用 GeometricShapeBox 等价表达平面薄壁。
//
//  ★ x 范围覆盖 [-DL_sponge, DL+DL_sponge]（含 sponge 全覆盖，用户决策）：
//  sponge 区有真实流体粒子（ChannelSpongeFluidBlock），入口/出口 corner 流体粒子
//  （sponge 外壳 + 紧邻壁面）必须有充足 wall contact 邻居——wall x 覆盖 sponge 让 corner
//  粒子 wall contact 不缺。wall x 端面粒子（x=±DL_sponge 处）法向退化为 ±x，但这些端面
//  wall 粒子落在 sponge 外侧（流体粒子稀疏区），其 wall contact 贡献由 FarField 主导，
//  端面退化影响被 sponge buffer 吸收。
//
//  ★ z 范围覆盖 [-DZ_buffer, DW+DZ_buffer]（z 缓冲层外延，用户决策）：
//  z 是周期方向，fluid 域 z 保持 [0,DW] 周期不变（不破坏周期长度 DW）。wall z 外延
//  DZ_buffer = z_buffer_factor * dp（默认 4×dp）到周期 ghost 区，让流体 z 端面粒子
//  （z≈0/DW）的 wall contact 邻居经周期 CellLinkedList 映射到壁面内部粒子（法向 ±y），
//  消除 z 端面 wall 粒子法向退化（±z）对 z 端面流体粒子的影响。
//
//  命名 "ChannelWalls" 供 defineBodyLevelSetShape 与 NormalDirectionFromBodyShape 引用。
//----------------------------------------------------------------------
class ChannelWallsBlock : public ComplexShape
{
  public:
    explicit ChannelWallsBlock(const std::string &shape_name, const SimulationConfig &cfg)
        : ComplexShape(shape_name)
    {
        const Real dp = cfg.global_resolution;
        const Real DL_sponge = cfg.sponge_width_factor * dp;
        const Real DZ_buffer = cfg.z_buffer_factor * dp;
        // x 含 sponge 全覆盖 [-DL_sponge, DL+DL_sponge]；z 含缓冲层 [-DZ_buffer, DW+DZ_buffer]。
        const Vecd wall_halfsize(0.5 * (cfg.DL + 2.0 * DL_sponge), 0.5 * dp,
                                 0.5 * cfg.DW + DZ_buffer);
        // 下壁：中心 (DL/2, 0, DW/2) —— x 中心物理核心中点，含 sponge 对称。
        const Vecd lower_center(0.5 * cfg.DL, 0.0, 0.5 * cfg.DW);
        add<GeometricShapeBox>(Transform(lower_center), wall_halfsize);
        // 上壁：中心 (DL/2, DH, DW/2)
        const Vecd upper_center(0.5 * cfg.DL, cfg.DH, 0.5 * cfg.DW);
        add<GeometricShapeBox>(Transform(upper_center), wall_halfsize);
    }
};

//----------------------------------------------------------------------
//	Initial condition: parabolic profile across the whole domain.
//  Matching the inlet reduces startup transients.
//----------------------------------------------------------------------
class EulerianChannelInitialCondition : public fluid_dynamics::FluidInitialCondition
{
  public:
    explicit EulerianChannelInitialCondition(SPHBody &sph_body, const SimulationConfig &cfg)
        : fluid_dynamics::FluidInitialCondition(sph_body),
          cfg_(cfg),
          U_max_(inletUMax(cfg)),
          // registerStateVariableData is idempotent: returns the existing pointer
          // if the variable was already registered by an earlier dynamics, and
          // registers it otherwise. This decouples IC construction order from
          // other fluid dynamics (mirrors test_3d_FVM_incompressible_channel_flow).
          rho_(particles_->registerStateVariableData<Real>("Density")),
          p_(particles_->registerStateVariableData<Real>("Pressure")),
          mass_(particles_->registerStateVariableData<Real>("Mass")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          mom_(particles_->registerStateVariableData<Vecd>("Momentum")) {};

    void update(size_t index_i, Real dt)
    {
        rho_[index_i] = cfg_.rho0_f;
        p_[index_i] = 0.0;
        const Real u_x = inletProfileVelocity(pos_[index_i][1], cfg_.DH, U_max_);
        vel_[index_i][0] = u_x;
        vel_[index_i][1] = 0.0;
        vel_[index_i][2] = 0.0;
        mass_[index_i] = rho_[index_i] * Vol_[index_i];
        mom_[index_i] = mass_[index_i] * vel_[index_i];
    }

  protected:
    const SimulationConfig &cfg_;
    Real U_max_;
    Real *rho_, *p_, *mass_, *Vol_;
    Vecd *mom_;
};

//----------------------------------------------------------------------
//  2D LG paradigm: FarFieldBoundary 非反射开边界（3D 槽道语义 + cfg 本地化）。
//  继承 fluid_dynamics::NonReflectiveBoundaryCorrection，复用基类邻域加权累加
//  （interaction: inner_weight_summation / rho_average / vel_normal_average /
//   vel_tangential_average，过滤 indicator==1 邻居）。
//
//  构造本地化（契约 H1）：2D LG 的 FarFieldBoundary 真实签名只接 inner_relation，
//  远场量靠全局变量 rho0_f/c_f/U_f 硬编码；3D channel 用 SimulationConfig 无这些
//  全局变量，故构造签名扩展为 (inner_relation, cfg)，从 cfg 注入远场目标。
//
//  ★ 出流远场 vel 用抛物线（用户决策，对齐 2D LG 入流抛物线 + 出流压力加权混合）：
//  vel_farfield 不再是均匀 (U_max,0,0)，而是逐粒子按 pos_i[1] 的抛物线 u_x(y)。
//  sponge 硬约束粒子直接强制抛物线 vel；出流/入流凸组合粒子用抛物线 vel_farfield
//  作远场目标。rho_farfield=rho0 恒定（出流压力加权：rho = rho_avg·w + rho0·(1-w)）。
//
//  update 四分支（对齐 cylinder_lg_calculation.hpp:1493-1564，自写不委托基类）：
//  ① sponge outer buffer 硬约束：x<0 或 x>DL 强制远场值（rho=rho0, vel=抛物线, p=EOS）。
//  ② 壁面附近粒子（y≤wall_band 或 y≥DH-wall_band）跳过，交给 wall contact。
//  ③ 非边界粒子（indicator!=1 且 smeared_surface!=1）跳过。
//  ④ 边界粒子按入流/出流×亚声速/超声速凸组合（出流切向零梯度透传 vel_tangential_average_）。
//
//  依赖前置（构造顺序硬约束）：基类构造时 getVariableDataByName 取
//  Indicator / SmearedSurface / NormalDirection，必须由 FreeSurfaceIndicationComplex
//  + SmearedSurfaceIndication + NormalDirectionFromBodyShape 先注册（main 中构造顺序保证）。
//  SmearedSurface 由 SmearedSurfaceIndication 隐式注册，无需显式 registerStateVariableData。
//----------------------------------------------------------------------
class FarFieldBoundary : public fluid_dynamics::NonReflectiveBoundaryCorrection
{
  public:
    explicit FarFieldBoundary(BaseInnerRelation &inner_relation, const SimulationConfig &cfg)
        : fluid_dynamics::NonReflectiveBoundaryCorrection(inner_relation),
          DL_(cfg.DL),
          DH_(cfg.DH),
          U_max_(inletUMax(cfg)),
          wall_band_(static_cast<Real>(cfg.boundary_n_layers) * cfg.global_resolution)
    {
        // 远场目标基线：rho_farfield=rho0，sound_speed=c_f。
        // vel_farfield_ 仅作基类 interaction 的占位；逐粒子远场速度由 farfieldVelocity(i) 给出。
        rho_farfield_ = cfg.rho0_f;
        sound_speed_ = cfg.c_f;
        vel_farfield_ = Vecd(U_max_, 0.0, 0.0);
        // 共享层状态同步 helper（rho/p/mass/mom 在写入 vel/rho 后一次性拉齐）。
        boundary_state_ = makeEulerianWeaklyCompressibleBoundaryState(*particles_);
    };
    virtual ~FarFieldBoundary() {};

    // 远场抛物线速度 u_x(y) = U_max * 4 * eta * (1 - eta)，eta = y / DH。
    // 本 case 自写 update() 的辅助方法（非基类虚函数，基类上游版用 vel_farfield_ 成员）。
    Vecd farfieldVelocity(size_t index_i) const
    {
        const Real u_x = inletProfileVelocity(pos_[index_i][1], DH_, U_max_);
        return Vecd(u_x, 0.0, 0.0);
    }

    // 2D LG 范式 FarField update 四分支（对齐 cylinder_lg_calculation.hpp:1493-1564，
    // 自写出流/入流凸组合。基类 update 已通过 farfieldVelocity(i) 取逐粒子远场剖面，
    // 但本 case 额外需要 ① sponge 硬约束、② 壁面跳过、NaN clamp 等 case-specific 防护，
    // 故仍自写 update（防护留在工况内，不上提基类））。
    void update(size_t index_i, Real dt = 0.0)
    {
        const Vecd &pos_i = pos_[index_i];
        // ① sponge outer buffer 硬约束：物理域外（x<0 或 x>DL）强制远场值（抛物线 vel）。
        if (pos_i[0] < Real(0.0) || pos_i[0] > DL_)
        {
            rho_[index_i] = rho_farfield_;
            vel_[index_i] = farfieldVelocity(index_i);
            syncEulerianWeaklyCompressibleState(boundary_state_, index_i);
            return;
        }
        // ② 壁面附近粒子跳过（交给 wall contact + inner 积分器，防 FarField vel 与 wall 无滑移冲突）。
        if (pos_i[1] <= wall_band_ || pos_i[1] >= DH_ - wall_band_)
        {
            return;
        }
        // ③ 非边界粒子跳过。
        if (indicator_[index_i] != 1 && smeared_surface_[index_i] != 1)
        {
            return;
        }
        // ④ 边界粒子按入流/出流×亚声速/超声速凸组合。
        const Vecd &n_i = n_[index_i];
        const Vecd vel_ff = farfieldVelocity(index_i);
        const Real velocity_farfield_normal = vel_ff.dot(n_i);
        const Real velocity_boundary_normal = vel_[index_i].dot(n_i);
        const Real abs_vbn = std::fabs(velocity_boundary_normal);
        // 凸组合权重 w = clamp(inner_weight_summation, 0, 1)，NaN 守卫。
        Real w = inner_weight_summation_[index_i];
        if (!std::isfinite(static_cast<double>(w)))
        {
            w = Real(0.0);
        }
        else
        {
            w = w < Real(0.0) ? Real(0.0) : (w > Real(1.0) ? Real(1.0) : w);
        }
        const Real one_minus_w = Real(1.0) - w;

        // 入流判定（严格对齐基类 NonReflectiveBoundaryCorrection::interaction 的
        // `n_x<=0 || |n_y|>|n_x|`，3D 增 |n_z|>|n_x|）：保证子类 update 分支与基类
        // interaction 填充的邻域平均值（vel_tangential_average_ 等）读写匹配，避免
        // corner 粒子（n_y 或 n_z 主导但 n_x>0）读陈旧值。当前进入④的粒子法向 ±x 主导，
        // 此判定与简化版 n_x<=0 等价，但防御未来 wall_band/face_band 调整。
        const Real abs_nx = std::fabs(n_i[0]);
        const bool inflow = n_i[0] <= Real(0.0) ||
                            std::fabs(n_i[1]) > abs_nx ||
                            std::fabs(n_i[2]) > abs_nx;

        if (abs_vbn >= sound_speed_)
        {
            // 超声速：入流强制远场值；出流用内部平均（零梯度透传）。
            if (inflow)
            {
                rho_[index_i] = rho_farfield_;
                vel_[index_i] = vel_ff;
            }
            else
            {
                rho_[index_i] = rho_average_[index_i] + TinyReal;
                vel_[index_i] = vel_average_[index_i];
            }
            syncEulerianWeaklyCompressibleState(boundary_state_, index_i);
            return;
        }
        // 亚声速凸组合：rho = rho_avg·w + rho0·(1-w)；vel 法向凸组合。
        rho_[index_i] = rho_average_[index_i] * w + rho_farfield_ * one_minus_w;
        const Real vel_normal =
            vel_normal_average_[index_i] * w + velocity_farfield_normal * one_minus_w;
        if (inflow)
        {
            // 入流：切向用远场抛物线切向（vel_ff - vel_ff·n·n），强制抛物线剖面。
            vel_[index_i] = vel_normal * n_i + (vel_ff - velocity_farfield_normal * n_i);
        }
        else
        {
            // 出流：切向零梯度透传 vel_tangential_average_（对齐 2D LG 出流分支）。
            vel_[index_i] = vel_normal * n_i + vel_tangential_average_[index_i];
        }
        syncEulerianWeaklyCompressibleState(boundary_state_, index_i);
    }

  protected:
    Real DL_;          ///< 物理域流向长度，sponge 区判定（x<0 或 x>DL）。
    Real DH_;          ///< 物理通道高度，壁面附近判定 + 抛物线 eta。
    Real U_max_;       ///< 抛物线峰值 U_max = 1.5×u_bulk。
    Real wall_band_;   ///< 壁面附近带宽，该带内粒子交给 wall contact。
    EulerianWeaklyCompressibleBoundaryState boundary_state_; ///< rho/p/mass/Vol/vel/mom/EOS 聚合
};

//----------------------------------------------------------------------
//  2D LG paradigm: sponge 边界只读诊断（契约 M2 / H3 验证）。
//  统计 indicator 分布（sponge 外壳 indicator=1 数量、壁面附近 indicator=0 验证）
//  + NormalDirection 分布（入口 n_x<0、出口 n_x>0 验证，H3 几何保证）。
//  不修改任何 state，纯 cout。NormalDirection 未注册时安静 no-op。
//----------------------------------------------------------------------
inline void reportSpongeBoundaryDiagnostics(BaseParticles &particles, const SimulationConfig &cfg,
                                            const std::string &stage)
{
    Vecd *normal = particles.getVariableDataByName<Vecd>("NormalDirection");
    if (normal == nullptr)
    {
        std::cout << "[ChannelDiag][Sponge] " << stage
                  << ": NormalDirection not registered yet, skip." << std::endl;
        return;
    }
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    int *indicator = particles.getVariableDataByName<int>("Indicator");
    const Real DL_sponge = cfg.sponge_width_factor * cfg.global_resolution;
    const Real dp = cfg.global_resolution;
    const Real wall_band = static_cast<Real>(cfg.boundary_n_layers) * dp;

    const size_t total_particles = particles.TotalRealParticles();
    size_t indicator_one = 0;
    size_t inlet_sponge = 0, outlet_sponge = 0, wall_near = 0, interior = 0;
    size_t inlet_nx_negative = 0, outlet_nx_positive = 0;
    size_t corner_misjudged = 0; // H3: |n_y| > |n_x| 的外壳粒子（会被基类误判入流）

    for (size_t i = 0; i != total_particles; ++i)
    {
        if (indicator != nullptr && indicator[i] == 1)
        {
            ++indicator_one;
        }
        const Real x = pos[i][0];
        const Real y = pos[i][1];
        if (x < -0.5 * dp)
        {
            ++inlet_sponge; // 入口 sponge 区（x < 0）
        }
        else if (x > cfg.DL + 0.5 * dp)
        {
            ++outlet_sponge; // 出口 sponge 区（x > DL）
        }
        else if (y <= wall_band || y >= cfg.DH - wall_band)
        {
            ++wall_near;
        }
        else
        {
            ++interior;
        }

        // H3 验证：sponge 外壳粒子的法向应严格沿 ±x。
        const Real nx = normal[i][0];
        const Real ny = normal[i][1];
        if (std::fabs(nx) > 0.5) // 外壳粒子（法向显著）
        {
            if (nx < 0.0)
            {
                ++inlet_nx_negative;
            }
            else
            {
                ++outlet_nx_positive;
            }
            if (std::fabs(ny) > std::fabs(nx))
            {
                ++corner_misjudged;
            }
        }
    }

    std::cout << "[ChannelDiag][Sponge] " << stage
              << ": total=" << total_particles
              << ", indicator=1 count=" << indicator_one
              << ", inlet_sponge=" << inlet_sponge
              << ", outlet_sponge=" << outlet_sponge
              << ", wall_near=" << wall_near
              << ", interior=" << interior
              << ", DL_sponge=" << DL_sponge
              << "\n              inlet_nx<0=" << inlet_nx_negative
              << ", outlet_nx>0=" << outlet_nx_positive
              << ", corner_misjudged(|n_y|>|n_x|)=" << corner_misjudged
              << std::endl;
}

//----------------------------------------------------------------------
//  Core diagnostics (Task 0): 2D LG 风格诊断机制移植。
//  参照 tests/2d_examples/test_2d_eulerian_flow_around_cylinder_LG/
//  cylinder_lg_calculation.hpp:589-1027 的 5 件函数。
//  3D 槽道可用的变量子集：pos / Vol / rho / mass / p / vel / mom /
//  boundary_type / indicator / (n_ 在 Task 1 后可用)。
//----------------------------------------------------------------------

/** @brief 判断 Real 数值是否为有限值（NaN/Inf 防护）。 */
inline bool isFiniteRealEulerianChannel(Real value)
{
    return std::isfinite(static_cast<double>(value));
}

/** @brief 判断 SPHinXsys 向量每个分量是否为有限值。 */
inline bool isFiniteVectorEulerianChannel(const Vecd &value)
{
    for (int axis = 0; axis != Dimensions; ++axis)
    {
        if (!isFiniteRealEulerianChannel(value[axis]))
        {
            return false;
        }
    }
    return true;
}

/**
 * @brief 扫描流体粒子的 state quintuple（rho/p/mass/vel/mom）的有限性，
 *        第一个非有限粒子输出 stage 标签 + pos + 全字段。
 * @return true 当所有粒子 state 全有限；false 当首个非有限粒子出现。
 */
inline bool reportFluidFiniteState(BaseParticles &particles, const std::string &stage,
                                   size_t *first_bad_particle = nullptr)
{
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    Real *pressure = particles.getVariableDataByName<Real>("Pressure");
    Vecd *velocity = particles.getVariableDataByName<Vecd>("Velocity");
    Vecd *momentum = particles.getVariableDataByName<Vecd>("Momentum");
    const size_t total_particles = particles.TotalRealParticles();

    Real rho_min = std::numeric_limits<Real>::max();
    Real rho_max = -std::numeric_limits<Real>::max();
    Real p_min = std::numeric_limits<Real>::max();
    Real p_max = -std::numeric_limits<Real>::max();
    Real mass_min = std::numeric_limits<Real>::max();
    Real mass_max = -std::numeric_limits<Real>::max();
    Real speed_max = Real(0.0);
    size_t rho_count = 0, p_count = 0, mass_count = 0;
    size_t first_bad = total_particles;

    for (size_t i = 0; i != total_particles; ++i)
    {
        const bool finite =
            isFiniteVectorEulerianChannel(pos[i]) && isFiniteRealEulerianChannel(vol[i]) &&
            isFiniteRealEulerianChannel(rho[i]) && isFiniteRealEulerianChannel(mass[i]) &&
            isFiniteRealEulerianChannel(pressure[i]) &&
            isFiniteVectorEulerianChannel(velocity[i]) && isFiniteVectorEulerianChannel(momentum[i]);
        if (!finite && first_bad == total_particles)
        {
            first_bad = i;
        }
        if (isFiniteRealEulerianChannel(rho[i]))
        {
            rho_min = SMIN(rho_min, rho[i]);
            rho_max = SMAX(rho_max, rho[i]);
            ++rho_count;
        }
        if (isFiniteRealEulerianChannel(pressure[i]))
        {
            p_min = SMIN(p_min, pressure[i]);
            p_max = SMAX(p_max, pressure[i]);
            ++p_count;
        }
        if (isFiniteRealEulerianChannel(mass[i]))
        {
            mass_min = SMIN(mass_min, mass[i]);
            mass_max = SMAX(mass_max, mass[i]);
            ++mass_count;
        }
        if (isFiniteVectorEulerianChannel(velocity[i]))
        {
            speed_max = SMAX(speed_max, velocity[i].norm());
        }
    }

    if (first_bad != total_particles)
    {
        std::cout << "[ChannelDiag][State] " << stage
                  << ": rho=[" << (rho_count > 0 ? rho_min : Real(0.0)) << ", "
                  << (rho_count > 0 ? rho_max : Real(0.0)) << "]"
                  << ", p=[" << (p_count > 0 ? p_min : Real(0.0)) << ", "
                  << (p_count > 0 ? p_max : Real(0.0)) << "]"
                  << ", mass=[" << (mass_count > 0 ? mass_min : Real(0.0)) << ", "
                  << (mass_count > 0 ? mass_max : Real(0.0)) << "]"
                  << ", max|u|=" << speed_max
                  << ", first_bad_particle=" << first_bad << std::endl;
        std::cout << "[ChannelDiag][BadParticle] " << stage
                  << ": i=" << first_bad
                  << ", pos=(" << pos[first_bad][0] << ", " << pos[first_bad][1] << ", "
                  << pos[first_bad][2] << ")"
                  << ", Vol=" << vol[first_bad]
                  << ", rho=" << rho[first_bad]
                  << ", mass=" << mass[first_bad]
                  << ", p=" << pressure[first_bad]
                  << ", vel=(" << velocity[first_bad][0] << ", " << velocity[first_bad][1]
                  << ", " << velocity[first_bad][2] << ")"
                  << ", mom=(" << momentum[first_bad][0] << ", " << momentum[first_bad][1]
                  << ", " << momentum[first_bad][2] << ")" << std::endl;
    }
    if (first_bad_particle != nullptr)
    {
        *first_bad_particle = first_bad;
    }
    return first_bad == total_particles;
}

/**
 * @brief 单个非反射边界粒子的诊断快照记录。
 *        若 n_/indicator_ 未注册（Task 2 之前），active_boundary=false 直接返回。
 */
struct FarFieldBoundaryDiagnosticRecord
{
    size_t index = 0;
    size_t inner_neighbor_count = 0;
    bool active_boundary = false;
    bool inflow = false;
    bool subsonic = true;
    Real velocity_boundary_normal = 0.0;
    Real velocity_farfield_normal = 0.0;
    Real inner_weight_summation = 0.0;
    Real bounded_inner_weight_summation = 0.0;
    Real rho_average = 0.0;
    Real vel_normal_average = 0.0;
    Real rho_predicted = 0.0;
    Vecd vel_average = Vecd::Zero();
    Vecd vel_tangential_average = Vecd::Zero();
};

/**
 * @brief 为单个边界粒子重建远场分支预测（入流/出流/亚声速/超音速）。
 *        内部权重求和、邻域 ρ/vel 加权平均，用于定位 prediction 风险。
 *        若 NormalDirection 未注册则直接返回 inactive record。
 */
inline FarFieldBoundaryDiagnosticRecord computeFarFieldBoundaryDiagnostic(
    BaseParticles &particles, BaseInnerRelation &inner, size_t particle_index,
    const SimulationConfig &cfg)
{
    FarFieldBoundaryDiagnosticRecord record;
    record.index = particle_index;

    // NormalDirection 由 NormalDirectionFromBodyShape 写入。
    Vecd *normal = particles.getVariableDataByName<Vecd>("NormalDirection");
    if (normal == nullptr)
    {
        return record;
    }

    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real c_f = cfg.c_f;
    const Real U_max = inletUMax(cfg);
    const Real rho0_f = cfg.rho0_f;

    // 3D 槽道法向判定：x 法向 → 入流 (n_x<=0) 或 出流 (n_x>0)；y 法向 → 壁面，跳过。
    // 必须在读取 record.inflow 之前赋值，否则 inflow 默认为 false 导致
    // farfield_u_x 三元表达式恒取 0 分支（review BUG-1 修复）。
    record.inflow = normal[particle_index][0] <= Real(0.0);
    record.active_boundary =
        std::fabs(normal[particle_index][0]) > Real(0.5) ||
        std::fabs(normal[particle_index][1]) > Real(0.5);
    record.velocity_boundary_normal = vel[particle_index].dot(normal[particle_index]);
    record.subsonic = std::fabs(record.velocity_boundary_normal) < c_f;

    if (!record.active_boundary)
    {
        return record;
    }

    // 3D 槽道远场速度：入流段按 y 决定的抛物线 vel_farfield_x(y)，
    // 出流段 vel_farfield = 0（实际 correction 用 vel_self 零梯度透射）。
    // 此处只在 diagnostic 中用于显示，不参与 state 改写。
    const Real farfield_u_x =
        record.inflow ? inletProfileVelocity(pos[particle_index][1], cfg.DH, U_max) : Real(0.0);
    const Vecd farfield_velocity(farfield_u_x, Real(0.0), Real(0.0));
    record.velocity_farfield_normal = farfield_velocity.dot(normal[particle_index]);

    Real rho_sum = 0.0;
    Real vel_normal_sum = 0.0;
    Vecd vel_sum = Vecd::Zero();
    Vecd vel_tangential_sum = Vecd::Zero();
    const Neighborhood &nbh = inner.inner_configuration_[particle_index];
    for (size_t n = 0; n != nbh.current_size_; ++n)
    {
        const size_t j = nbh.j_[n];
        ++record.inner_neighbor_count;
        rho_sum += rho[j];
        vel_sum += vel[j];
        const Real v_n_j = vel[j].dot(normal[particle_index]);
        vel_normal_sum += v_n_j;
        vel_tangential_sum += vel[j] - v_n_j * normal[particle_index];
        record.inner_weight_summation += nbh.W_ij_[n] * vol[j];
    }
    // 孤立边界粒子（inner 邻居为 0）兜底：默认远场值，避免 risk scan 误报。
    if (record.inner_neighbor_count == 0)
    {
        record.rho_average = rho0_f;
        record.vel_normal_average = Real(0.0);
        record.vel_average = Vecd::Zero();
        record.vel_tangential_average = Vecd::Zero();
        record.bounded_inner_weight_summation = Real(0.0);
        record.rho_predicted = rho0_f; // 显式赋正值，让 risk scan 跳过（review EDGE-1 修复）
        return record;
    }
    const Real inv_n = Real(1.0) / static_cast<Real>(record.inner_neighbor_count);
    record.rho_average = rho_sum * inv_n;
    record.vel_normal_average = vel_normal_sum * inv_n;
    record.vel_average = vel_sum * inv_n;
    record.vel_tangential_average = vel_tangential_sum * inv_n;
    record.bounded_inner_weight_summation =
        SMIN(Real(1.0), SMAX(Real(0.0), record.inner_weight_summation));
    // 预测 ρ：subsonic → ρ_avg·w + ρ0·(1-w)（与 base update 一致）
    record.rho_predicted =
        record.subsonic
            ? record.rho_average * record.inner_weight_summation +
                  rho0_f * (Real(1.0) - record.inner_weight_summation)
            : (record.inflow ? rho0_f : record.rho_average + TinyReal);
    return record;
}

/**
 * @brief 输出单个粒子在指定阶段的远场诊断快照（含几何 + 邻域统计 + 预测 ρ）。
 *        若 NormalDirection 未注册则 no-op（避免 ghost 边界误报）。
 */
inline void writeFarFieldBoundaryDiagnostic(BaseParticles &particles, BaseInnerRelation &inner,
                                            size_t particle_index, const SimulationConfig &cfg,
                                            const std::string &stage)
{
    Vecd *normal = particles.getVariableDataByName<Vecd>("NormalDirection");
    if (normal == nullptr)
    {
        return; // NormalDirection 未注册，安静跳过
    }
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    Real *p = particles.getVariableDataByName<Real>("Pressure");
    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const FarFieldBoundaryDiagnosticRecord record =
        computeFarFieldBoundaryDiagnostic(particles, inner, particle_index, cfg);
    const std::string direction = record.inflow ? "inflow" : "outflow";
    const std::string acoustic = record.subsonic ? "subsonic" : "supersonic";
    std::cout << "[ChannelDiag][FarField] " << stage
              << ": i=" << particle_index
              << ", active=" << (record.active_boundary ? "true" : "false")
              << ", branch=" << direction << "/" << acoustic
              << ", pos=(" << pos[particle_index][0] << ", " << pos[particle_index][1]
              << ", " << pos[particle_index][2] << ")"
              << ", n=(" << normal[particle_index][0] << ", "
              << normal[particle_index][1] << ", " << normal[particle_index][2] << ")"
              << ", u_n=" << record.velocity_boundary_normal
              << ", u_farfield_n=" << record.velocity_farfield_normal
              << ", inner_neighbors=" << record.inner_neighbor_count
              << ", w=" << record.inner_weight_summation
              << ", w_bounded=" << record.bounded_inner_weight_summation
              << ", rho_average=" << record.rho_average
              << ", rho_predicted=" << record.rho_predicted
              << ", rho_current=" << rho[particle_index]
              << ", mass_current=" << mass[particle_index]
              << ", p_current=" << p[particle_index]
              << ", vel=(" << vel[particle_index][0] << ", " << vel[particle_index][1]
              << ", " << vel[particle_index][2] << ")" << std::endl;
}

/**
 * @brief 扫描全部流体粒子的远场预测，当 min_raw_predicted_rho<=0 时输出风险诊断。
 *        不修改任何状态，纯观测。NormalDirection 未注册时安静 no-op。
 */
inline void reportFarFieldBoundaryRiskScan(BaseParticles &particles, BaseInnerRelation &inner,
                                          const SimulationConfig &cfg, const std::string &stage)
{
    Vecd *normal = particles.getVariableDataByName<Vecd>("NormalDirection");
    if (normal == nullptr)
    {
        return;
    }
    const size_t total_particles = particles.TotalRealParticles();
    size_t active_count = 0;
    size_t min_pred_idx = total_particles;
    Real min_raw_predicted_rho = std::numeric_limits<Real>::max();
    for (size_t i = 0; i != total_particles; ++i)
    {
        const FarFieldBoundaryDiagnosticRecord record =
            computeFarFieldBoundaryDiagnostic(particles, inner, i, cfg);
        if (!record.active_boundary)
        {
            continue;
        }
        ++active_count;
        if (record.rho_predicted < min_raw_predicted_rho)
        {
            min_raw_predicted_rho = record.rho_predicted;
            min_pred_idx = i;
        }
    }
    if (active_count == 0 || min_raw_predicted_rho > Real(0.0))
    {
        return; // 无 active 边界或预测全正，不报警
    }
    std::cout << "[ChannelDiag][FarFieldScan] " << stage
              << ": active=" << active_count
              << ", min_raw_predicted_rho=" << min_raw_predicted_rho
              << ", min_pred_i=" << min_pred_idx << std::endl;
    if (min_pred_idx != total_particles)
    {
        writeFarFieldBoundaryDiagnostic(particles, inner, min_pred_idx, cfg,
                                       stage + " min_predicted");
    }
}

} // namespace SPH

#endif // EULERIAN_CHANNEL_GEOMETRY_H
