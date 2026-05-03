#ifndef CYLINDER_LG_CALCULATION_HPP
#define CYLINDER_LG_CALCULATION_HPP

/**
 * @file cylinder_lg_calculation.hpp
 * @brief 二维 Eulerian cylinder LG 算例的局部几何、诊断函数和自定义动力学。
 */

#include "cylinder_lg_data.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

namespace SPH
{
//----------------------------------------------------------------------
//	Define geometries and body shapes
//----------------------------------------------------------------------
/** @brief 构造包含 sponge padding 的流体外边界多边形。 */
inline std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge, -DH_sponge));
    water_block_shape.push_back(Vecd(-DL_sponge, DH + DH_sponge));
    water_block_shape.push_back(Vecd(DL, DH + DH_sponge));
    water_block_shape.push_back(Vecd(DL, -DH_sponge));
    water_block_shape.push_back(Vecd(-DL_sponge, -DH_sponge));

    return water_block_shape;
}

/** @brief 构造标记自适应加密区域的矩形多边形。 */
inline std::vector<Vecd> createLocalRefinementRegionShape()
{
    std::vector<Vecd> refinement_region_shape;
    refinement_region_shape.push_back(Vecd(local_refinement_region_x_min, local_refinement_region_y_min));
    refinement_region_shape.push_back(Vecd(local_refinement_region_x_min, local_refinement_region_y_max));
    refinement_region_shape.push_back(Vecd(local_refinement_region_x_max, local_refinement_region_y_max));
    refinement_region_shape.push_back(Vecd(local_refinement_region_x_max, local_refinement_region_y_min));
    refinement_region_shape.push_back(Vecd(local_refinement_region_x_min, local_refinement_region_y_min));

    return refinement_region_shape;
}

/** @brief 创建用于自适应粒子生成的命名 MultiPolygonShape。 */
inline MultiPolygonShape createLocalRefinementRegion()
{
    return MultiPolygonShape(MultiPolygon(createLocalRefinementRegionShape()), "LocalRefinementRegion");
}

/** @brief 判断粒子位置是否位于局部加密矩形内。 */
inline bool isInsideLocalRefinementRegion(const Vecd &position)
{
    return position[0] >= local_refinement_region_x_min &&
           position[0] <= local_refinement_region_x_max &&
           position[1] >= local_refinement_region_y_min &&
           position[1] <= local_refinement_region_y_max;
}

/** @brief 判断粒子位置是否位于入口/出口 sponge buffer 中。 */
inline bool isInsideLocalRefinementOuterBuffer(const Vecd &position)
{
    return position[0] < Real(0.0) || position[1] < Real(0.0) || position[1] > DH;
}

/** @brief 输出加密区内的粒子数、体积和 smoothing-ratio 统计。 */
inline void reportLocalRefinementParticleStats(SPHBody &fluid_body)
{
    BaseParticles &particles = fluid_body.getBaseParticles();
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *h_ratio = particles.getVariableDataByName<Real>("SmoothingLengthRatio");
    const size_t total_particles = particles.TotalRealParticles();

    size_t region_particles = 0;
    Real region_min_vol = std::numeric_limits<Real>::max();
    Real region_max_vol = Real(0.0);
    Real region_min_h_ratio = std::numeric_limits<Real>::max();
    Real region_max_h_ratio = Real(0.0);
    for (size_t i = 0; i != total_particles; ++i)
    {
        if (!isInsideLocalRefinementRegion(pos[i]))
        {
            continue;
        }
        ++region_particles;
        region_min_vol = SMIN(region_min_vol, vol[i]);
        region_max_vol = SMAX(region_max_vol, vol[i]);
        region_min_h_ratio = SMIN(region_min_h_ratio, h_ratio[i]);
        region_max_h_ratio = SMAX(region_max_h_ratio, h_ratio[i]);
    }

    std::cout << "[LocalRefinement] particles total=" << total_particles
              << ", region=" << region_particles
              << ", region_vol_min=" << (region_particles > 0 ? region_min_vol : Real(0.0))
              << ", region_vol_max=" << (region_particles > 0 ? region_max_vol : Real(0.0))
              << ", region_h_ratio_min=" << (region_particles > 0 ? region_min_h_ratio : Real(0.0))
              << ", region_h_ratio_max=" << (region_particles > 0 ? region_max_h_ratio : Real(0.0))
              << std::endl;
}

/** @brief 判断 Real 数值是否为有限值。 */
inline bool isFiniteReal(Real value)
{
    return std::isfinite(static_cast<double>(value));
}

/** @brief 判断 SPHinXsys 向量的每个分量是否为有限值。 */
inline bool isFiniteVector(const Vecd &value)
{
    for (int axis = 0; axis != Dimensions; ++axis)
    {
        if (!isFiniteReal(value[axis]))
        {
            return false;
        }
    }
    return true;
}

/** @brief 向标准输出写出紧凑的向量诊断文本。 */
inline void writeVectorForDiagnostic(const Vecd &value)
{
    std::cout << "(";
    for (int axis = 0; axis != Dimensions; ++axis)
    {
        if (axis != 0)
        {
            std::cout << ", ";
        }
        std::cout << value[axis];
    }
    std::cout << ")";
}

/** @brief 只用有限诊断值更新 min/max 范围。 */
inline void updateFiniteRange(Real value, Real &min_value, Real &max_value, size_t &finite_count)
{
    if (!isFiniteReal(value))
    {
        return;
    }
    min_value = SMIN(min_value, value);
    max_value = SMAX(max_value, value);
    ++finite_count;
}

/** @brief 用于监控局部加密守恒漂移的质量和动量分桶统计。 */
struct LocalRefinementConservationStats
{
    size_t total_count = 0;
    size_t region_count = 0;
    size_t outer_buffer_count = 0;
    size_t main_domain_count = 0;
    size_t coarse_main_count = 0;
    Real total_mass = 0.0;
    Real region_mass = 0.0;
    Real outer_buffer_mass = 0.0;
    Real main_domain_mass = 0.0;
    Real coarse_main_mass = 0.0;
    Vecd total_momentum = Vecd::Zero();
    Vecd region_momentum = Vecd::Zero();
    Vecd outer_buffer_momentum = Vecd::Zero();
    Vecd main_domain_momentum = Vecd::Zero();
    Vecd coarse_main_momentum = Vecd::Zero();
    Real min_rho = std::numeric_limits<Real>::max();
    Real max_rho = -std::numeric_limits<Real>::max();
    Real min_mass = std::numeric_limits<Real>::max();
    Real max_mass = -std::numeric_limits<Real>::max();
    size_t min_rho_index = 0;
    size_t max_rho_index = 0;
    size_t min_mass_index = 0;
    size_t max_mass_index = 0;
};

/** @brief 将单个粒子的守恒量累加到指定诊断分桶。 */
inline void accumulateConservationBucket(
    Real mass, const Vecd &momentum, size_t &count, Real &mass_sum, Vecd &momentum_sum)
{
    ++count;
    mass_sum += mass;
    momentum_sum += momentum;
}

/** @brief 计算全域、加密区、外层 buffer 和非加密主区域的守恒统计。 */
inline LocalRefinementConservationStats computeLocalRefinementConservationStats(SPHBody &fluid_body)
{
    LocalRefinementConservationStats stats;

    BaseParticles &particles = fluid_body.getBaseParticles();
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    Vecd *momentum = particles.getVariableDataByName<Vecd>("Momentum");
    const size_t total_particles = particles.TotalRealParticles();

    for (size_t i = 0; i != total_particles; ++i)
    {
        accumulateConservationBucket(
            mass[i], momentum[i], stats.total_count, stats.total_mass, stats.total_momentum);

        const bool in_region = isInsideLocalRefinementRegion(pos[i]);
        const bool in_outer_buffer = isInsideLocalRefinementOuterBuffer(pos[i]);
        if (in_region)
        {
            accumulateConservationBucket(
                mass[i], momentum[i], stats.region_count, stats.region_mass, stats.region_momentum);
        }
        if (in_outer_buffer)
        {
            accumulateConservationBucket(
                mass[i], momentum[i], stats.outer_buffer_count,
                stats.outer_buffer_mass, stats.outer_buffer_momentum);
        }
        else
        {
            accumulateConservationBucket(
                mass[i], momentum[i], stats.main_domain_count,
                stats.main_domain_mass, stats.main_domain_momentum);
            if (!in_region)
            {
                accumulateConservationBucket(
                    mass[i], momentum[i], stats.coarse_main_count,
                    stats.coarse_main_mass, stats.coarse_main_momentum);
            }
        }

        if (rho[i] < stats.min_rho)
        {
            stats.min_rho = rho[i];
            stats.min_rho_index = i;
        }
        if (rho[i] > stats.max_rho)
        {
            stats.max_rho = rho[i];
            stats.max_rho_index = i;
        }
        if (mass[i] < stats.min_mass)
        {
            stats.min_mass = mass[i];
            stats.min_mass_index = i;
        }
        if (mass[i] > stats.max_mass)
        {
            stats.max_mass = mass[i];
            stats.max_mass_index = i;
        }
    }

    if (total_particles == 0)
    {
        stats.min_rho = Real(0.0);
        stats.max_rho = Real(0.0);
        stats.min_mass = Real(0.0);
        stats.max_mass = Real(0.0);
    }
    return stats;
}

/** @brief 计算相对参考分桶的带符号质量漂移。 */
inline Real relativeMassDrift(Real value, Real reference)
{
    const Real scale = SMAX(std::fabs(reference), TinyReal);
    return (value - reference) / scale;
}

/** @brief 使用动量尺度或质量-速度尺度计算归一化动量漂移。 */
inline Real relativeMomentumDrift(const Vecd &value, const Vecd &reference, Real reference_mass)
{
    const Real scale = SMAX(SMAX(reference.norm(), std::fabs(reference_mass) * U_f), TinyReal);
    return (value - reference).norm() / scale;
}

/** @brief 输出一个守恒分桶及其相对参考状态的漂移。 */
inline void writeConservationBucket(
    const std::string &name, size_t count, Real mass, const Vecd &momentum,
    size_t reference_count, Real reference_mass, const Vecd &reference_momentum)
{
    std::cout << ", " << name << "_count=" << count
              << ", " << name << "_mass=" << mass
              << ", " << name << "_mass_rel_drift=" << relativeMassDrift(mass, reference_mass)
              << ", " << name << "_momentum=";
    writeVectorForDiagnostic(momentum);
    std::cout << ", " << name << "_momentum_rel_drift="
              << relativeMomentumDrift(momentum, reference_momentum, reference_mass);
    if (count != reference_count)
    {
        std::cout << ", " << name << "_count_initial=" << reference_count;
    }
}

/** @brief 在指定仿真阶段输出守恒分桶以及 rho/mass 范围。 */
inline void reportLocalRefinementConservationStats(
    SPHBody &fluid_body, const std::string &stage,
    const LocalRefinementConservationStats &reference)
{
    const LocalRefinementConservationStats stats =
        computeLocalRefinementConservationStats(fluid_body);

    BaseParticles &particles = fluid_body.getBaseParticles();
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");

    std::cout << "[LocalRefinement][Conservation] " << stage;
    writeConservationBucket(
        "total", stats.total_count, stats.total_mass, stats.total_momentum,
        reference.total_count, reference.total_mass, reference.total_momentum);
    writeConservationBucket(
        "region", stats.region_count, stats.region_mass, stats.region_momentum,
        reference.region_count, reference.region_mass, reference.region_momentum);
    writeConservationBucket(
        "outer_buffer", stats.outer_buffer_count, stats.outer_buffer_mass, stats.outer_buffer_momentum,
        reference.outer_buffer_count, reference.outer_buffer_mass, reference.outer_buffer_momentum);
    writeConservationBucket(
        "main_domain", stats.main_domain_count, stats.main_domain_mass, stats.main_domain_momentum,
        reference.main_domain_count, reference.main_domain_mass, reference.main_domain_momentum);
    writeConservationBucket(
        "coarse_main", stats.coarse_main_count, stats.coarse_main_mass, stats.coarse_main_momentum,
        reference.coarse_main_count, reference.coarse_main_mass, reference.coarse_main_momentum);
    std::cout << std::endl;

    std::cout << "[LocalRefinement][ConservationRange] " << stage
              << ": rho=[" << stats.min_rho << "@" << stats.min_rho_index << ", "
              << stats.max_rho << "@" << stats.max_rho_index << "]"
              << ", mass=[" << stats.min_mass << "@" << stats.min_mass_index << ", "
              << stats.max_mass << "@" << stats.max_mass_index << "]"
              << ", min_rho_pos=";
    writeVectorForDiagnostic(pos[stats.min_rho_index]);
    std::cout << ", min_mass_pos=";
    writeVectorForDiagnostic(pos[stats.min_mass_index]);
    std::cout << std::endl;
}

/** @brief 记录密度松弛中最负通量贡献的详细数据。 */
struct MassFluxContributionDiagnostic
{
    bool found = false;
    size_t body_index = 0;
    size_t neighbor_position = 0;
    size_t neighbor_index = 0;
    Real contribution = 0.0;
    Real r_ij = 0.0;
    Real dW_ij = 0.0;
    Real dW_ijV_j = 0.0;
    Real Vol_j = 0.0;
    Real rho_j = 0.0;
    Real p_j = 0.0;
    Real interface_rho = 0.0;
    Real interface_p = 0.0;
    Vecd e_ij = Vecd::Zero();
    Vecd vel_j = Vecd::Zero();
    Vecd interface_vel = Vecd::Zero();
    Vecd wall_normal = Vecd::Zero();
    Vecd wall_vel_ave = Vecd::Zero();
};

/** @brief 将密度松弛质量通量拆分为 inner 和 wall-contact 贡献。 */
struct MassFluxDiagnosticRecord
{
    size_t inner_neighbor_count = 0;
    size_t wall_neighbor_count = 0;
    Real inner_sum = 0.0;
    Real wall_sum = 0.0;
    Real total_sum = 0.0;
    MassFluxContributionDiagnostic most_negative_inner;
    MassFluxContributionDiagnostic most_negative_wall;
};

/** @brief 保留当前观察到的最负质量通量贡献。 */
inline void updateMostNegativeMassFluxContribution(
    Real contribution, size_t body_index, size_t neighbor_position, size_t neighbor_index,
    Real r_ij, Real dW_ij, Real dW_ijV_j, Real Vol_j, Real rho_j, Real p_j,
    const Vecd &e_ij, const Vecd &vel_j, const FluidStateOut &interface_state,
    MassFluxContributionDiagnostic &most_negative)
{
    if (most_negative.found && contribution >= most_negative.contribution)
    {
        return;
    }

    most_negative.found = true;
    most_negative.body_index = body_index;
    most_negative.neighbor_position = neighbor_position;
    most_negative.neighbor_index = neighbor_index;
    most_negative.contribution = contribution;
    most_negative.r_ij = r_ij;
    most_negative.dW_ij = dW_ij;
    most_negative.dW_ijV_j = dW_ijV_j;
    most_negative.Vol_j = Vol_j;
    most_negative.rho_j = rho_j;
    most_negative.p_j = p_j;
    most_negative.e_ij = e_ij;
    most_negative.vel_j = vel_j;
    most_negative.interface_rho = interface_state.rho_;
    most_negative.interface_p = interface_state.p_;
    most_negative.interface_vel = interface_state.vel_;
}

/** @brief 输出一个已保存的质量通量贡献及其几何和界面状态细节。 */
inline void writeMassFluxContributionDiagnostic(
    const MassFluxContributionDiagnostic &contribution, const std::string &label)
{
    std::cout << "[LocalRefinement][DensityFluxDiagnostic] " << label
              << ": found=" << (contribution.found ? "true" : "false");
    if (!contribution.found)
    {
        std::cout << std::endl;
        return;
    }

    std::cout << ", body=" << contribution.body_index
              << ", n=" << contribution.neighbor_position
              << ", j=" << contribution.neighbor_index
              << ", contribution=" << contribution.contribution
              << ", r=" << contribution.r_ij
              << ", dW=" << contribution.dW_ij
              << ", dW_Vj=" << contribution.dW_ijV_j
              << ", Vol_j=" << contribution.Vol_j
              << ", rho_j=" << contribution.rho_j
              << ", p_j=" << contribution.p_j
              << ", e=";
    writeVectorForDiagnostic(contribution.e_ij);
    std::cout << ", vel_j=";
    writeVectorForDiagnostic(contribution.vel_j);
    std::cout << ", interface_rho=" << contribution.interface_rho
              << ", interface_p=" << contribution.interface_p
              << ", interface_vel=";
    writeVectorForDiagnostic(contribution.interface_vel);
    if (contribution.wall_normal.norm() > Real(0.0) || contribution.wall_vel_ave.norm() > Real(0.0))
    {
        std::cout << ", wall_normal=";
        writeVectorForDiagnostic(contribution.wall_normal);
        std::cout << ", wall_vel_ave=";
        writeVectorForDiagnostic(contribution.wall_vel_ave);
    }
    std::cout << std::endl;
}

/** @brief 为单个粒子重新计算密度松弛通量项，用于定位问题。 */
inline MassFluxDiagnosticRecord computeDensityRelaxationMassFluxDiagnostic(
    SPHBody &fluid_body, BaseInnerRelation &inner_relation,
    BaseContactRelation &wall_contact_relation, size_t particle_index)
{
    MassFluxDiagnosticRecord record;

    BaseParticles &particles = fluid_body.getBaseParticles();
    Fluid &fluid = DynamicCast<Fluid>(&fluid_body, particles.getBaseMaterial());
    AcousticRiemannSolver riemann_solver(fluid, fluid, Real(15.0));

    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *pressure = particles.getVariableDataByName<Real>("Pressure");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Vecd *velocity = particles.getVariableDataByName<Vecd>("Velocity");
    FluidStateIn state_i(rho[particle_index], velocity[particle_index], pressure[particle_index]);

    const Neighborhood &inner_neighborhood = inner_relation.inner_configuration_[particle_index];
    record.inner_neighbor_count = inner_neighborhood.current_size_;
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        const size_t index_j = inner_neighborhood.j_[n];
        const Vecd &e_ij = inner_neighborhood.e_ij_[n];
        const Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * vol[index_j];
        FluidStateIn state_j(rho[index_j], velocity[index_j], pressure[index_j]);
        const FluidStateOut interface_state = riemann_solver.InterfaceState(state_i, state_j, e_ij);
        const Real contribution =
            -2.0 * vol[particle_index] *
            (interface_state.rho_ * interface_state.vel_).dot(e_ij) * dW_ijV_j;

        record.inner_sum += contribution;
        updateMostNegativeMassFluxContribution(
            contribution, 0, n, index_j, inner_neighborhood.r_ij_[n],
            inner_neighborhood.dW_ij_[n], dW_ijV_j, vol[index_j], rho[index_j],
            pressure[index_j], e_ij, velocity[index_j], interface_state,
            record.most_negative_inner);
    }

    for (size_t k = 0; k != wall_contact_relation.contact_configuration_.size(); ++k)
    {
        BaseParticles *wall_particles = wall_contact_relation.contact_particles_[k];
        Solid &solid_material = DynamicCast<Solid>(&fluid_body, wall_particles->getBaseMaterial());
        Real *wall_vol = wall_particles->getVariableDataByName<Real>("VolumetricMeasure");
        Vecd *wall_normal = wall_particles->getVariableDataByName<Vecd>("NormalDirection");
        Vecd *wall_velocity_average = solid_material.AverageVelocity(wall_particles);
        const Neighborhood &wall_neighborhood =
            wall_contact_relation.contact_configuration_[k][particle_index];
        record.wall_neighbor_count += wall_neighborhood.current_size_;
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            const size_t index_j = wall_neighborhood.j_[n];
            const Vecd &e_ij = wall_neighborhood.e_ij_[n];
            const Real dW_ijV_j = wall_neighborhood.dW_ij_[n] * wall_vol[index_j];
            Vecd vel_j_in_wall = 2.0 * wall_velocity_average[index_j] - state_i.vel_;
            Real rho_j_in_wall = state_i.rho_;
            Real p_j_in_wall = state_i.p_;
            FluidStateIn state_j(rho_j_in_wall, vel_j_in_wall, p_j_in_wall);
            const FluidStateOut interface_state =
                riemann_solver.InterfaceState(state_i, state_j, wall_normal[index_j]);
            const Real contribution =
                -2.0 * vol[particle_index] *
                (interface_state.rho_ * interface_state.vel_).dot(e_ij) * dW_ijV_j;

            record.wall_sum += contribution;
            updateMostNegativeMassFluxContribution(
                contribution, k, n, index_j, wall_neighborhood.r_ij_[n],
                wall_neighborhood.dW_ij_[n], dW_ijV_j, wall_vol[index_j],
                state_i.rho_, state_i.p_, e_ij, vel_j_in_wall, interface_state,
                record.most_negative_wall);
            record.most_negative_wall.wall_normal = wall_normal[index_j];
            record.most_negative_wall.wall_vel_ave = wall_velocity_average[index_j];
        }
    }

    record.total_sum = record.inner_sum + record.wall_sum;
    return record;
}

/** @brief 输出单个粒子在指定阶段的密度松弛诊断。 */
inline void reportDensityRelaxationMassFluxDiagnostic(
    SPHBody &fluid_body, BaseInnerRelation &inner_relation,
    BaseContactRelation &wall_contact_relation, size_t particle_index,
    Real dt, const std::string &stage)
{
    BaseParticles &particles = fluid_body.getBaseParticles();
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *h_ratio = particles.getVariableDataByName<Real>("SmoothingLengthRatio");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    Real *pressure = particles.getVariableDataByName<Real>("Pressure");
    Vecd *velocity = particles.getVariableDataByName<Vecd>("Velocity");
    Vecd *momentum = particles.getVariableDataByName<Vecd>("Momentum");
    Real *mass_change_rate = particles.getVariableDataByName<Real>("MassChangeRate");
    int *indicator = particles.getVariableDataByName<int>("Indicator");
    int *smeared_surface = particles.getVariableDataByName<int>("SmearedSurface");
    Vecd *normal = particles.getVariableDataByName<Vecd>("NormalDirection");

    const MassFluxDiagnosticRecord record =
        computeDensityRelaxationMassFluxDiagnostic(
            fluid_body, inner_relation, wall_contact_relation, particle_index);
    const Real previous_mass_estimate = mass[particle_index] - mass_change_rate[particle_index] * dt;
    const Real predicted_mass_from_split = previous_mass_estimate + record.total_sum * dt;
    const Real stored_minus_split = mass_change_rate[particle_index] - record.total_sum;

    std::cout << "[LocalRefinement][DensityFluxDiagnostic] " << stage
              << ": i=" << particle_index
              << ", in_region=" << (isInsideLocalRefinementRegion(pos[particle_index]) ? "true" : "false")
              << ", pos=";
    writeVectorForDiagnostic(pos[particle_index]);
    std::cout << ", Vol=" << vol[particle_index]
              << ", h_ratio=" << h_ratio[particle_index]
              << ", rho=" << rho[particle_index]
              << ", mass=" << mass[particle_index]
              << ", p=" << pressure[particle_index]
              << ", vel=";
    writeVectorForDiagnostic(velocity[particle_index]);
    std::cout << ", mom=";
    writeVectorForDiagnostic(momentum[particle_index]);
    std::cout << ", dt=" << dt
              << ", MassChangeRate_stored=" << mass_change_rate[particle_index]
              << ", previous_mass_estimate=" << previous_mass_estimate
              << ", predicted_mass_from_split=" << predicted_mass_from_split
              << ", stored_minus_split=" << stored_minus_split
              << ", indicator=" << indicator[particle_index]
              << ", smeared_surface=" << smeared_surface[particle_index]
              << ", normal=";
    writeVectorForDiagnostic(normal[particle_index]);
    std::cout << std::endl;

    std::cout << "[LocalRefinement][DensityFluxDiagnostic] " << stage
              << ": inner_neighbors=" << record.inner_neighbor_count
              << ", wall_neighbors=" << record.wall_neighbor_count
              << ", sum_inner=" << record.inner_sum
              << ", sum_wall=" << record.wall_sum
              << ", sum_total=" << record.total_sum
              << ", inner_fraction="
              << (std::fabs(record.total_sum) > TinyReal ? record.inner_sum / record.total_sum : Real(0.0))
              << ", wall_fraction="
              << (std::fabs(record.total_sum) > TinyReal ? record.wall_sum / record.total_sum : Real(0.0))
              << std::endl;
    writeMassFluxContributionDiagnostic(record.most_negative_inner, stage + " most_negative_inner");
    writeMassFluxContributionDiagnostic(record.most_negative_wall, stage + " most_negative_wall");
}

/** @brief 扫描流体状态变量中的非有限值，并输出首个异常粒子。 */
inline bool reportFluidFiniteState(SPHBody &fluid_body, const std::string &stage,
                                   size_t *first_bad_particle = nullptr)
{
    BaseParticles &particles = fluid_body.getBaseParticles();
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *h_ratio = enable_local_refinement
                        ? particles.getVariableDataByName<Real>("SmoothingLengthRatio")
                        : nullptr;
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    Real *pressure = particles.getVariableDataByName<Real>("Pressure");
    Vecd *velocity = particles.getVariableDataByName<Vecd>("Velocity");
    Vecd *momentum = particles.getVariableDataByName<Vecd>("Momentum");
    Vecd *momentum_change_rate = particles.getVariableDataByName<Vecd>("MomentumChangeRate");
    Vecd *force = particles.getVariableDataByName<Vecd>("Force");
    Vecd *force_prior = particles.getVariableDataByName<Vecd>("ForcePrior");
    Vecd *viscous_force = particles.getVariableDataByName<Vecd>("ViscousForce");
    Vecd *previous_viscous_force = particles.getVariableDataByName<Vecd>("PreviousViscousForce");
    const size_t total_particles = particles.TotalRealParticles();

    Real rho_min = std::numeric_limits<Real>::max();
    Real rho_max = -std::numeric_limits<Real>::max();
    Real pressure_min = std::numeric_limits<Real>::max();
    Real pressure_max = -std::numeric_limits<Real>::max();
    Real mass_min = std::numeric_limits<Real>::max();
    Real mass_max = -std::numeric_limits<Real>::max();
    Real speed_max = Real(0.0);
    Real momentum_norm_max = Real(0.0);
    Real force_norm_max = Real(0.0);
    Real force_prior_norm_max = Real(0.0);
    Real viscous_force_norm_max = Real(0.0);
    Real previous_viscous_force_norm_max = Real(0.0);
    size_t first_bad = total_particles;
    size_t rho_count = 0;
    size_t pressure_count = 0;
    size_t mass_count = 0;
    for (size_t i = 0; i != total_particles; ++i)
    {
        const bool h_ratio_finite = h_ratio == nullptr || isFiniteReal(h_ratio[i]);
        const bool finite =
            isFiniteVector(pos[i]) && isFiniteReal(vol[i]) && h_ratio_finite &&
            isFiniteReal(rho[i]) && isFiniteReal(mass[i]) && isFiniteReal(pressure[i]) &&
            isFiniteVector(velocity[i]) && isFiniteVector(momentum[i]) &&
            isFiniteVector(momentum_change_rate[i]) && isFiniteVector(force[i]) &&
            isFiniteVector(force_prior[i]) && isFiniteVector(viscous_force[i]) &&
            isFiniteVector(previous_viscous_force[i]);
        if (!finite && first_bad == total_particles)
        {
            first_bad = i;
        }
        updateFiniteRange(rho[i], rho_min, rho_max, rho_count);
        updateFiniteRange(pressure[i], pressure_min, pressure_max, pressure_count);
        updateFiniteRange(mass[i], mass_min, mass_max, mass_count);
        if (isFiniteVector(velocity[i]))
        {
            speed_max = SMAX(speed_max, velocity[i].norm());
        }
        if (isFiniteVector(momentum[i]))
        {
            momentum_norm_max = SMAX(momentum_norm_max, momentum[i].norm());
        }
        if (isFiniteVector(force[i]))
        {
            force_norm_max = SMAX(force_norm_max, force[i].norm());
        }
        if (isFiniteVector(force_prior[i]))
        {
            force_prior_norm_max = SMAX(force_prior_norm_max, force_prior[i].norm());
        }
        if (isFiniteVector(viscous_force[i]))
        {
            viscous_force_norm_max = SMAX(viscous_force_norm_max, viscous_force[i].norm());
        }
        if (isFiniteVector(previous_viscous_force[i]))
        {
            previous_viscous_force_norm_max = SMAX(previous_viscous_force_norm_max, previous_viscous_force[i].norm());
        }
    }

    if (first_bad != total_particles)
    {
        std::cout << "[LocalRefinement][State] " << stage
                  << ": rho=[" << (rho_count > 0 ? rho_min : Real(0.0)) << ", " << (rho_count > 0 ? rho_max : Real(0.0)) << "]"
                  << ", p=[" << (pressure_count > 0 ? pressure_min : Real(0.0)) << ", " << (pressure_count > 0 ? pressure_max : Real(0.0)) << "]"
                  << ", mass=[" << (mass_count > 0 ? mass_min : Real(0.0)) << ", " << (mass_count > 0 ? mass_max : Real(0.0)) << "]"
                  << ", max|u|=" << speed_max
                  << ", max|mom|=" << momentum_norm_max
                  << ", max|force|=" << force_norm_max
                  << ", max|force_prior|=" << force_prior_norm_max
                  << ", max|viscous_force|=" << viscous_force_norm_max
                  << ", max|previous_viscous_force|=" << previous_viscous_force_norm_max
                  << ", first_bad_particle=" << first_bad
                  << std::endl;
        std::cout << "[LocalRefinement][BadParticle] " << stage
                  << ": i=" << first_bad
                  << ", in_region=" << (isInsideLocalRefinementRegion(pos[first_bad]) ? "true" : "false")
                  << ", pos=";
        writeVectorForDiagnostic(pos[first_bad]);
        std::cout << ", Vol=" << vol[first_bad];
        if (h_ratio != nullptr)
        {
            std::cout << ", h_ratio=" << h_ratio[first_bad];
        }
        std::cout << ", rho=" << rho[first_bad]
                  << ", mass=" << mass[first_bad]
                  << ", p=" << pressure[first_bad]
                  << ", vel=";
        writeVectorForDiagnostic(velocity[first_bad]);
        std::cout << ", mom=";
        writeVectorForDiagnostic(momentum[first_bad]);
        std::cout << ", dmom_dt=";
        writeVectorForDiagnostic(momentum_change_rate[first_bad]);
        std::cout << ", force=";
        writeVectorForDiagnostic(force[first_bad]);
        std::cout << ", force_prior=";
        writeVectorForDiagnostic(force_prior[first_bad]);
        std::cout << ", viscous_force=";
        writeVectorForDiagnostic(viscous_force[first_bad]);
        std::cout << ", previous_viscous_force=";
        writeVectorForDiagnostic(previous_viscous_force[first_bad]);
        std::cout << std::endl;
        std::cout << "[LocalRefinement][BadParticleFinite] " << stage
                  << ": pos=" << isFiniteVector(pos[first_bad])
                  << ", Vol=" << isFiniteReal(vol[first_bad])
                  << ", h_ratio=" << (h_ratio == nullptr || isFiniteReal(h_ratio[first_bad]))
                  << ", rho=" << isFiniteReal(rho[first_bad])
                  << ", mass=" << isFiniteReal(mass[first_bad])
                  << ", p=" << isFiniteReal(pressure[first_bad])
                  << ", vel=" << isFiniteVector(velocity[first_bad])
                  << ", mom=" << isFiniteVector(momentum[first_bad])
                  << ", dmom_dt=" << isFiniteVector(momentum_change_rate[first_bad])
                  << ", force=" << isFiniteVector(force[first_bad])
                  << ", force_prior=" << isFiniteVector(force_prior[first_bad])
                  << ", viscous_force=" << isFiniteVector(viscous_force[first_bad])
                  << ", previous_viscous_force=" << isFiniteVector(previous_viscous_force[first_bad])
                  << std::endl;
    }
    if (first_bad_particle != nullptr)
    {
        *first_bad_particle = first_bad;
    }
    return first_bad == total_particles;
}

/** @brief 扫描密度和质量中的非正值，并输出首个失败粒子。 */
inline bool reportFluidPositiveState(SPHBody &fluid_body, const std::string &stage,
                                     size_t *first_non_positive_particle = nullptr,
                                     bool ignore_local_refinement_outer_buffer = false)
{
    BaseParticles &particles = fluid_body.getBaseParticles();
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    const size_t total_particles = particles.TotalRealParticles();

    for (size_t i = 0; i != total_particles; ++i)
    {
        if (ignore_local_refinement_outer_buffer && isInsideLocalRefinementOuterBuffer(pos[i]))
        {
            continue;
        }
        if (rho[i] > Real(0.0) && mass[i] > Real(0.0))
        {
            continue;
        }

        std::cout << "[LocalRefinement][PositiveState] " << stage
                  << ": first_non_positive i=" << i
                  << ", in_region=" << (isInsideLocalRefinementRegion(pos[i]) ? "true" : "false")
                  << ", pos=";
        writeVectorForDiagnostic(pos[i]);
        std::cout << ", Vol=" << vol[i]
                  << ", rho=" << rho[i]
                  << ", mass=" << mass[i]
                  << std::endl;
        if (first_non_positive_particle != nullptr)
        {
            *first_non_positive_particle = i;
        }
        return false;
    }
    if (first_non_positive_particle != nullptr)
    {
        *first_non_positive_particle = total_particles;
    }
    return true;
}

/** @brief 单个非反射远场边界粒子的诊断快照。 */
struct FarFieldBoundaryDiagnosticRecord
{
    size_t index = 0;
    size_t inner_neighbor_count = 0;
    bool active_boundary = false;
    bool inflow = false;
    bool subsonic = false;
    bool predicted = false;
    Real velocity_boundary_normal = 0.0;
    Real velocity_farfield_normal = 0.0;
    Real inner_weight_summation = 0.0;
    Real bounded_inner_weight_summation = 0.0;
    Real rho_average = 0.0;
    Real vel_normal_average = 0.0;
    Real rho_predicted = 0.0;
    Real rho_predicted_bounded = 0.0;
    Real neighbor_rho_min = std::numeric_limits<Real>::max();
    Real neighbor_rho_max = -std::numeric_limits<Real>::max();
    Real neighbor_vol_min = std::numeric_limits<Real>::max();
    Real neighbor_vol_max = Real(0.0);
    Vecd vel_average = Vecd::Zero();
    Vecd vel_tangential_average = Vecd::Zero();
};

/** @brief 为远场诊断累积邻居 rho/volume 范围。 */
inline void updateFarFieldNeighborRange(Real rho, Real vol, FarFieldBoundaryDiagnosticRecord &record)
{
    record.neighbor_rho_min = SMIN(record.neighbor_rho_min, rho);
    record.neighbor_rho_max = SMAX(record.neighbor_rho_max, rho);
    record.neighbor_vol_min = SMIN(record.neighbor_vol_min, vol);
    record.neighbor_vol_max = SMAX(record.neighbor_vol_max, vol);
}

/** @brief 重建单个边界粒子的远场分支预测。 */
inline FarFieldBoundaryDiagnosticRecord computeFarFieldBoundaryDiagnostic(
    SPHBody &fluid_body, BaseInnerRelation &inner_relation, size_t particle_index)
{
    FarFieldBoundaryDiagnosticRecord record;
    record.index = particle_index;

    BaseParticles &particles = fluid_body.getBaseParticles();
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Vecd *velocity = particles.getVariableDataByName<Vecd>("Velocity");
    int *indicator = particles.getVariableDataByName<int>("Indicator");
    int *smeared_surface = particles.getVariableDataByName<int>("SmearedSurface");
    Vecd *normal = particles.getVariableDataByName<Vecd>("NormalDirection");

    record.active_boundary =
        indicator[particle_index] == 1 || smeared_surface[particle_index] == 1;
    if (!record.active_boundary)
    {
        return record;
    }

    const Vecd farfield_velocity(U_f, 0.0);
    record.velocity_farfield_normal = farfield_velocity.dot(normal[particle_index]);
    record.velocity_boundary_normal = velocity[particle_index].dot(normal[particle_index]);
    record.inflow =
        normal[particle_index][0] <= 0.0 ||
        std::fabs(normal[particle_index][1]) > std::fabs(normal[particle_index][0]);
    record.subsonic = std::fabs(record.velocity_boundary_normal) < c_f;

    Real rho_summation = 0.0;
    Real vel_normal_summation = 0.0;
    Vecd vel_summation = Vecd::Zero();
    Vecd vel_tangential_summation = Vecd::Zero();
    const Neighborhood &inner_neighborhood = inner_relation.inner_configuration_[particle_index];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        const size_t index_j = inner_neighborhood.j_[n];
        if (indicator[index_j] == 1)
        {
            continue;
        }

        ++record.inner_neighbor_count;
        rho_summation += rho[index_j];
        vel_summation += velocity[index_j];
        vel_normal_summation += velocity[index_j].dot(normal[particle_index]);
        vel_tangential_summation +=
            velocity[index_j] - velocity[index_j].dot(normal[particle_index]) * normal[particle_index];
        record.inner_weight_summation += inner_neighborhood.W_ij_[n] * vol[index_j];
        updateFarFieldNeighborRange(rho[index_j], vol[index_j], record);
    }

    const Real neighbor_count = static_cast<Real>(record.inner_neighbor_count);
    record.rho_average = rho_summation / (neighbor_count + TinyReal);
    record.vel_normal_average = vel_normal_summation / (neighbor_count + TinyReal);
    record.vel_average = vel_summation / (neighbor_count + TinyReal);
    record.vel_tangential_average = vel_tangential_summation / (neighbor_count + TinyReal);
    record.bounded_inner_weight_summation =
        SMIN(Real(1.0), SMAX(Real(0.0), record.inner_weight_summation));

    record.predicted = true;
    if (record.inflow)
    {
        record.rho_predicted = record.subsonic
                                   ? record.rho_average * record.inner_weight_summation +
                                         rho0_f * (1.0 - record.inner_weight_summation)
                                   : rho0_f;
    }
    else
    {
        record.rho_predicted = record.subsonic
                                   ? record.rho_average * record.inner_weight_summation +
                                         rho0_f * (1.0 - record.inner_weight_summation)
                                   : record.rho_average + TinyReal;
    }
    record.rho_predicted_bounded = record.subsonic
                                       ? record.rho_average * record.bounded_inner_weight_summation +
                                             rho0_f * (1.0 - record.bounded_inner_weight_summation)
                                       : record.rho_predicted;
    return record;
}

/** @brief 输出单个粒子的远场诊断快照。 */
inline void writeFarFieldBoundaryDiagnostic(
    SPHBody &fluid_body, BaseInnerRelation &inner_relation, size_t particle_index,
    const std::string &stage)
{
    BaseParticles &particles = fluid_body.getBaseParticles();
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    Real *pressure = particles.getVariableDataByName<Real>("Pressure");
    Vecd *velocity = particles.getVariableDataByName<Vecd>("Velocity");
    Vecd *momentum = particles.getVariableDataByName<Vecd>("Momentum");
    int *indicator = particles.getVariableDataByName<int>("Indicator");
    int *smeared_surface = particles.getVariableDataByName<int>("SmearedSurface");
    Vecd *normal = particles.getVariableDataByName<Vecd>("NormalDirection");

    const FarFieldBoundaryDiagnosticRecord record =
        computeFarFieldBoundaryDiagnostic(fluid_body, inner_relation, particle_index);
    const std::string direction = record.inflow ? "inflow" : "outflow";
    const std::string acoustic = record.subsonic ? "subsonic" : "supersonic";

    std::cout << "[LocalRefinement][FarFieldDiagnostic] " << stage
              << ": i=" << particle_index
              << ", active_boundary=" << (record.active_boundary ? "true" : "false")
              << ", indicator=" << indicator[particle_index]
              << ", smeared_surface=" << smeared_surface[particle_index]
              << ", in_region=" << (isInsideLocalRefinementRegion(pos[particle_index]) ? "true" : "false")
              << ", branch=" << direction << "/" << acoustic
              << ", pos=";
    writeVectorForDiagnostic(pos[particle_index]);
    std::cout << ", n=";
    writeVectorForDiagnostic(normal[particle_index]);
    std::cout << ", u_n=" << record.velocity_boundary_normal
              << ", u_farfield_n=" << record.velocity_farfield_normal
              << ", inner_neighbors=" << record.inner_neighbor_count
              << ", w=" << record.inner_weight_summation
              << ", w_bounded=" << record.bounded_inner_weight_summation
              << ", rho_average=" << record.rho_average
              << ", rho_predicted=" << record.rho_predicted
              << ", rho_predicted_bounded=" << record.rho_predicted_bounded
              << ", rho_current=" << rho[particle_index]
              << ", mass_current=" << mass[particle_index]
              << ", p_current=" << pressure[particle_index]
              << ", Vol=" << vol[particle_index]
              << ", vel=";
    writeVectorForDiagnostic(velocity[particle_index]);
    std::cout << ", mom=";
    writeVectorForDiagnostic(momentum[particle_index]);
    std::cout << std::endl;

    std::cout << "[LocalRefinement][FarFieldDiagnostic] " << stage
              << ": neighbor_rho=["
              << (record.inner_neighbor_count > 0 ? record.neighbor_rho_min : Real(0.0))
              << ", "
              << (record.inner_neighbor_count > 0 ? record.neighbor_rho_max : Real(0.0))
              << "], neighbor_vol=["
              << (record.inner_neighbor_count > 0 ? record.neighbor_vol_min : Real(0.0))
              << ", "
              << (record.inner_neighbor_count > 0 ? record.neighbor_vol_max : Real(0.0))
              << "], vel_normal_average=" << record.vel_normal_average
              << ", vel_average=";
    writeVectorForDiagnostic(record.vel_average);
    std::cout << ", vel_tangential_average=";
    writeVectorForDiagnostic(record.vel_tangential_average);
    std::cout << std::endl;
}

/** @brief 扫描远场粒子，并在预测密度存在风险时展开诊断。 */
inline void reportFarFieldBoundaryRiskScan(
    SPHBody &fluid_body, BaseInnerRelation &inner_relation, const std::string &stage)
{
    BaseParticles &particles = fluid_body.getBaseParticles();
    const size_t total_particles = particles.TotalRealParticles();

    size_t active_boundary_count = 0;
    size_t min_predicted_index = total_particles;
    size_t max_weight_index = total_particles;
    Real min_raw_predicted_rho = std::numeric_limits<Real>::max();
    Real min_bounded_predicted_rho = std::numeric_limits<Real>::max();
    Real max_inner_weight = -std::numeric_limits<Real>::max();

    for (size_t i = 0; i != total_particles; ++i)
    {
        const FarFieldBoundaryDiagnosticRecord record =
            computeFarFieldBoundaryDiagnostic(fluid_body, inner_relation, i);
        if (!record.active_boundary || !record.predicted)
        {
            continue;
        }
        ++active_boundary_count;
        if (record.rho_predicted < min_raw_predicted_rho)
        {
            min_raw_predicted_rho = record.rho_predicted;
            min_predicted_index = i;
        }
        min_bounded_predicted_rho =
            SMIN(min_bounded_predicted_rho, record.rho_predicted_bounded);
        if (record.inner_weight_summation > max_inner_weight)
        {
            max_inner_weight = record.inner_weight_summation;
            max_weight_index = i;
        }
    }

    const bool has_active_boundary = active_boundary_count > 0;
    const bool risky_prediction =
        has_active_boundary &&
        (min_raw_predicted_rho <= Real(0.0));
    if (!risky_prediction)
    {
        return;
    }

    std::cout << "[LocalRefinement][FarFieldScan] " << stage
              << ": active_boundary=" << active_boundary_count
              << ", min_raw_predicted_rho=" << min_raw_predicted_rho
              << ", min_bounded_predicted_rho=" << min_bounded_predicted_rho
              << ", min_predicted_i=" << min_predicted_index
              << ", max_w=" << max_inner_weight
              << ", max_w_i=" << max_weight_index << std::endl;
    if (min_predicted_index != total_particles)
    {
        writeFarFieldBoundaryDiagnostic(
            fluid_body, inner_relation, min_predicted_index, stage + " min_predicted");
    }
    if (max_weight_index != total_particles && max_weight_index != min_predicted_index)
    {
        writeFarFieldBoundaryDiagnostic(
            fluid_body, inner_relation, max_weight_index, stage + " max_weight");
    }
}

/** @brief 输出单个粒子周围 inner 和 wall-contact 邻居质量，用于黏性力诊断。 */
inline void reportViscousForceNeighborhoodDiagnostics(
    SPHBody &fluid_body, BaseInnerRelation &inner_relation,
    BaseContactRelation &wall_contact_relation, size_t particle_index,
    const std::string &stage)
{
    BaseParticles &particles = fluid_body.getBaseParticles();
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Vecd *velocity = particles.getVariableDataByName<Vecd>("Velocity");
    Real *h_ratio = particles.getVariableDataByName<Real>("SmoothingLengthRatio");

    const Neighborhood &inner_neighborhood = inner_relation.inner_configuration_[particle_index];
    std::cout << "[LocalRefinement][ViscousNeighborhood] " << stage
              << ": i=" << particle_index
              << ", inner_neighbors=" << inner_neighborhood.current_size_;

    Real min_inner_r = std::numeric_limits<Real>::max();
    size_t first_bad_inner = inner_neighborhood.current_size_;
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        const size_t j = inner_neighborhood.j_[n];
        min_inner_r = SMIN(min_inner_r, inner_neighborhood.r_ij_[n]);
        const bool finite_pair =
            isFiniteReal(inner_neighborhood.r_ij_[n]) &&
            isFiniteReal(inner_neighborhood.dW_ij_[n]) &&
            isFiniteVector(inner_neighborhood.e_ij_[n]) &&
            isFiniteReal(vol[j]) && isFiniteVector(velocity[j]);
        if (!finite_pair && first_bad_inner == inner_neighborhood.current_size_)
        {
            first_bad_inner = n;
        }
    }
    std::cout << ", min_inner_r="
              << (inner_neighborhood.current_size_ > 0 ? min_inner_r : Real(0.0));
    if (first_bad_inner != inner_neighborhood.current_size_)
    {
        const size_t j = inner_neighborhood.j_[first_bad_inner];
        std::cout << ", first_bad_inner_n=" << first_bad_inner
                  << ", j=" << j
                  << ", r=" << inner_neighborhood.r_ij_[first_bad_inner]
                  << ", current_distance=" << (pos[particle_index] - pos[j]).norm()
                  << ", dW=" << inner_neighborhood.dW_ij_[first_bad_inner]
                  << ", e=";
        writeVectorForDiagnostic(inner_neighborhood.e_ij_[first_bad_inner]);
        std::cout << ", pos_j=";
        writeVectorForDiagnostic(pos[j]);
        std::cout << ", h_ratio_j=" << h_ratio[j]
                  << ", Vol_j=" << vol[j] << ", vel_j=";
        writeVectorForDiagnostic(velocity[j]);
    }
    std::cout << std::endl;

    for (size_t k = 0; k != wall_contact_relation.contact_configuration_.size(); ++k)
    {
        const Neighborhood &contact_neighborhood =
            wall_contact_relation.contact_configuration_[k][particle_index];
        BaseParticles *wall_particles = wall_contact_relation.contact_particles_[k];
        Vecd *wall_pos = wall_particles->getVariableDataByName<Vecd>("Position");
        Real *wall_vol = wall_particles->getVariableDataByName<Real>("VolumetricMeasure");

        Real min_contact_r = std::numeric_limits<Real>::max();
        size_t first_bad_contact = contact_neighborhood.current_size_;
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            const size_t j = contact_neighborhood.j_[n];
            min_contact_r = SMIN(min_contact_r, contact_neighborhood.r_ij_[n]);
            const bool finite_pair =
                isFiniteReal(contact_neighborhood.r_ij_[n]) &&
                isFiniteReal(contact_neighborhood.dW_ij_[n]) &&
                isFiniteVector(contact_neighborhood.e_ij_[n]) &&
                isFiniteReal(wall_vol[j]);
            if (!finite_pair && first_bad_contact == contact_neighborhood.current_size_)
            {
                first_bad_contact = n;
            }
        }

        std::cout << "[LocalRefinement][ViscousNeighborhood] " << stage
                  << ": contact_body=" << k
                  << ", contact_neighbors=" << contact_neighborhood.current_size_
                  << ", min_contact_r="
                  << (contact_neighborhood.current_size_ > 0 ? min_contact_r : Real(0.0));
        if (first_bad_contact != contact_neighborhood.current_size_)
        {
            const size_t j = contact_neighborhood.j_[first_bad_contact];
            std::cout << ", first_bad_contact_n=" << first_bad_contact
                      << ", j=" << j
                      << ", r=" << contact_neighborhood.r_ij_[first_bad_contact]
                      << ", dW=" << contact_neighborhood.dW_ij_[first_bad_contact]
                      << ", e=";
            writeVectorForDiagnostic(contact_neighborhood.e_ij_[first_bad_contact]);
            std::cout << ", wall_pos_j=";
            writeVectorForDiagnostic(wall_pos[j]);
            std::cout << ", wall_Vol_j=" << wall_vol[j];
        }
        std::cout << std::endl;
    }

    std::cout << "[LocalRefinement][ViscousNeighborhood] " << stage
              << ": pos_i=";
    writeVectorForDiagnostic(pos[particle_index]);
    std::cout << ", h_ratio_i=" << h_ratio[particle_index]
              << ", Vol_i=" << vol[particle_index]
              << ", vel_i=";
    writeVectorForDiagnostic(velocity[particle_index]);
    std::cout << std::endl;
}

/** @brief 校验 inner 和 contact relation 条目是否包含有限的距离/核数据。 */
inline bool reportRelationFiniteState(
    SPHBody &fluid_body, BaseInnerRelation &inner_relation,
    BaseContactRelation &wall_contact_relation, const std::string &stage)
{
    BaseParticles &particles = fluid_body.getBaseParticles();
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *h_ratio = particles.getVariableDataByName<Real>("SmoothingLengthRatio");
    const size_t total_particles = particles.TotalRealParticles();

    for (size_t i = 0; i != total_particles; ++i)
    {
        const Neighborhood &inner_neighborhood = inner_relation.inner_configuration_[i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            const size_t j = inner_neighborhood.j_[n];
            const bool finite_pair =
                j < total_particles &&
                isFiniteReal(inner_neighborhood.r_ij_[n]) &&
                isFiniteReal(inner_neighborhood.dW_ij_[n]) &&
                isFiniteVector(inner_neighborhood.e_ij_[n]);
            if (!finite_pair)
            {
                std::cout << "[LocalRefinement][Relation] " << stage
                          << ": first_bad_inner i=" << i
                          << ", n=" << n
                          << ", j=" << j
                          << ", r=" << inner_neighborhood.r_ij_[n]
                          << ", current_distance="
                          << (j < total_particles ? (pos[i] - pos[j]).norm() : Real(0.0))
                          << ", dW=" << inner_neighborhood.dW_ij_[n]
                          << ", e=";
                writeVectorForDiagnostic(inner_neighborhood.e_ij_[n]);
                std::cout << ", pos_i=";
                writeVectorForDiagnostic(pos[i]);
                if (j < total_particles)
                {
                    std::cout << ", pos_j=";
                    writeVectorForDiagnostic(pos[j]);
                    std::cout << ", h_ratio_i=" << h_ratio[i]
                              << ", h_ratio_j=" << h_ratio[j];
                }
                std::cout << std::endl;
                return false;
            }
        }
    }

    for (size_t k = 0; k != wall_contact_relation.contact_configuration_.size(); ++k)
    {
        BaseParticles *wall_particles = wall_contact_relation.contact_particles_[k];
        const size_t total_wall_particles = wall_particles->TotalRealParticles();
        for (size_t i = 0; i != total_particles; ++i)
        {
            const Neighborhood &contact_neighborhood = wall_contact_relation.contact_configuration_[k][i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                const size_t j = contact_neighborhood.j_[n];
                const bool finite_pair =
                    j < total_wall_particles &&
                    isFiniteReal(contact_neighborhood.r_ij_[n]) &&
                    isFiniteReal(contact_neighborhood.dW_ij_[n]) &&
                    isFiniteVector(contact_neighborhood.e_ij_[n]);
                if (!finite_pair)
                {
                    std::cout << "[LocalRefinement][Relation] " << stage
                              << ": first_bad_contact body=" << k
                              << ", i=" << i
                              << ", n=" << n
                              << ", j=" << j
                              << ", r=" << contact_neighborhood.r_ij_[n]
                              << ", dW=" << contact_neighborhood.dW_ij_[n]
                              << ", e=";
                    writeVectorForDiagnostic(contact_neighborhood.e_ij_[n]);
                    std::cout << std::endl;
                    return false;
                }
            }
        }
    }

    return true;
}

/** @brief 根据密度和速度场重新计算压力、质量和动量。 */
inline void syncEulerianStateFromDensityAndVelocity(SPHBody &fluid_body)
{
    BaseParticles &particles = fluid_body.getBaseParticles();
    Fluid &fluid = DynamicCast<Fluid>(&fluid_body, particles.getBaseMaterial());
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    Real *pressure = particles.getVariableDataByName<Real>("Pressure");
    Vecd *velocity = particles.getVariableDataByName<Vecd>("Velocity");
    Vecd *momentum = particles.registerStateVariableData<Vecd>("Momentum");
    const size_t total_particles = particles.TotalRealParticles();
    for (size_t i = 0; i != total_particles; ++i)
    {
        pressure[i] = fluid.getPressure(rho[i]);
        mass[i] = rho[i] * vol[i];
        momentum[i] = mass[i] * velocity[i];
    }
}

/** @brief 在普通 reload/lattice 生成后初始化 Eulerian 守恒量和力累加量。 */
inline void initializeEulerianStateFromDensityAndVelocity(SPHBody &fluid_body)
{
    BaseParticles &particles = fluid_body.getBaseParticles();
    Vecd *momentum_change_rate = particles.registerStateVariableData<Vecd>("MomentumChangeRate");
    Vecd *force = particles.registerStateVariableData<Vecd>("Force");
    Vecd *force_prior = particles.registerStateVariableData<Vecd>("ForcePrior");
    Vecd *viscous_force = particles.registerStateVariableData<Vecd>("ViscousForce");
    Vecd *previous_viscous_force = particles.registerStateVariableData<Vecd>("PreviousViscousForce");
    const size_t total_particles = particles.TotalRealParticles();
    syncEulerianStateFromDensityAndVelocity(fluid_body);
    for (size_t i = 0; i != total_particles; ++i)
    {
        momentum_change_rate[i] = Vecd::Zero();
        force[i] = Vecd::Zero();
        force_prior[i] = Vecd::Zero();
        viscous_force[i] = Vecd::Zero();
        previous_viscous_force[i] = Vecd::Zero();
    }
}

/** @brief 从已经包含密度和速度的 remapped reload 初始化守恒量。 */
inline void initializeEulerianStateFromRemappedReload(SPHBody &fluid_body)
{
    BaseParticles &particles = fluid_body.getBaseParticles();
    Fluid &fluid = DynamicCast<Fluid>(&fluid_body, particles.getBaseMaterial());
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *mass = particles.registerStateVariableData<Real>("Mass");
    Real *pressure = particles.getVariableDataByName<Real>("Pressure");
    Real *density_change_rate = particles.registerStateVariableData<Real>("DensityChangeRate");
    Real *mass_change_rate = particles.registerStateVariableData<Real>("MassChangeRate");
    Vecd *velocity = particles.getVariableDataByName<Vecd>("Velocity");
    Vecd *momentum = particles.registerStateVariableData<Vecd>("Momentum");
    Vecd *momentum_change_rate = particles.registerStateVariableData<Vecd>("MomentumChangeRate");
    Vecd *force = particles.registerStateVariableData<Vecd>("Force");
    Vecd *force_prior = particles.registerStateVariableData<Vecd>("ForcePrior");
    Vecd *viscous_force = particles.registerStateVariableData<Vecd>("ViscousForce");
    Vecd *previous_viscous_force = particles.registerStateVariableData<Vecd>("PreviousViscousForce");
    const size_t total_particles = particles.TotalRealParticles();

    for (size_t i = 0; i != total_particles; ++i)
    {
        if (!isFiniteReal(rho[i]) || rho[i] <= Real(0.0))
        {
            throw std::runtime_error("Remapped reload has non-positive or non-finite Density.");
        }
        if (!isFiniteVector(velocity[i]))
        {
            throw std::runtime_error("Remapped reload has non-finite Velocity.");
        }

        const Real eos_pressure = fluid.getPressure(rho[i]);
        if (!isFiniteReal(pressure[i]))
        {
            pressure[i] = eos_pressure;
        }
        mass[i] = rho[i] * vol[i];
        momentum[i] = mass[i] * velocity[i];
        density_change_rate[i] = Real(0.0);
        mass_change_rate[i] = Real(0.0);
        momentum_change_rate[i] = Vecd::Zero();
        force[i] = Vecd::Zero();
        force_prior[i] = Vecd::Zero();
        viscous_force[i] = Vecd::Zero();
        previous_viscous_force[i] = Vecd::Zero();
    }
}

/** @brief 支持任意局部粒子间距因子的 AdaptiveWithinShape。 */
class AdaptiveWithinShapeBySpacingFactor : public AdaptiveWithinShape
{
  public:
    AdaptiveWithinShapeBySpacingFactor(Real global_resolution, Real h_spacing_ratio,
                                       Real refinement_to_global, int local_refinement_level,
                                       Real spacing_factor)
        : AdaptiveWithinShape(global_resolution, h_spacing_ratio, refinement_to_global,
                              local_refinement_level)
    {
        if (!std::isfinite(spacing_factor) || !(spacing_factor > Real(0.0)) ||
            spacing_factor > Real(1.0))
        {
            throw std::runtime_error(
                "Config value must be in (0, 1]: localrefinement.spacing_factor");
        }
        spacing_min_ = spacing_ref_ * spacing_factor;
        Vol_min_ = std::pow(spacing_min_, Dimensions);
        h_ratio_max_ = spacing_ref_ / spacing_min_;
        finest_spacing_bound_ = spacing_min_ + Eps;
        coarsest_spacing_bound_ = spacing_ref_ - Eps;
        max_cut_off_radius_ = kernel_ptr_->KernelSize() * h_ref_;
    }
};

// ---------------------------------------------------------------------------
//  Deterministic replacements for rand_uniform-driven routines
//
//  Why: SPHinXsys 的 rand_uniform 基于 system_clock::now() 重建 RNG，跨硬件不
//  可重现。local-refinement 加密几何阶段依赖 rand_uniform 决定 lattice 接受
//  率和 relaxation 初始扰动，超算 vs 本地会生成粒子数不同的 reload，直接导
//  致 Stage 2 dt 坍缩。此处在 case 级用 SplitMix64 风格的确定性哈希替代。
// ---------------------------------------------------------------------------

inline std::uint64_t cylinderLgSplitMix64(std::uint64_t x)
{
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

template <typename... Keys>
inline Real cylinderLgDeterministicUniform01(std::uint64_t seed, Keys... keys)
{
    std::uint64_t state = seed;
    const std::uint64_t key_values[] = {static_cast<std::uint64_t>(keys)...};
    for (std::uint64_t k : key_values)
    {
        state = cylinderLgSplitMix64(state ^ cylinderLgSplitMix64(k));
    }
    // 53-bit mantissa → [0, 1)
    return Real(state >> 11) * (Real(1.0) / Real(1ULL << 53));
}

/** @brief Tag：使用确定性哈希替代 lattice accept/reject 的粒子生成器。 */
struct DeterministicLattice
{
};

/** @brief 特化：adaptive lattice 生成器，用 cell index hash 取代 rand_uniform。 */
template <>
class ParticleGenerator<BaseParticles, DeterministicLattice>
    : public ParticleGenerator<BaseParticles>,
      public GeneratingMethod<Lattice>
{
  public:
    ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles, Shape &target_shape)
        : ParticleGenerator<BaseParticles>(sph_body, base_particles),
          GeneratingMethod<Lattice>(sph_body),
          target_shape_(target_shape) {}

    ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles)
        : ParticleGenerator(sph_body, base_particles, sph_body.getInitialShape()) {}

    ~ParticleGenerator() override = default;

    void prepareGeometricData() override
    {
        Mesh mesh(domain_bounds_, lattice_spacing_, 0);
        Real particle_volume = lattice_spacing_ * lattice_spacing_;
        Arrayi number_of_lattices = mesh.AllCells();
        for (int i = 0; i < number_of_lattices[0]; ++i)
        {
            for (int j = 0; j < number_of_lattices[1]; ++j)
            {
                Vecd particle_position = mesh.CellPositionFromIndex(Arrayi(i, j));
                if (initial_shape_.checkContain(particle_position))
                {
                    addDeterministicPosition(particle_position, particle_volume, i, j);
                }
            }
        }
    }

  protected:
    Shape &target_shape_;
    static constexpr std::uint64_t kLatticeSeed = 0xC0FFEE1234567ULL;

    void addDeterministicPosition(const Vecd &position, Real volume, int i, int j)
    {
        Real local_particle_spacing = sph_adaptation_.getLocalSpacing(target_shape_, position);
        Real local_particle_volume_ratio =
            std::pow(lattice_spacing_ / local_particle_spacing, Dimensions);
        const Real u = cylinderLgDeterministicUniform01(
            kLatticeSeed,
            static_cast<std::uint64_t>(static_cast<std::uint32_t>(i)),
            static_cast<std::uint64_t>(static_cast<std::uint32_t>(j)));
        if (u < local_particle_volume_ratio)
        {
            ParticleGenerator<BaseParticles>::addPositionAndVolumetricMeasure(
                position, volume / local_particle_volume_ratio);
        }
    }
};

/** @brief 确定性 relaxation 初始扰动，替换 RandomizeParticlePosition 以消除时钟随机源。 */
class DeterministicRandomizeParticlePosition : public LocalDynamics
{
  public:
    explicit DeterministicRandomizeParticlePosition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          randomize_scale_(sph_body.getSPHAdaptation().MinimumSpacing()) {}

    void update(size_t index_i, Real dt = 0.0)
    {
        Vecd &pos_n_i = pos_[index_i];
        for (int k = 0; k < pos_n_i.size(); ++k)
        {
            const Real u = cylinderLgDeterministicUniform01(
                kRandomizeSeed,
                static_cast<std::uint64_t>(index_i),
                static_cast<std::uint64_t>(static_cast<std::uint32_t>(k)));
            pos_n_i[k] += dt * (Real(2.0) * u - Real(1.0)) * randomize_scale_;
        }
    }

  private:
    Vecd *pos_;
    Real randomize_scale_;
    static constexpr std::uint64_t kRandomizeSeed = 0xBADC0DEABCDEFULL;
};

/** @brief 流体域形状：带 sponge 扩展的通道，并扣除圆柱区域。 */
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon outer_boundary(createWaterBlockShape());
        add<MultiPolygonShape>(outer_boundary, "OuterBoundary");
        MultiPolygon circle(cylinder_center, cylinder_radius, 100);
        subtract<MultiPolygonShape>(circle);
    }
};

/** @brief 用于 wall contact 和受力积分的固体圆柱几何。 */
class Cylinder : public MultiPolygonShape
{
  public:
    explicit Cylinder(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        multi_polygon_.addACircle(cylinder_center, cylinder_radius, 100, GeometricOps::add);
    }
};

/** @brief 非反射边界修正，并对局部加密 sponge 粒子做专门处理。 */
class FarFieldBoundary : public fluid_dynamics::NonReflectiveBoundaryCorrection
{
  public:
    explicit FarFieldBoundary(BaseInnerRelation &inner_relation)
        : fluid_dynamics::NonReflectiveBoundaryCorrection(inner_relation)
    {
        rho_farfield_ = rho0_f;
        sound_speed_ = c_f;
        vel_farfield_ = Vecd(U_f, 0.0);
    };
    virtual ~FarFieldBoundary() {};

    /** @brief 按局部分支和声学条件更新单个远场粒子状态。 */
    void update(size_t index_i, Real dt = 0.0)
    {
        if (enable_local_refinement && isInsideLocalRefinementOuterBuffer(pos_[index_i]))
        {
            rho_[index_i] = rho_farfield_;
            p_[index_i] = fluid_.getPressure(rho_[index_i]);
            vel_[index_i] = vel_farfield_;
            mass_[index_i] = rho_[index_i] * Vol_[index_i];
            mom_[index_i] = mass_[index_i] * vel_[index_i];
            return;
        }

        if (!enable_local_refinement || (indicator_[index_i] != 1 && smeared_surface_[index_i] != 1))
        {
            fluid_dynamics::NonReflectiveBoundaryCorrection::update(index_i, dt);
            return;
        }

        Real velocity_farfield_normal = vel_farfield_.dot(n_[index_i]);
        Real velocity_boundary_normal = vel_[index_i].dot(n_[index_i]);

        if (n_[index_i][0] <= 0.0 || std::fabs(n_[index_i][1]) > std::fabs(n_[index_i][0]))
        {
            if (std::fabs(velocity_boundary_normal) >= sound_speed_)
            {
                vel_[index_i] = vel_farfield_;
                rho_[index_i] = rho_farfield_;
                mass_[index_i] = rho_[index_i] * Vol_[index_i];
                mom_[index_i] = mass_[index_i] * vel_[index_i];
            }
            if (std::fabs(velocity_boundary_normal) < sound_speed_)
            {
                const Real bounded_weight =
                    SMIN(Real(1.0), SMAX(Real(0.0), inner_weight_summation_[index_i]));
                rho_[index_i] =
                    rho_average_[index_i] * bounded_weight + rho_farfield_ * (1.0 - bounded_weight);
                p_[index_i] = fluid_.getPressure(rho_[index_i]);
                Real vel_normal =
                    vel_normal_average_[index_i] * bounded_weight +
                    velocity_farfield_normal * (1.0 - bounded_weight);
                vel_[index_i] = vel_normal * n_[index_i] +
                                (vel_farfield_ - velocity_farfield_normal * n_[index_i]);
                mass_[index_i] = rho_[index_i] * Vol_[index_i];
                mom_[index_i] = mass_[index_i] * vel_[index_i];
            }
        }
        else
        {
            if (std::fabs(velocity_boundary_normal) >= sound_speed_)
            {
                rho_[index_i] = rho_average_[index_i] + TinyReal;
                vel_[index_i] = vel_average_[index_i];
                mass_[index_i] = rho_[index_i] * Vol_[index_i];
                mom_[index_i] = mass_[index_i] * vel_[index_i];
            }

            if (std::fabs(velocity_boundary_normal) < sound_speed_)
            {
                const Real bounded_weight =
                    SMIN(Real(1.0), SMAX(Real(0.0), inner_weight_summation_[index_i]));
                rho_[index_i] =
                    rho_average_[index_i] * bounded_weight + rho_farfield_ * (1.0 - bounded_weight);
                p_[index_i] = fluid_.getPressure(rho_[index_i]);
                Real vel_normal =
                    vel_normal_average_[index_i] * bounded_weight +
                    velocity_farfield_normal * (1.0 - bounded_weight);
                vel_[index_i] = vel_normal * n_[index_i] + vel_tangential_average_[index_i];
                mass_[index_i] = rho_[index_i] * Vol_[index_i];
                mom_[index_i] = mass_[index_i] * vel_[index_i];
            }
        }
    }
};

/** @brief 调整 wall normal 方向，使其与 fluid-wall 粒子对方向一致。 */
inline Vecd orientedWallNormal(const Vecd &e_ij, const Vecd &wall_normal)
{
    return SGN(e_ij.dot(wall_normal)) * wall_normal;
}

/** @brief 将力投影到 wall normal 方向，并容忍退化法向量。 */
inline Vecd projectForceToWallNormal(const Vecd &force, const Vecd &wall_normal)
{
    const Real normal_norm = wall_normal.norm();
    if (normal_norm <= Real(1.0e-12))
    {
        return force;
    }
    const Vecd unit_normal = wall_normal / normal_norm;
    return force.dot(unit_normal) * unit_normal;
}

/** @brief 根据 pressure-contact 模式在完整向量力和法向力之间切换。 */
inline Vecd applyPressureContactNormalMode(const Vecd &force, const Vecd &wall_normal,
                                           bool pressure_contact_only_normal)
{
    return pressure_contact_only_normal ? projectForceToWallNormal(force, wall_normal) : force;
}

/** @brief 只保留壁面法向压力贡献的 Eulerian wall-contact 压力积分。 */
template <class RiemannSolverType>
class EulerianPressureIntegration1stHalfWallContactNormalOnly
    : public fluid_dynamics::BaseEulerianIntegrationWithWall
{
  public:
    /** @brief 构造 wall-contact 积分器并注册压力壁面力输出变量。 */
    explicit EulerianPressureIntegration1stHalfWallContactNormalOnly(
        BaseContactRelation &wall_contact_relation, bool pressure_contact_only_normal,
        Real limiter_parameter = 15.0)
        : fluid_dynamics::BaseEulerianIntegrationWithWall(wall_contact_relation),
          riemann_solver_(this->fluid_, this->fluid_, limiter_parameter),
          wall_force_pressure_(this->particles_->registerStateVariableData<Vecd>("WallForcePressure")),
          wall_force_pressure_impulse_(this->particles_->registerStateVariableData<Vecd>("WallForcePressureImpulse")),
          pressure_contact_only_normal_(pressure_contact_only_normal)
    {
    }

    /** @brief 累加单个流体粒子的壁面压力和对流通量贡献。 */
    void interaction(size_t index_i, Real dt = 0.0)
    {
        FluidStateIn state_i(this->rho_[index_i], this->vel_[index_i], this->p_[index_i]);
        Vecd momentum_change_rate = Vecd::Zero();
        Vecd pressure_force = Vecd::Zero();
        for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
        {
            Vecd *n_k = this->wall_n_[k];
            Real *Vol_k = this->wall_Vol_[k];
            Vecd *vel_ave_k = this->wall_vel_ave_[k];
            Neighborhood &wall_neighborhood = (*this->contact_configuration_[k])[index_i];
            for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
            {
                const size_t index_j = wall_neighborhood.j_[n];
                const Vecd &e_ij = wall_neighborhood.e_ij_[n];
                const Real dW_ijV_j = wall_neighborhood.dW_ij_[n] * Vol_k[index_j];

                Real rho_j_in_wall = state_i.rho_;
                Real p_j_in_wall = state_i.p_;
                Vecd vel_j_in_wall = 2.0 * vel_ave_k[index_j] - state_i.vel_;
                FluidStateIn state_j(rho_j_in_wall, vel_j_in_wall, p_j_in_wall);
                FluidStateOut interface_state = riemann_solver_.InterfaceState(state_i, state_j, n_k[index_j]);
                const Matd convect_flux = interface_state.rho_ * interface_state.vel_ * interface_state.vel_.transpose();
                const Vecd convective_force =
                    -2.0 * this->Vol_[index_i] * convect_flux * e_ij * dW_ijV_j;
                const Vecd raw_pressure_force =
                    -2.0 * this->Vol_[index_i] * interface_state.p_ * e_ij * dW_ijV_j;
                const Vecd face_to_fluid_n = orientedWallNormal(e_ij, n_k[index_j]);
                const Vecd projected_pressure_force = applyPressureContactNormalMode(
                    raw_pressure_force, face_to_fluid_n, pressure_contact_only_normal_);

                momentum_change_rate += convective_force + projected_pressure_force;
                pressure_force += projected_pressure_force;
            }
        }
        this->dmom_dt_[index_i] += momentum_change_rate;
        wall_force_pressure_[index_i] = pressure_force;
        wall_force_pressure_impulse_[index_i] += pressure_force * dt;
    }

  protected:
    RiemannSolverType riemann_solver_;
    Vecd *wall_force_pressure_;
    Vecd *wall_force_pressure_impulse_;
    bool pressure_contact_only_normal_;
};

/** @brief 将 inner Eulerian 压力积分与法向 wall-contact 压力积分组合。 */
template <class RiemannSolverType>
class EulerianPressureIntegration1stHalfWithWallNormalOnly
    : public fluid_dynamics::EulerianIntegration1stHalf<Inner<>, RiemannSolverType>
{
    using InnerIntegration = fluid_dynamics::EulerianIntegration1stHalf<Inner<>, RiemannSolverType>;

  public:
    /** @brief 构造 inner 和 wall-contact 两部分的一阶压力积分器。 */
    explicit EulerianPressureIntegration1stHalfWithWallNormalOnly(
        BaseInnerRelation &inner_relation, BaseContactRelation &wall_contact_relation,
        bool pressure_contact_only_normal, Real limiter_parameter = 15.0)
        : InnerIntegration(inner_relation, limiter_parameter),
          wall_contact_integration_(wall_contact_relation, pressure_contact_only_normal, limiter_parameter)
    {
    }

    using InnerIntegration::update;

    /** @brief 保持与 SPHinXsys dynamics 接口一致；本算例无需额外初始化。 */
    void initialization(size_t index_i, Real dt = 0.0) {}

    /** @brief 先执行 inner 压力积分，再执行壁面压力积分。 */
    void interaction(size_t index_i, Real dt = 0.0)
    {
        InnerIntegration::interaction(index_i, dt);
        wall_contact_integration_.interaction(index_i, dt);
    }

  protected:
    EulerianPressureIntegration1stHalfWallContactNormalOnly<RiemannSolverType>
        wall_contact_integration_;
};

using EulerianPressureRelaxationWithWallNormalOnly =
    EulerianPressureIntegration1stHalfWithWallNormalOnly<AcousticRiemannSolver>;

/** @brief 从流体侧累积壁面黏性力贡献，供 force output 使用。 */
class FluidWallViscousForceRecorder : public LocalDynamics, public DataDelegateContact
{
  public:
    /** @brief 构造黏性力记录器并缓存 wall contact 所需的体积和平均速度。 */
    explicit FluidWallViscousForceRecorder(BaseContactRelation &wall_contact_relation)
        : LocalDynamics(wall_contact_relation.getSPHBody()),
          DataDelegateContact(wall_contact_relation),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
          wall_force_viscous_(particles_->registerStateVariableData<Vecd>("WallForceViscous")),
          mu_(particles_), kernel_correction_(particles_),
          smoothing_length_(getSPHAdaptation().ReferenceSmoothingLength())
    {
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            Solid &solid_material = DynamicCast<Solid>(this, contact_particles_[k]->getBaseMaterial());
            wall_vel_ave_.push_back(solid_material.AverageVelocity(contact_particles_[k]));
            wall_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        }
    }

    /** @brief 计算单个流体粒子与壁面邻居产生的黏性力。 */
    void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd viscous_force = Vecd::Zero();
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Vecd *vel_ave_k = wall_vel_ave_[k];
            Real *wall_Vol_k = wall_Vol_[k];
            const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                const size_t index_j = contact_neighborhood.j_[n];
                const Real r_ij = contact_neighborhood.r_ij_[n];
                const Vecd &e_ij = contact_neighborhood.e_ij_[n];
                const Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * wall_Vol_k[index_j];

                const Vecd vel_derivative =
                    2.0 * (vel_[index_i] - vel_ave_k[index_j]) /
                    (r_ij + 0.01 * smoothing_length_);
                viscous_force +=
                    2.0 * e_ij.dot(kernel_correction_(index_i) * e_ij) *
                    mu_(index_i, index_i) * vel_derivative * dW_ijV_j;
            }
        }

        wall_force_viscous_[index_i] = viscous_force * Vol_[index_i];
    }

  protected:
    Real *Vol_;
    Vecd *vel_;
    Vecd *wall_force_viscous_;
    fluid_dynamics::FixedViscosity mu_;
    LinearGradientCorrection kernel_correction_;
    Real smoothing_length_;
    StdVec<Vecd *> wall_vel_ave_;
    StdVec<Real *> wall_Vol_;
};

/** @brief 在每个输出窗口开始时清零流体侧壁面力累加量。 */
class FluidWallForceReset : public LocalDynamics
{
  public:
    /** @brief 获取压力、压力冲量和黏性力变量指针。 */
    explicit FluidWallForceReset(SPHBody &fluid_body)
        : LocalDynamics(fluid_body),
          wall_force_pressure_(particles_->getVariableDataByName<Vecd>("WallForcePressure")),
          wall_force_pressure_impulse_(particles_->getVariableDataByName<Vecd>("WallForcePressureImpulse")),
          wall_force_viscous_(particles_->getVariableDataByName<Vecd>("WallForceViscous"))
    {
    }

    /** @brief 清零单个粒子的壁面力相关状态。 */
    void update(size_t index_i, Real dt = 0.0)
    {
        wall_force_pressure_[index_i] = Vecd::Zero();
        wall_force_pressure_impulse_[index_i] = Vecd::Zero();
        wall_force_viscous_[index_i] = Vecd::Zero();
    }

  protected:
    Vecd *wall_force_pressure_;
    Vecd *wall_force_pressure_impulse_;
    Vecd *wall_force_viscous_;
};

/** @brief 将流体侧壁面压力力和黏性力写入 output/fluid_wall_force.dat。 */
class FluidWallForceRecording
{
  public:
    /** @brief 构造输出器，并在续算时决定是否追加已有文件。 */
    explicit FluidWallForceRecording(SPHBody &fluid_body, bool is_restart = false)
        : fluid_body_(fluid_body),
          physical_time_(*fluid_body.getSPHSystem().getSystemVariableDataByName<Real>("PhysicalTime")),
          file_path_(std::filesystem::current_path() / "output" / "fluid_wall_force.dat")
    {
        header_written_ = cylinder_lg::shouldAppendRestartOutputFile(file_path_, is_restart);
    }

    /** @brief 汇总所有流体粒子的壁面力并写出当前输出步。 */
    void writeToFile(size_t iteration_step, Real averaging_time)
    {
        std::filesystem::create_directories(file_path_.parent_path());
        std::ofstream out_file(file_path_, header_written_ ? std::ios::app : std::ios::out);
        if (!header_written_)
        {
            out_file << "\"run_time\"   \"iteration\"   \"pressure_force_x\"   "
                     << "\"pressure_force_y\"   \"viscous_force_x\"   \"viscous_force_y\"\n";
            header_written_ = true;
        }

        Vecd pressure_force = Vecd::Zero();
        Vecd viscous_force = Vecd::Zero();
        BaseParticles &particles = fluid_body_.getBaseParticles();
        Vecd *wall_force_pressure_impulse = particles.getVariableDataByName<Vecd>("WallForcePressureImpulse");
        Vecd *wall_force_viscous = particles.getVariableDataByName<Vecd>("WallForceViscous");
        const size_t total_particles = particles.TotalRealParticles();
        for (size_t i = 0; i != total_particles; ++i)
        {
            pressure_force += wall_force_pressure_impulse[i];
            viscous_force += wall_force_viscous[i];
        }
        if (averaging_time > Real(0.0))
        {
            pressure_force /= averaging_time;
        }

        out_file << std::fixed << std::setprecision(12)
                 << physical_time_ << "   " << iteration_step << "   "
                 << pressure_force[0] << "   " << pressure_force[1] << "   "
                 << viscous_force[0] << "   " << viscous_force[1] << "\n";
    }

  protected:
    SPHBody &fluid_body_;
    Real &physical_time_;
    std::filesystem::path file_path_;
    bool header_written_{false};
};
} // namespace SPH

#endif // CYLINDER_LG_CALCULATION_HPP
