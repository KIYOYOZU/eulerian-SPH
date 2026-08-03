#include "sphinxsys.h"

/**
 * @file local_refinement_remap.cpp
 * @brief 将粗网格 restart 状态映射到局部加密 reload 粒子，并输出守恒检查报告。
 */

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace
{
using SPH::Real;
using tinyxml2::XMLElement;
using tinyxml2::XMLDocument;

/** @brief 重映射工具内部使用的二维标量向量。 */
struct Vec2
{
    Real x = 0.0;
    Real y = 0.0;
};

/** @brief 从粗网格 restart XML 读取的源粒子状态。 */
struct CoarseParticle
{
    unsigned int original_id = 0;
    Vec2 position;
    Real volume = 0.0;
    Real density = 0.0;
    Real pressure = 0.0;
    Vec2 velocity;
    Real mass = 0.0;
    Vec2 momentum;
};

/** @brief 从目标 reload XML 读取并写回的加密粒子状态。 */
struct TargetParticle
{
    unsigned int original_id = 0;
    Vec2 position;
    Real volume = 0.0;
    Real smoothing_length_ratio = 1.0;
    Real density = 0.0;
    Real pressure = 0.0;
    Vec2 velocity;
    Real mass = 0.0;
    Vec2 momentum;
};

/** @brief 用于比较 remap 前后守恒量的分桶统计。 */
struct BucketStats
{
    size_t count = 0;
    Real mass = 0.0;
    Vec2 momentum;
};

/** @brief remap 后守恒检查的固定容许值。 */
constexpr Real REMAP_MASS_TOL_TOTAL = 1.0e-6;
constexpr Real REMAP_MASS_TOL_REGION = 1.0e-3;
constexpr Real REMAP_MOMENTUM_TOL_TOTAL = 1.0e-4;
constexpr Real REMAP_MOMENTUM_TOL_REGION = 1.0e-3;

/** @brief 命令行参数。 */
struct Options
{
    fs::path source_restart;
    fs::path target_water_reload;
    fs::path target_cylinder_reload;
    fs::path output_dir;
    fs::path report_path;
    Real rho0 = 1.225;
    Real sound_speed = 0.6;
    Real velocity_scale = 0.06;
    Real region_x_min = 0.3;
    Real region_x_max = 0.55;
    Real region_y_min = 0.25;
    Real region_y_max = 0.55;
    int neighbors = 8;
};

/** @brief 读取需要下一项作为值的命令行参数。 */
std::string requireArgument(int argc, char *argv[], int &i, const std::string &name)
{
    if (i + 1 >= argc)
    {
        throw std::runtime_error("Missing value after " + name);
    }
    ++i;
    return argv[i];
}

/** @brief 判断字符串是否以指定前缀开头。 */
bool startsWith(const std::string &value, const std::string &prefix)
{
    return value.rfind(prefix, 0) == 0;
}

/** @brief 解析 Real 命令行参数，并拒绝非法尾随字符。 */
Real parseReal(const std::string &value, const std::string &name)
{
    size_t parsed = 0;
    try
    {
        const Real result = static_cast<Real>(std::stod(value, &parsed));
        if (parsed != value.size())
        {
            throw std::invalid_argument("trailing characters");
        }
        return result;
    }
    catch (const std::exception &)
    {
        throw std::runtime_error("Invalid real value for " + name + ": " + value);
    }
}

/** @brief 解析整数命令行参数，并拒绝非法尾随字符。 */
int parseInt(const std::string &value, const std::string &name)
{
    size_t parsed = 0;
    try
    {
        const int result = std::stoi(value, &parsed);
        if (parsed != value.size())
        {
            throw std::invalid_argument("trailing characters");
        }
        return result;
    }
    catch (const std::exception &)
    {
        throw std::runtime_error("Invalid integer value for " + name + ": " + value);
    }
}

/** @brief 解析并校验 coarse-to-refined remap 所需的命令行参数。 */
Options parseOptions(int argc, char *argv[])
{
    Options options;
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        auto read_option = [&](const std::string &long_name) -> std::string
        {
            const std::string with_equal = long_name + "=";
            if (startsWith(arg, with_equal))
            {
                return arg.substr(with_equal.size());
            }
            if (arg == long_name)
            {
                return requireArgument(argc, argv, i, long_name);
            }
            return "";
        };

        std::string value;
        if (!(value = read_option("--source-restart")).empty())
        {
            options.source_restart = value;
        }
        else if (!(value = read_option("--target-water-reload")).empty())
        {
            options.target_water_reload = value;
        }
        else if (!(value = read_option("--target-cylinder-reload")).empty())
        {
            options.target_cylinder_reload = value;
        }
        else if (!(value = read_option("--output-dir")).empty())
        {
            options.output_dir = value;
        }
        else if (!(value = read_option("--report")).empty())
        {
            options.report_path = value;
        }
        else if (!(value = read_option("--rho0")).empty())
        {
            options.rho0 = parseReal(value, "--rho0");
        }
        else if (!(value = read_option("--sound-speed")).empty())
        {
            options.sound_speed = parseReal(value, "--sound-speed");
        }
        else if (!(value = read_option("--velocity-scale")).empty())
        {
            options.velocity_scale = parseReal(value, "--velocity-scale");
        }
        else if (!(value = read_option("--region-x-min")).empty())
        {
            options.region_x_min = parseReal(value, "--region-x-min");
        }
        else if (!(value = read_option("--region-x-max")).empty())
        {
            options.region_x_max = parseReal(value, "--region-x-max");
        }
        else if (!(value = read_option("--region-y-min")).empty())
        {
            options.region_y_min = parseReal(value, "--region-y-min");
        }
        else if (!(value = read_option("--region-y-max")).empty())
        {
            options.region_y_max = parseReal(value, "--region-y-max");
        }
        else if (!(value = read_option("--neighbors")).empty())
        {
            options.neighbors = parseInt(value, "--neighbors");
        }
        else
        {
            throw std::runtime_error("Unknown option: " + arg);
        }
    }

    if (options.source_restart.empty() || options.target_water_reload.empty() ||
        options.target_cylinder_reload.empty() || options.output_dir.empty() ||
        options.report_path.empty())
    {
        throw std::runtime_error(
            "Usage: local_refinement_remap --source-restart <Restart_*.xml> "
            "--target-water-reload <WaterBlock_rld.xml> --target-cylinder-reload <Cylinder_rld.xml> "
            "--output-dir <case/reload> --report <report.csv>");
    }
    if (options.neighbors <= 0)
    {
        throw std::runtime_error("--neighbors must be positive.");
    }
    if (!(options.region_x_min < options.region_x_max) ||
        !(options.region_y_min < options.region_y_max))
    {
        throw std::runtime_error("Invalid refinement region bounds.");
    }
    if (!(options.rho0 > 0.0) || !(options.sound_speed > 0.0) || !(options.velocity_scale > 0.0))
    {
        throw std::runtime_error("--rho0, --sound-speed and --velocity-scale must be positive.");
    }
    return options;
}

/** @brief 加载 XML 文件，失败时抛出包含路径的错误。 */
void loadXmlDocument(XMLDocument &doc, const fs::path &path)
{
    const auto result = doc.LoadFile(path.string().c_str());
    if (result != tinyxml2::XML_SUCCESS)
    {
        throw std::runtime_error("Failed to load XML file: " + path.string());
    }
}

/** @brief 读取 XML 必需属性，缺失时抛出带路径的错误。 */
const char *requiredAttribute(XMLElement *element, const char *name, const fs::path &path)
{
    const char *value = element->Attribute(name);
    if (value == nullptr)
    {
        throw std::runtime_error("Missing XML attribute '" + std::string(name) + "' in " + path.string());
    }
    return value;
}

/** @brief 从 XML 属性读取 Real 值。 */
Real readReal(XMLElement *element, const char *name, const fs::path &path)
{
    return parseReal(requiredAttribute(element, name, path), name);
}

/** @brief 从 XML 属性读取无符号整数值。 */
unsigned int readUnsigned(XMLElement *element, const char *name, const fs::path &path)
{
    size_t parsed = 0;
    const std::string value = requiredAttribute(element, name, path);
    try
    {
        const unsigned long parsed_value = std::stoul(value, &parsed);
        if (parsed != value.size())
        {
            throw std::invalid_argument("trailing characters");
        }
        return static_cast<unsigned int>(parsed_value);
    }
    catch (const std::exception &)
    {
        throw std::runtime_error("Invalid unsigned value for " + std::string(name) + ": " + value);
    }
}

/** @brief 从 XML 属性读取逗号分隔的二维向量。 */
Vec2 readVec2(XMLElement *element, const char *name, const fs::path &path)
{
    const std::string value = requiredAttribute(element, name, path);
    const size_t comma_pos = value.find(',');
    if (comma_pos == std::string::npos)
    {
        throw std::runtime_error("Invalid Vec2 value for " + std::string(name) + ": " + value);
    }
    return {parseReal(value.substr(0, comma_pos), name),
            parseReal(value.substr(comma_pos + 1), name)};
}

/** @brief 以 Real 精度格式化数值用于 XML 写回。 */
std::string formatReal(Real value)
{
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<Real>::digits10) << value;
    return out.str();
}

/** @brief 将二维向量格式化为 SPHinXsys XML 属性文本。 */
std::string formatVec2(const Vec2 &value)
{
    return formatReal(value.x) + ", " + formatReal(value.y);
}

/** @brief 二维向量加法。 */
Vec2 operator+(const Vec2 &a, const Vec2 &b)
{
    return {a.x + b.x, a.y + b.y};
}

/** @brief 二维向量减法。 */
Vec2 operator-(const Vec2 &a, const Vec2 &b)
{
    return {a.x - b.x, a.y - b.y};
}

/** @brief 标量乘二维向量。 */
Vec2 operator*(Real scale, const Vec2 &value)
{
    return {scale * value.x, scale * value.y};
}

/** @brief 二维向量除以标量。 */
Vec2 operator/(const Vec2 &value, Real scale)
{
    return {value.x / scale, value.y / scale};
}

/** @brief 计算二维向量范数。 */
Real norm(const Vec2 &value)
{
    return std::sqrt(value.x * value.x + value.y * value.y);
}

/** @brief 计算两个二维点之间的平方距离。 */
Real squaredDistance(const Vec2 &a, const Vec2 &b)
{
    const Real dx = a.x - b.x;
    const Real dy = a.y - b.y;
    return dx * dx + dy * dy;
}

/** @brief 空间哈希网格中的整数单元坐标。 */
struct CellCoord
{
    long long ix = 0;
    long long iy = 0;

    bool operator==(const CellCoord &other) const
    {
        return ix == other.ix && iy == other.iy;
    }
};

/** @brief 为 CellCoord 提供 unordered_map 哈希。 */
struct CellCoordHash
{
    size_t operator()(const CellCoord &coord) const
    {
        const size_t h1 = std::hash<long long>{}(coord.ix);
        const size_t h2 = std::hash<long long>{}(coord.iy);
        return h1 ^ (h2 + size_t(0x9e3779b97f4a7c15ULL) + (h1 << 6) + (h1 >> 2));
    }
};

/** @brief 基于粗网格粒子位置的空间索引，用于近邻查找。 */
struct SourceSpatialIndex
{
    const std::vector<CoarseParticle> &source;
    Real cell_size;
    std::unordered_map<CellCoord, std::vector<size_t>, CellCoordHash> cells;

    /** @brief 按给定 cell size 构建源粒子空间索引。 */
    SourceSpatialIndex(const std::vector<CoarseParticle> &source_particles, Real spacing)
        : source(source_particles), cell_size(spacing)
    {
        if (!(cell_size > 0.0))
        {
            throw std::runtime_error("Spatial index cell size must be positive.");
        }
        for (size_t i = 0; i != source.size(); ++i)
        {
            cells[cellIndex(source[i].position)].push_back(i);
        }
    }

    /** @brief 返回指定位置所在的哈希网格单元。 */
    CellCoord cellIndex(const Vec2 &position) const
    {
        return {static_cast<long long>(std::floor(position.x / cell_size)),
                static_cast<long long>(std::floor(position.y / cell_size))};
    }

    /** @brief 查找距离指定位置最近的若干源粒子索引。 */
    std::vector<size_t> nearest(const Vec2 &position, int count) const
    {
        const CellCoord center = cellIndex(position);
        std::vector<std::pair<Real, size_t>> distances;
        distances.reserve(static_cast<size_t>(count) * 4);

        int ring = 0;
        while ((distances.size() < static_cast<size_t>(count) || ring < 2) && ring <= 128)
        {
            for (long long ix = center.ix - ring; ix <= center.ix + ring; ++ix)
            {
                for (long long iy = center.iy - ring; iy <= center.iy + ring; ++iy)
                {
                    if (ring != 0 &&
                        ix != center.ix - ring && ix != center.ix + ring &&
                        iy != center.iy - ring && iy != center.iy + ring)
                    {
                        continue;
                    }
                    const auto cell_iter = cells.find(CellCoord{ix, iy});
                    if (cell_iter == cells.end())
                    {
                        continue;
                    }
                    for (size_t index : cell_iter->second)
                    {
                        distances.emplace_back(squaredDistance(source[index].position, position), index);
                    }
                }
            }
            ++ring;
        }

        if (distances.empty())
        {
            throw std::runtime_error("Spatial index found no source neighbors.");
        }

        const size_t neighbor_count = std::min<size_t>(static_cast<size_t>(count), distances.size());
        if (neighbor_count < distances.size())
        {
            std::nth_element(
                distances.begin(),
                distances.begin() + static_cast<std::ptrdiff_t>(neighbor_count),
                distances.end());
        }
        distances.resize(neighbor_count);
        std::sort(distances.begin(), distances.end());

        std::vector<size_t> result;
        result.reserve(neighbor_count);
        for (const auto &entry : distances)
        {
            result.push_back(entry.second);
        }
        return result;
    }
};

/** @brief 判断位置是否落入用户指定的局部加密区域。 */
bool inRegion(const Options &options, const Vec2 &position)
{
    return position.x >= options.region_x_min && position.x <= options.region_x_max &&
           position.y >= options.region_y_min && position.y <= options.region_y_max;
}

/** @brief 根据线性状态方程由密度重建压力。 */
Real pressureFromDensity(Real density, Real rho0, Real sound_speed)
{
    return sound_speed * sound_speed * (density - rho0);
}

/** @brief 从粗网格 restart XML 中读取 WaterBlock 粒子状态。 */
std::vector<CoarseParticle> readCoarseWaterParticles(const fs::path &restart_path)
{
    XMLDocument doc;
    loadXmlDocument(doc, restart_path);
    XMLElement *root = doc.FirstChildElement("restart_data");
    if (root == nullptr)
    {
        throw std::runtime_error("Restart XML has no restart_data root: " + restart_path.string());
    }

    XMLElement *body = root->FirstChildElement("body");
    while (body != nullptr)
    {
        const char *name = body->Attribute("name");
        if (name != nullptr && std::string(name) == "WaterBlock")
        {
            break;
        }
        body = body->NextSiblingElement("body");
    }
    if (body == nullptr)
    {
        throw std::runtime_error("Restart XML has no WaterBlock body: " + restart_path.string());
    }

    std::vector<CoarseParticle> particles;
    for (XMLElement *element = body->FirstChildElement("particle");
         element != nullptr; element = element->NextSiblingElement("particle"))
    {
        CoarseParticle particle;
        particle.original_id = readUnsigned(element, "OriginalID", restart_path);
        particle.volume = readReal(element, "VolumetricMeasure", restart_path);
        particle.density = readReal(element, "Density", restart_path);
        particle.mass = readReal(element, "Mass", restart_path);
        particle.pressure = readReal(element, "Pressure", restart_path);
        particle.position = readVec2(element, "Position", restart_path);
        particle.velocity = readVec2(element, "Velocity", restart_path);
        particle.momentum = readVec2(element, "Momentum", restart_path);
        particles.push_back(particle);
    }
    if (particles.empty())
    {
        throw std::runtime_error("Restart XML WaterBlock body has no particles: " + restart_path.string());
    }
    return particles;
}

/** @brief 从目标 reload XML 中读取粒子位置、体积和 smoothing ratio。 */
std::vector<TargetParticle> readTargetParticles(const fs::path &reload_path)
{
    XMLDocument doc;
    loadXmlDocument(doc, reload_path);
    XMLElement *root = doc.FirstChildElement("particles");
    if (root == nullptr)
    {
        throw std::runtime_error("Reload XML has no particles root: " + reload_path.string());
    }

    std::vector<TargetParticle> particles;
    for (XMLElement *element = root->FirstChildElement("particle");
         element != nullptr; element = element->NextSiblingElement("particle"))
    {
        TargetParticle particle;
        particle.original_id = readUnsigned(element, "OriginalID", reload_path);
        particle.volume = readReal(element, "VolumetricMeasure", reload_path);
        particle.smoothing_length_ratio = readReal(element, "SmoothingLengthRatio", reload_path);
        particle.position = readVec2(element, "Position", reload_path);
        particles.push_back(particle);
    }
    if (particles.empty())
    {
        throw std::runtime_error("Reload XML has no particles: " + reload_path.string());
    }
    return particles;
}

/** @brief 复制文件并覆盖目标路径，用于同步 cylinder reload。 */
void copyFile(const fs::path &source, const fs::path &destination)
{
    std::error_code ec;
    fs::copy_file(source, destination, fs::copy_options::overwrite_existing, ec);
    if (ec)
    {
        throw std::runtime_error("Failed to copy " + source.string() + " to " + destination.string() + ": " + ec.message());
    }
}

/** @brief 根据源粒子最小体积估算粗网格粒子间距。 */
Real estimateSourceSpacing(const std::vector<CoarseParticle> &source)
{
    Real min_volume = std::numeric_limits<Real>::max();
    for (const CoarseParticle &particle : source)
    {
        if (particle.volume > 0.0)
        {
            min_volume = std::min(min_volume, particle.volume);
        }
    }
    if (min_volume == std::numeric_limits<Real>::max())
    {
        throw std::runtime_error("Cannot estimate source spacing from non-positive volumes.");
    }
    return std::sqrt(min_volume);
}

/** @brief 用反距离加权将源粒子的密度、压力和速度插值到目标粒子。 */
void interpolateState(const std::vector<CoarseParticle> &source, const SourceSpatialIndex &source_index,
                      std::vector<TargetParticle> &target,
                      const Options &options)
{
    constexpr Real exact_distance_sq = 1.0e-24;
    for (TargetParticle &particle : target)
    {
        const std::vector<size_t> neighbors = source_index.nearest(particle.position, options.neighbors);
        Real weight_sum = 0.0;
        Real density = 0.0;
        Real pressure = 0.0;
        Vec2 velocity;

        if (!neighbors.empty() &&
            squaredDistance(source[neighbors.front()].position, particle.position) < exact_distance_sq)
        {
            const CoarseParticle &coarse = source[neighbors.front()];
            particle.density = coarse.density;
            particle.pressure = coarse.pressure;
            particle.velocity = coarse.velocity;
            continue;
        }

        for (size_t index : neighbors)
        {
            const CoarseParticle &coarse = source[index];
            const Real d2 = squaredDistance(coarse.position, particle.position);
            const Real weight = Real(1.0) / (std::sqrt(d2) + Real(1.0e-12));
            weight_sum += weight;
            density += weight * coarse.density;
            pressure += weight * coarse.pressure;
            velocity = velocity + weight * coarse.velocity;
        }

        if (!(weight_sum > 0.0))
        {
            throw std::runtime_error("IDW interpolation produced zero weight.");
        }
        particle.density = density / weight_sum;
        particle.pressure = pressure / weight_sum;
        particle.velocity = velocity / weight_sum;
    }
}

/** @brief 统计粗网格源粒子在加密区或非加密区内的质量和动量。 */
BucketStats computeCoarseStats(const std::vector<CoarseParticle> &particles, const Options &options, bool region_only)
{
    BucketStats stats;
    for (const CoarseParticle &particle : particles)
    {
        if (inRegion(options, particle.position) != region_only)
        {
            continue;
        }
        ++stats.count;
        stats.mass += particle.mass;
        stats.momentum = stats.momentum + particle.momentum;
    }
    return stats;
}

/** @brief 统计目标粒子在加密区或非加密区内的质量和动量。 */
BucketStats computeTargetStats(const std::vector<TargetParticle> &particles, const Options &options, bool region_only)
{
    BucketStats stats;
    for (const TargetParticle &particle : particles)
    {
        if (inRegion(options, particle.position) != region_only)
        {
            continue;
        }
        ++stats.count;
        stats.mass += particle.mass;
        stats.momentum = stats.momentum + particle.momentum;
    }
    return stats;
}

/** @brief 根据目标粒子的密度、体积和速度重建质量、动量和缺失压力。 */
void rebuildConservedState(std::vector<TargetParticle> &particles, const Options &options)
{
    for (TargetParticle &particle : particles)
    {
        particle.mass = particle.density * particle.volume;
        particle.momentum = particle.mass * particle.velocity;
        if (!std::isfinite(static_cast<double>(particle.pressure)))
        {
            particle.pressure = pressureFromDensity(particle.density, options.rho0, options.sound_speed);
        }
    }
}

/** @brief 对单个分桶缩放密度并平移速度，使质量和平均动量匹配源状态。 */
void matchBucketConservation(std::vector<TargetParticle> &target, const BucketStats &source_stats,
                             const Options &options, bool region_only)
{
    BucketStats target_stats = computeTargetStats(target, options, region_only);
    if (source_stats.count == 0 || target_stats.count == 0)
    {
        return;
    }
    if (!(source_stats.mass > 0.0) || !(target_stats.mass > 0.0))
    {
        throw std::runtime_error("Cannot conserve remap bucket with non-positive mass.");
    }

    const Real density_scale = source_stats.mass / target_stats.mass;
    for (TargetParticle &particle : target)
    {
        if (inRegion(options, particle.position) == region_only)
        {
            particle.density *= density_scale;
            particle.pressure = pressureFromDensity(particle.density, options.rho0, options.sound_speed);
            particle.mass = particle.density * particle.volume;
            particle.momentum = particle.mass * particle.velocity;
        }
    }

    target_stats = computeTargetStats(target, options, region_only);
    const Vec2 source_mean_velocity = source_stats.momentum / source_stats.mass;
    const Vec2 target_mean_velocity = target_stats.momentum / target_stats.mass;
    const Vec2 velocity_shift = source_mean_velocity - target_mean_velocity;
    for (TargetParticle &particle : target)
    {
        if (inRegion(options, particle.position) == region_only)
        {
            particle.velocity = particle.velocity + velocity_shift;
            particle.momentum = particle.mass * particle.velocity;
        }
    }
}

/** @brief 计算 remap 后相对 remap 前的质量误差。 */
Real relativeMassError(Real after, Real before)
{
    return std::abs(after - before) / std::max(std::abs(before), Real(1.0e-30));
}

/** @brief 按质量-速度尺度归一化计算动量误差。 */
Real relativeMomentumError(const Vec2 &after, const Vec2 &before, Real mass_scale, Real velocity_scale)
{
    const Real scale = std::max(std::max(norm(before), std::abs(mass_scale) * velocity_scale), Real(1.0e-30));
    return norm(after - before) / scale;
}

/** @brief 将 remap 后的 WaterBlock 粒子状态写成 reload XML。 */
void writeWaterReload(const fs::path &path, const std::vector<TargetParticle> &particles)
{
    XMLDocument doc;
    XMLElement *root = doc.NewElement("particles");
    doc.InsertFirstChild(root);
    for (const TargetParticle &particle : particles)
    {
        XMLElement *element = doc.NewElement("particle");
        element->SetAttribute("OriginalID", std::to_string(particle.original_id).c_str());
        element->SetAttribute("VolumetricMeasure", formatReal(particle.volume).c_str());
        element->SetAttribute("SmoothingLengthRatio", formatReal(particle.smoothing_length_ratio).c_str());
        element->SetAttribute("Density", formatReal(particle.density).c_str());
        element->SetAttribute("Mass", formatReal(particle.mass).c_str());
        element->SetAttribute("Pressure", formatReal(particle.pressure).c_str());
        element->SetAttribute("Velocity", formatVec2(particle.velocity).c_str());
        element->SetAttribute("Momentum", formatVec2(particle.momentum).c_str());
        element->SetAttribute("MassChangeRate", formatReal(Real(0.0)).c_str());
        element->SetAttribute("MomentumChangeRate", formatVec2(Vec2{}).c_str());
        element->SetAttribute("Force", formatVec2(Vec2{}).c_str());
        element->SetAttribute("ForcePrior", formatVec2(Vec2{}).c_str());
        element->SetAttribute("PreviousViscousForce", formatVec2(Vec2{}).c_str());
        element->SetAttribute("Position", formatVec2(particle.position).c_str());
        root->InsertEndChild(element);
    }

    const auto result = doc.SaveFile(path.string().c_str());
    if (result != tinyxml2::XML_SUCCESS)
    {
        throw std::runtime_error("Failed to write remapped reload: " + path.string());
    }
}

/** @brief 写出 total、region 和 coarse_non_region 三个分桶的守恒误差报告。 */
void writeReport(const fs::path &path, const BucketStats &source_total, const BucketStats &target_total,
                 const BucketStats &source_region, const BucketStats &target_region,
                 const BucketStats &source_coarse, const BucketStats &target_coarse,
                 const Options &options)
{
    std::ofstream report(path, std::ios::out | std::ios::trunc);
    if (!report.is_open())
    {
        throw std::runtime_error("Failed to open remap report: " + path.string());
    }
    report << "bucket,source_count,target_count,source_mass,target_mass,relative_mass_error,"
              "source_momentum_x,source_momentum_y,target_momentum_x,target_momentum_y,relative_momentum_error\n";

    auto write_row = [&](const std::string &name, const BucketStats &source, const BucketStats &target)
    {
        report << name << ","
               << source.count << "," << target.count << ","
               << std::setprecision(16) << source.mass << "," << target.mass << ","
               << relativeMassError(target.mass, source.mass) << ","
               << source.momentum.x << "," << source.momentum.y << ","
               << target.momentum.x << "," << target.momentum.y << ","
               << relativeMomentumError(target.momentum, source.momentum, source.mass, options.velocity_scale)
               << "\n";
    };

    write_row("total", source_total, target_total);
    write_row("region", source_region, target_region);
    write_row("coarse_non_region", source_coarse, target_coarse);
}

/** @brief 按固定阈值检查 remap 后的质量和动量守恒误差。 */
void enforceThresholds(const BucketStats &source_total, const BucketStats &target_total,
                       const BucketStats &source_region, const BucketStats &target_region,
                       const Options &options)
{
    const Real mass_total = relativeMassError(target_total.mass, source_total.mass);
    const Real mass_region = relativeMassError(target_region.mass, source_region.mass);
    const Real momentum_total =
        relativeMomentumError(target_total.momentum, source_total.momentum, source_total.mass, options.velocity_scale);
    const Real momentum_region =
        relativeMomentumError(target_region.momentum, source_region.momentum, source_region.mass, options.velocity_scale);

    if (mass_total > REMAP_MASS_TOL_TOTAL ||
        mass_region > REMAP_MASS_TOL_REGION ||
        momentum_total > REMAP_MOMENTUM_TOL_TOTAL ||
        momentum_region > REMAP_MOMENTUM_TOL_REGION)
    {
        std::ostringstream message;
        message << "Remap conservation check failed: mass_total=" << mass_total
                << " (tol " << REMAP_MASS_TOL_TOTAL << "), mass_region=" << mass_region
                << " (tol " << REMAP_MASS_TOL_REGION << "), momentum_total=" << momentum_total
                << " (tol " << REMAP_MOMENTUM_TOL_TOTAL << "), momentum_region=" << momentum_region
                << " (tol " << REMAP_MOMENTUM_TOL_REGION << ")";
        throw std::runtime_error(message.str());
    }
}
} // namespace

/** @brief 重映射工具入口：读取源/目标粒子、插值、守恒校正并写出 reload/report。 */
int main(int argc, char *argv[])
{
    try
    {
        const Options options = parseOptions(argc, argv);
        fs::create_directories(options.output_dir);
        if (options.report_path.has_parent_path())
        {
            fs::create_directories(options.report_path.parent_path());
        }

        std::vector<CoarseParticle> source_particles = readCoarseWaterParticles(options.source_restart);
        std::vector<TargetParticle> target_particles = readTargetParticles(options.target_water_reload);

        const Real source_spacing = estimateSourceSpacing(source_particles);
        SourceSpatialIndex source_index(source_particles, source_spacing * Real(2.5));
        interpolateState(source_particles, source_index, target_particles, options);
        rebuildConservedState(target_particles, options);

        const BucketStats source_region = computeCoarseStats(source_particles, options, true);
        const BucketStats source_coarse = computeCoarseStats(source_particles, options, false);
        matchBucketConservation(target_particles, source_region, options, true);
        matchBucketConservation(target_particles, source_coarse, options, false);
        rebuildConservedState(target_particles, options);

        BucketStats source_total = source_region;
        source_total.count += source_coarse.count;
        source_total.mass += source_coarse.mass;
        source_total.momentum = source_total.momentum + source_coarse.momentum;

        const BucketStats target_region = computeTargetStats(target_particles, options, true);
        const BucketStats target_coarse = computeTargetStats(target_particles, options, false);
        BucketStats target_total = target_region;
        target_total.count += target_coarse.count;
        target_total.mass += target_coarse.mass;
        target_total.momentum = target_total.momentum + target_coarse.momentum;

        enforceThresholds(source_total, target_total, source_region, target_region, options);

        writeWaterReload(options.output_dir / "WaterBlock_rld.xml", target_particles);
        copyFile(options.target_cylinder_reload, options.output_dir / "Cylinder_rld.xml");
        writeReport(options.report_path, source_total, target_total, source_region, target_region,
                    source_coarse, target_coarse, options);

        std::cout << "[LocalRefinement][Remap] source_particles=" << source_particles.size()
                  << ", target_particles=" << target_particles.size()
                  << ", report=" << options.report_path.string() << std::endl;
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "[LocalRefinement][Remap][ERROR] " << e.what() << std::endl;
        return 1;
    }
}
