#ifndef CYLINDER_LG_DATA_HPP
#define CYLINDER_LG_DATA_HPP

/**
 * @file cylinder_lg_data.hpp
 * @brief 二维 Eulerian cylinder LG 算例的配置、续算输出辅助函数和运行时参数。
 */

#include "sphinxsys.h"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <functional>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <process.h>
#endif

namespace SPH
{
namespace cylinder_lg
{
namespace fs = std::filesystem;

/**
 * @brief 在构造任何 SPH 对象前，从 config.ini 读取的算例配置。
 */
struct SimulationConfig
{
    std::string config_file;

    Real dl;
    Real dh;
    Real global_resolution;
    Real sponge_width_factor;
    Real cylinder_center_x;
    Real cylinder_center_y;
    Real cylinder_radius;

    Real rho0_f;
    Real u_f;
    Real sound_speed_factor;
    Real re;

    Real end_time;
    Real output_interval;
    int screen_output_interval;
    Real acoustic_cfl;
    Real riemann_limiter_parameter;
    bool run_particle_relaxation;
    bool reload_particles;
    bool pressure_contact_only_normal;
    bool enable_restart;
    int restart_step;
    int restart_output_factor;
    int restart_keep_last_n;
    bool enable_staged_refinement;
    Real coarse_stage_fraction;
    bool enable_local_refinement;
    bool load_remapped_state;
    Real remap_physical_time;
    bool local_refinement_runtime_diagnostics;
    int local_refinement_level;
    Real local_refinement_spacing_factor;
    Real local_refinement_region_x_min;
    Real local_refinement_region_x_max;
    Real local_refinement_region_y_min;
    Real local_refinement_region_y_max;

    /** @brief 返回弱可压流体使用的声速。 */
    Real soundSpeed() const
    {
        return sound_speed_factor * u_f;
    }

    /** @brief 根据密度、来流速度、圆柱直径和雷诺数返回动力黏度。 */
    Real dynamicViscosity() const
    {
        return rho0_f * u_f * (2.0 * cylinder_radius) / re;
    }
};

namespace detail
{
using IniSection = std::unordered_map<std::string, std::string>;
using IniConfig = std::unordered_map<std::string, IniSection>;

/** @brief 返回去除首尾 ASCII 空白字符后的字符串副本。 */
inline std::string trimCopy(const std::string &value)
{
    const auto begin = value.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos)
    {
        return "";
    }
    const auto end = value.find_last_not_of(" \t\r\n");
    return value.substr(begin, end - begin + 1);
}

/** @brief 返回小写字符串副本，用于配置节名和键名的大小写无关匹配。 */
inline std::string toLowerCopy(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

/** @brief 去除 INI 注释，同时保留引号内部的注释符号。 */
inline std::string stripComment(const std::string &line)
{
    bool in_quotes = false;
    for (size_t i = 0; i < line.size(); ++i)
    {
        const char c = line[i];
        if (c == '"')
        {
            in_quotes = !in_quotes;
        }
        if (!in_quotes && (c == '#' || c == ';'))
        {
            return line.substr(0, i);
        }
    }
    return line;
}

/** @brief 去除首行可能存在的 UTF-8 BOM。 */
inline std::string stripUtf8Bom(const std::string &line)
{
    if (line.size() >= 3 &&
        static_cast<unsigned char>(line[0]) == 0xEF &&
        static_cast<unsigned char>(line[1]) == 0xBB &&
        static_cast<unsigned char>(line[2]) == 0xBF)
    {
        return line.substr(3);
    }
    return line;
}

/** @brief 将分节的 config.ini 解析为规范化后的 section/key 映射。 */
inline IniConfig loadIniFile(const fs::path &path)
{
    std::ifstream input(path);
    if (!input.is_open())
    {
        throw std::runtime_error("Failed to open config file: " + path.string());
    }

    IniConfig ini;
    std::string current_section;
    std::string line;
    while (std::getline(input, line))
    {
        line = stripUtf8Bom(line);
        const std::string cleaned = trimCopy(stripComment(line));
        if (cleaned.empty())
        {
            continue;
        }
        if (cleaned.front() == '[' && cleaned.back() == ']')
        {
            current_section = toLowerCopy(trimCopy(cleaned.substr(1, cleaned.size() - 2)));
            continue;
        }

        const size_t equal_pos = cleaned.find('=');
        if (equal_pos == std::string::npos)
        {
            throw std::runtime_error("Invalid config line without '=': " + cleaned);
        }

        const std::string key = toLowerCopy(trimCopy(cleaned.substr(0, equal_pos)));
        const std::string value = trimCopy(cleaned.substr(equal_pos + 1));
        if (current_section.empty())
        {
            throw std::runtime_error("Config key outside any section: " + key);
        }
        if (key.empty())
        {
            throw std::runtime_error("Config contains an empty key in section [" + current_section + "]");
        }
        ini[current_section][key] = value;
    }
    return ini;
}

/** @brief 读取必需的 INI 原始值；缺失时抛出带上下文的错误。 */
inline std::string getIniValue(const IniConfig &ini,
                               const std::string &section,
                               const std::string &key)
{
    const std::string section_key = toLowerCopy(section);
    const std::string item_key = toLowerCopy(key);
    auto section_iter = ini.find(section_key);
    if (section_iter == ini.end())
    {
        throw std::runtime_error("Config missing section [" + section + "]");
    }
    auto item_iter = section_iter->second.find(item_key);
    if (item_iter == section_iter->second.end())
    {
        throw std::runtime_error("Config missing key '" + key + "' in section [" + section + "]");
    }
    return item_iter->second;
}

/** @brief 读取可选的 INI 原始值；缺失时返回给定默认文本。 */
inline std::string getIniValueOrDefault(const IniConfig &ini,
                                        const std::string &section,
                                        const std::string &key,
                                        const std::string &default_value)
{
    const std::string section_key = toLowerCopy(section);
    const std::string item_key = toLowerCopy(key);
    auto section_iter = ini.find(section_key);
    if (section_iter == ini.end())
    {
        return default_value;
    }
    auto item_iter = section_iter->second.find(item_key);
    if (item_iter == section_iter->second.end())
    {
        return default_value;
    }
    return item_iter->second;
}

/** @brief 解析 Real 数值，并拒绝带有尾随字符的输入。 */
inline Real parseReal(const std::string &value, const std::string &name)
{
    std::string text = trimCopy(value);
    size_t parsed = 0;
    try
    {
        const Real result = static_cast<Real>(std::stod(text, &parsed));
        if (parsed != text.size())
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

/** @brief 解析整数数值，并拒绝带有尾随字符的输入。 */
inline int parseInt(const std::string &value, const std::string &name)
{
    std::string text = trimCopy(value);
    size_t parsed = 0;
    try
    {
        const int result = std::stoi(text, &parsed);
        if (parsed != text.size())
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

/** @brief 解析 config.ini 中常见的布尔文本形式。 */
inline bool parseBool(const std::string &value, const std::string &name)
{
    const std::string text = toLowerCopy(trimCopy(value));
    if (text == "true" || text == "1" || text == "yes" || text == "on")
    {
        return true;
    }
    if (text == "false" || text == "0" || text == "no" || text == "off")
    {
        return false;
    }
    throw std::runtime_error("Invalid boolean value for " + name + ": " + value);
}

/** @brief 从指定 section/key 读取必需的 Real 值。 */
inline Real getReal(const IniConfig &ini, const std::string &section, const std::string &key)
{
    return parseReal(getIniValue(ini, section, key), section + "." + key);
}

/** @brief 从指定 section/key 读取必需的整数值。 */
inline int getInt(const IniConfig &ini, const std::string &section, const std::string &key)
{
    return parseInt(getIniValue(ini, section, key), section + "." + key);
}

/** @brief 从指定 section/key 读取必需的布尔值。 */
inline bool getBool(const IniConfig &ini, const std::string &section, const std::string &key)
{
    return parseBool(getIniValue(ini, section, key), section + "." + key);
}

/** @brief 读取可选整数值；缺失时使用编译期默认值。 */
inline int getIntOrDefault(const IniConfig &ini, const std::string &section,
                           const std::string &key, int default_value)
{
    return parseInt(getIniValueOrDefault(ini, section, key, std::to_string(default_value)),
                    section + "." + key);
}

/** @brief 读取可选布尔值；缺失时使用编译期默认值。 */
inline bool getBoolOrDefault(const IniConfig &ini, const std::string &section,
                             const std::string &key, bool default_value)
{
    return parseBool(getIniValueOrDefault(ini, section, key, default_value ? "true" : "false"),
                     section + "." + key);
}

/** @brief 读取可选 Real 值；缺失时使用编译期默认值。 */
inline Real getRealOrDefault(const IniConfig &ini, const std::string &section,
                              const std::string &key, Real default_value)
{
    std::ostringstream default_stream;
    default_stream << std::setprecision(16) << default_value;
    return parseReal(getIniValueOrDefault(ini, section, key, default_stream.str()),
                      section + "." + key);
}

/** @brief 从内部阶段环境变量读取 Real 覆盖值；缺失时保持当前配置值。 */
inline Real getRealEnvOrDefault(const char *name, Real default_value)
{
    if (const char *value = std::getenv(name))
    {
        return parseReal(value, name);
    }
    return default_value;
}

/** @brief 从内部阶段环境变量读取整数覆盖值；缺失时保持当前配置值。 */
inline int getIntEnvOrDefault(const char *name, int default_value)
{
    if (const char *value = std::getenv(name))
    {
        return parseInt(value, name);
    }
    return default_value;
}

/** @brief 从内部阶段环境变量读取布尔覆盖值；缺失时保持当前配置值。 */
inline bool getBoolEnvOrDefault(const char *name, bool default_value)
{
    if (const char *value = std::getenv(name))
    {
        return parseBool(value, name);
    }
    return default_value;
}

/** @brief 校验标量配置值必须为正。 */
inline void requirePositive(Real value, const std::string &name)
{
    if (!(value > 0.0))
    {
        throw std::runtime_error("Config value must be positive: " + name);
    }
}

/** @brief 校验局部加密目标粒子间距因子，避免非法值先进入对数换算。 */
inline void requireSpacingFactorInRange(Real spacing_factor)
{
    if (!std::isfinite(spacing_factor) || !(spacing_factor > Real(0.0)) ||
        spacing_factor > Real(1.0))
    {
        throw std::runtime_error(
            "Config value must be in (0, 1]: localrefinement.spacing_factor");
    }
}

/** @brief 根据目标局部粒子间距因子估算 adaptive cell linked list 所需层级。 */
inline int meshLevelFromSpacingFactor(Real spacing_factor)
{
    requireSpacingFactorInRange(spacing_factor);
    if (spacing_factor >= Real(1.0))
    {
        return 0;
    }
    return static_cast<int>(
        std::ceil(std::log(Real(1.0) / spacing_factor) / std::log(Real(2.0))));
}

/** @brief 在全部配置加载完成后校验跨字段约束。 */
inline void validateConfig(SimulationConfig &cfg)
{
    requirePositive(cfg.dl, "geometry.dl");
    requirePositive(cfg.dh, "geometry.dh");
    requirePositive(cfg.global_resolution, "geometry.global_resolution");
    requirePositive(cfg.sponge_width_factor, "geometry.sponge_width_factor");
    requirePositive(cfg.cylinder_radius, "geometry.cylinder_radius");
    requirePositive(cfg.rho0_f, "fluid.rho0_f");
    requirePositive(cfg.u_f, "fluid.u_f");
    requirePositive(cfg.sound_speed_factor, "fluid.sound_speed_factor");
    requirePositive(cfg.re, "fluid.re");
    requirePositive(cfg.end_time, "simulation.end_time");
    requirePositive(cfg.output_interval, "simulation.output_interval");
    requirePositive(cfg.acoustic_cfl, "simulation.acoustic_cfl");
    requirePositive(cfg.riemann_limiter_parameter, "simulation.riemann_limiter_parameter");
    if (cfg.screen_output_interval <= 0)
    {
        throw std::runtime_error("Config value must be positive: simulation.screen_output_interval");
    }
    if (cfg.restart_step < -1)
    {
        throw std::runtime_error("Config value must be >= -1: restart.restart_step");
    }
    if (cfg.restart_output_factor <= 0)
    {
        throw std::runtime_error("Config value must be positive: restart.restart_output_factor");
    }
    if (cfg.restart_keep_last_n == 0)
    {
        std::cerr << "[WARNING] restart_keep_last_n = 0 is invalid; keeping all checkpoints instead.\n";
        cfg.restart_keep_last_n = -1;
    }
    if (cfg.enable_staged_refinement)
    {
        if (!(cfg.coarse_stage_fraction > 0.0 && cfg.coarse_stage_fraction < 1.0))
        {
            throw std::runtime_error(
                "Config value localrefinement.coarse_stage_fraction must be between 0 and 1.");
        }
    }
    if (cfg.enable_staged_refinement || cfg.enable_local_refinement)
    {
        if (!(cfg.local_refinement_region_x_min < cfg.local_refinement_region_x_max) ||
            !(cfg.local_refinement_region_y_min < cfg.local_refinement_region_y_max))
        {
            throw std::runtime_error("Invalid local refinement region bounds.");
        }
        requireSpacingFactorInRange(cfg.local_refinement_spacing_factor);
    }
    if (cfg.enable_local_refinement)
    {
        if (cfg.local_refinement_level <= 0)
        {
            throw std::runtime_error("Config value must be positive: localrefinement.refinement_level");
        }
        if (cfg.enable_restart && cfg.restart_step != 0)
        {
            throw std::runtime_error(
                "Local refinement cannot read an existing restart checkpoint; set restart_step=0 or enable_restart=false.");
        }
    }
    if (cfg.load_remapped_state)
    {
        if (!cfg.enable_local_refinement)
        {
            throw std::runtime_error(
                "localrefinement.load_remapped_state requires enable_local_refinement=true.");
        }
        if (cfg.run_particle_relaxation || !cfg.reload_particles)
        {
            throw std::runtime_error(
                "localrefinement.load_remapped_state requires run_particle_relaxation=false and reload_particles=true.");
        }
        if (!(cfg.remap_physical_time >= 0.0))
        {
            throw std::runtime_error("Config value must be non-negative: localrefinement.remap_physical_time");
        }
        if (!(cfg.end_time > cfg.remap_physical_time))
        {
            throw std::runtime_error(
                "simulation.end_time must be greater than localrefinement.remap_physical_time for remapped continuation.");
        }
    }
}
} // namespace detail

/** @brief 判断续算模式下的 reduced output 是否应追加写入已有文件。 */
inline bool shouldAppendRestartOutputFile(const fs::path &filepath, bool is_restart)
{
    if (!is_restart)
    {
        return false;
    }

    std::error_code ec;
    if (!fs::exists(filepath, ec) || ec)
    {
        return false;
    }
    if (!fs::is_regular_file(filepath, ec) || ec)
    {
        return false;
    }
    return fs::file_size(filepath, ec) > 0 && !ec;
}

/**
 * @brief 支持续算追加写入的 reduced quantity 记录器，避免重复写表头。
 */
template <class LocalReduceMethodType>
class RestartSafeReducedQuantityRecording : public ReducedQuantityRecording<LocalReduceMethodType>
{
  public:
    template <class DynamicsIdentifier, typename... Args>
    RestartSafeReducedQuantityRecording(bool is_restart, DynamicsIdentifier &identifier, Args &&...args)
        : ReducedQuantityRecording<LocalReduceMethodType>(identifier, std::forward<Args>(args)...)
    {
        if (shouldAppendRestartOutputFile(this->filefullpath_output_, is_restart))
        {
            this->header_written_ = true;
        }
    }
};

/** @brief 在 SPHSystem 清理输出前，临时移走已有 output 文件夹。 */
inline void preserveOutputFolderForRestart(bool is_restart_mode)
{
    if (!is_restart_mode || !fs::exists("./output"))
    {
        return;
    }

    try
    {
        if (fs::exists("./output_temp_preserve"))
        {
            std::cout << "[Restart Mode] Removing stale output_temp_preserve folder..." << std::endl;
            fs::remove_all("./output_temp_preserve");
        }
        fs::rename("./output", "./output_temp_preserve");
        std::cout << "[Restart Mode] Temporarily preserving existing output folder..." << std::endl;
    }
    catch (const fs::filesystem_error &e)
    {
        throw std::runtime_error(std::string("Failed to preserve output folder for restart: ") + e.what());
    }
}

/** @brief 在 SPHSystem 构造完成后恢复被临时保留的 output 文件夹。 */
inline void restoreOutputFolderForRestart(bool is_restart_mode)
{
    if (!is_restart_mode || !fs::exists("./output_temp_preserve"))
    {
        return;
    }

    try
    {
        if (fs::exists("./output"))
        {
            fs::remove_all("./output");
        }
        fs::rename("./output_temp_preserve", "./output");
        std::cout << "[Restart Mode] Restored existing output folder (continuing to write new results)" << std::endl;
    }
    catch (const fs::filesystem_error &e)
    {
        throw std::runtime_error(std::string("Failed to restore output folder for restart: ") + e.what());
    }
}

/** @brief 解析 checkpoint 文件名中固定宽度的 restart step 字段。 */
inline bool tryParseRestartStepToken(const std::string &step_token, int &step)
{
    if (step_token.size() != 10 ||
        !std::all_of(step_token.begin(), step_token.end(),
                     [](unsigned char c) { return std::isdigit(c) != 0; }))
    {
        return false;
    }

    try
    {
        step = std::stoi(step_token);
        return true;
    }
    catch (const std::exception &)
    {
        return false;
    }
}

/** @brief 从固定前缀、10 位数字和固定后缀组成的文件名中提取 restart step。 */
inline bool tryExtractFixedWidthRestartStep(const std::string &filename,
                                            const std::string &prefix,
                                            const std::string &suffix,
                                            int &step)
{
    if (filename.size() != prefix.size() + 10 + suffix.size() ||
        filename.rfind(prefix, 0) != 0 ||
        filename.rfind(suffix) != filename.size() - suffix.size())
    {
        return false;
    }
    return tryParseRestartStepToken(filename.substr(prefix.size(), 10), step);
}

/** @brief 识别主 restart 锚点文件，并返回其 checkpoint step。 */
inline bool tryExtractRestartAnchorStep(const std::string &filename, int &step)
{
    return tryExtractFixedWidthRestartStep(filename, "Restart_time_", ".dat", step) ||
           tryExtractFixedWidthRestartStep(filename, "Restart_time_", ".xml", step) ||
           tryExtractFixedWidthRestartStep(filename, "Restart_", ".xml", step);
}

/** @brief 根据 checkpoint step 构造标准 Restart_XXXXXXXXXX.xml 路径。 */
inline fs::path restartXmlPath(size_t restart_step, const std::string &restart_dir = "restart")
{
    std::ostringstream step_stream;
    step_stream << std::setw(10) << std::setfill('0') << restart_step;
    return fs::path(restart_dir) / ("Restart_" + step_stream.str() + ".xml");
}

/** @brief 打印 restart 兼容性错误并终止，避免继续生成污染结果。 */
[[noreturn]] inline void exitWithRestartError(const std::string &message)
{
    std::cerr << "[RESTART][ERROR] " << message << std::endl;
    std::exit(EXIT_FAILURE);
}

/** @brief 校验 restart XML 是否包含当前可执行程序要求的变量。 */
inline void requireRestartAttributes(size_t restart_step, const StdVec<std::string> &attributes,
                                     const std::string &body_name = "WaterBlock",
                                     const std::string &restart_dir = "restart")
{
    const fs::path restart_path = restartXmlPath(restart_step, restart_dir);
    if (!fs::exists(restart_path))
    {
        return;
    }

    XmlParser restart_xml("xml_restart");
    restart_xml.loadXmlFile(restart_path.string());
    tinyxml2::XMLElement *body_element = restart_xml.first_element_->FirstChildElement("body");
    while (body_element != nullptr)
    {
        const char *name = body_element->Attribute("name");
        if (name != nullptr && body_name == name)
        {
            break;
        }
        body_element = body_element->NextSiblingElement("body");
    }
    if (body_element == nullptr)
    {
        exitWithRestartError("Restart file is missing body '" + body_name + "': " + restart_path.string());
    }

    tinyxml2::XMLElement *particle_element = body_element->FirstChildElement("particle");
    if (particle_element == nullptr)
    {
        exitWithRestartError("Restart file has no particles for body '" + body_name + "': " + restart_path.string());
    }

    for (const auto &attribute : attributes)
    {
        if (particle_element->Attribute(attribute.c_str()) == nullptr)
        {
            exitWithRestartError(
                "Restart file is missing required " + body_name + " attribute '" + attribute +
                "'. Regenerate checkpoints with the current executable before continuing: " +
                restart_path.string());
        }
    }
}

/** @brief 识别属于同一 restart checkpoint 组的辅助文件。 */
inline bool tryExtractRestartRelatedStep(const std::string &filename, int &step)
{
    if (tryExtractRestartAnchorStep(filename, step))
    {
        return true;
    }

    const std::string rst_marker = "_rst_";
    const std::string xml_suffix = ".xml";
    if (filename.size() <= rst_marker.size() + 10 + xml_suffix.size() ||
        filename.rfind(xml_suffix) != filename.size() - xml_suffix.size())
    {
        return false;
    }

    const size_t rst_pos = filename.rfind(rst_marker);
    if (rst_pos == std::string::npos)
    {
        return false;
    }

    const size_t step_start = rst_pos + rst_marker.size();
    const size_t step_length = filename.size() - xml_suffix.size() - step_start;
    if (step_length != 10)
    {
        return false;
    }

    return tryParseRestartStepToken(filename.substr(step_start, step_length), step);
}

/** @brief 在 restart 目录中寻找最新的有效锚点 step。 */
inline int detectLatestRestartStep(const std::string &restart_dir = "restart")
{
    int max_step = 0;
    std::string max_step_filename;
    std::error_code ec;
    if (!fs::exists(restart_dir, ec) || !fs::is_directory(restart_dir, ec))
    {
        std::cout << "Restart folder not found: " << restart_dir << std::endl;
        return 0;
    }

    for (const auto &entry : fs::directory_iterator(restart_dir, ec))
    {
        if (ec || !entry.is_regular_file())
        {
            continue;
        }

        const std::string filename = entry.path().filename().string();
        int step = 0;
        if (tryExtractRestartAnchorStep(filename, step) && step > max_step)
        {
            max_step = step;
            max_step_filename = filename;
        }
    }

    if (max_step > 0)
    {
        std::cout << "Detected latest restart file: " << max_step_filename << std::endl;
    }
    else
    {
        std::cout << "No valid restart file found." << std::endl;
    }
    return max_step;
}

/** @brief 删除旧 checkpoint 组，仅保留最新的 N 个 restart 锚点。 */
inline void cleanupOldCheckpoints(const std::string &restart_dir, int keep_last_n)
{
    if (keep_last_n < 0)
    {
        return;
    }

    std::error_code ec;
    if (!fs::exists(restart_dir, ec) || !fs::is_directory(restart_dir, ec))
    {
        return;
    }

    std::vector<int> checkpoint_steps;
    std::unordered_set<int> checkpoint_step_set;
    std::unordered_map<int, std::vector<fs::path>> step_to_files;

    for (const auto &entry : fs::directory_iterator(restart_dir, ec))
    {
        if (ec || !entry.is_regular_file(ec))
        {
            continue;
        }

        const std::string filename = entry.path().filename().string();
        int step = 0;
        if (tryExtractRestartAnchorStep(filename, step))
        {
            if (checkpoint_step_set.insert(step).second)
            {
                checkpoint_steps.push_back(step);
            }
            step_to_files[step].push_back(entry.path());
        }
        else if (tryExtractRestartRelatedStep(filename, step))
        {
            step_to_files[step].push_back(entry.path());
        }
    }

    if (static_cast<int>(checkpoint_steps.size()) <= keep_last_n)
    {
        return;
    }

    std::sort(checkpoint_steps.begin(), checkpoint_steps.end(), std::greater<int>());
    std::unordered_set<int> steps_to_delete(checkpoint_steps.begin() + keep_last_n,
                                            checkpoint_steps.end());

    size_t deleted_count = 0;
    for (int step : steps_to_delete)
    {
        auto iter = step_to_files.find(step);
        if (iter == step_to_files.end())
        {
            continue;
        }
        for (const auto &file : iter->second)
        {
            ec.clear();
            if (fs::remove(file, ec))
            {
                ++deleted_count;
            }
        }
    }

    if (deleted_count > 0)
    {
        std::cout << "[RESTART] Removed " << steps_to_delete.size()
                  << " old checkpoint group(s), " << deleted_count
                  << " file(s); keeping latest " << keep_last_n << "." << std::endl;
    }
}

/** @brief 从 reduced-output 数据行中解析开头的时间列。 */
inline bool tryParseLeadingRealField(const std::string &line, Real &value)
{
    const std::string trimmed = detail::trimCopy(line);
    if (trimmed.empty())
    {
        return false;
    }

    size_t parsed = 0;
    try
    {
        value = static_cast<Real>(std::stod(trimmed, &parsed));
        return parsed > 0;
    }
    catch (const std::exception &)
    {
        return false;
    }
}

/** @brief 列出续算时需要按 checkpoint 时间截断的已有 reduced-output 文件。 */
inline std::vector<std::string> restartDatFilenamesToTruncate(const std::string &output_folder)
{
    const std::vector<std::string> candidates = {
        "Cylinder_TotalPressureForceFromFluid.dat",
        "Cylinder_TotalViscousForceFromFluid.dat",
        "WaterBlock_MaximumSpeed.dat",
        "fluid_wall_force.dat"};

    std::vector<std::string> existing_files;
    for (const auto &filename : candidates)
    {
        std::error_code ec;
        if (fs::exists(fs::path(output_folder) / filename, ec) && !ec)
        {
            existing_files.push_back(filename);
        }
    }
    return existing_files;
}

/** @brief 将 reduced-output 文件截断到不晚于 restart checkpoint 的数据行。 */
inline void truncateDatFileToTime(const std::string &filepath, Real checkpoint_time, Real eps = 1.0e-6)
{
    std::ifstream input(filepath);
    if (!input.is_open())
    {
        std::cout << "[TRUNCATE] File not found, skip: " << filepath << std::endl;
        return;
    }

    const fs::path source_path(filepath);
    const fs::path temp_path = source_path.string() + ".tmp";
    std::ofstream output(temp_path, std::ios::out | std::ios::trunc);
    if (!output.is_open())
    {
        throw std::runtime_error("Failed to open temp file for truncation: " + temp_path.string());
    }

    std::string line;
    size_t total_rows = 0;
    size_t kept_rows = 0;
    size_t removed_rows = 0;
    while (std::getline(input, line))
    {
        Real time = 0.0;
        if (tryParseLeadingRealField(line, time))
        {
            ++total_rows;
            if (time <= checkpoint_time + eps)
            {
                output << line << "\n";
                ++kept_rows;
            }
            else
            {
                ++removed_rows;
            }
        }
        else
        {
            output << line << "\n";
        }
    }

    input.close();
    output.close();

    std::error_code ec;
    fs::rename(temp_path, source_path, ec);
    if (ec)
    {
        fs::remove(source_path, ec);
        ec.clear();
        fs::rename(temp_path, source_path, ec);
    }
    if (ec)
    {
        throw std::runtime_error("Failed to replace truncated file: " + filepath);
    }

    std::cout << "[TRUNCATE] " << filepath << ": kept " << kept_rows << "/" << total_rows
              << " data rows, removed " << removed_rows << " rows after t="
              << checkpoint_time << "s." << std::endl;
}

/** @brief 判断 VTP 文件名是否属于本算例写出的物体。 */
inline bool hasCylinderLgVtpPrefix(const std::string &filename)
{
    return filename.rfind("WaterBlock_", 0) == 0 || filename.rfind("Cylinder_", 0) == 0;
}

/** @brief 从本算例的 SPHinXsys VTP 文件名中恢复物理时间。 */
inline bool tryParseVtpTimeFromFilename(const fs::path &filepath, Real &time)
{
    const std::string filename = filepath.filename().string();
    if (!hasCylinderLgVtpPrefix(filename))
    {
        return false;
    }

    const size_t underscore_pos = filename.find('_');
    const size_t dot_pos = filename.rfind(".vtp");
    if (underscore_pos == std::string::npos || dot_pos == std::string::npos ||
        dot_pos <= underscore_pos + 1)
    {
        return false;
    }

    const std::string token = filename.substr(underscore_pos + 1, dot_pos - underscore_pos - 1);
    if (token.empty() ||
        !std::all_of(token.begin(), token.end(),
                     [](unsigned char c) { return std::isdigit(c) != 0; }))
    {
        return false;
    }

    try
    {
        time = static_cast<Real>(std::stoull(token)) / 1.0e6;
        return true;
    }
    catch (const std::exception &)
    {
        return false;
    }
}

/** @brief 删除编码时间晚于 restart checkpoint 的可视化帧。 */
inline void cleanupFutureVtpFiles(const std::string &output_folder, Real checkpoint_time,
                                  Real eps = 1.0e-6)
{
    const fs::path output_dir(output_folder);
    std::error_code ec;
    if (!fs::exists(output_dir, ec) || !fs::is_directory(output_dir, ec))
    {
        return;
    }

    int scanned_files = 0;
    int kept_files = 0;
    int removed_files = 0;
    int skipped_files = 0;
    for (const auto &entry : fs::directory_iterator(output_dir, ec))
    {
        if (ec || !entry.is_regular_file(ec))
        {
            continue;
        }
        const fs::path filepath = entry.path();
        if (filepath.extension() != ".vtp" || !hasCylinderLgVtpPrefix(filepath.filename().string()))
        {
            continue;
        }

        ++scanned_files;
        Real frame_time = 0.0;
        if (!tryParseVtpTimeFromFilename(filepath, frame_time))
        {
            ++skipped_files;
            continue;
        }

        if (frame_time > checkpoint_time + eps)
        {
            ec.clear();
            fs::remove(filepath, ec);
            if (ec)
            {
                throw std::runtime_error("Failed to remove future VTP file: " + filepath.string());
            }
            ++removed_files;
        }
        else
        {
            ++kept_files;
        }
    }

    std::cout << "[CLEANUP][VTP] scanned=" << scanned_files
              << ", kept=" << kept_files
              << ", removed=" << removed_files
              << ", skipped=" << skipped_files
              << ", checkpoint_time=" << checkpoint_time << "s" << std::endl;
}

/** @brief 优先从运行目录解析 config.ini，找不到时回退到源码目录。 */
inline fs::path resolveDefaultConfigPath(const fs::path &source_dir = fs::path(__FILE__).parent_path(),
                                         const fs::path &cwd = fs::current_path())
{
    std::error_code ec;
    if (const char *config_from_env = std::getenv("CYLINDER_LG_CONFIG"))
    {
        const fs::path explicit_config = fs::absolute(fs::path(config_from_env), ec);
        const fs::path candidate = ec ? fs::path(config_from_env) : explicit_config;
        ec.clear();
        if (fs::exists(candidate, ec) && !ec)
        {
            return candidate;
        }
        throw std::runtime_error("CYLINDER_LG_CONFIG points to a missing config file: " +
                                 candidate.string());
    }

    const fs::path cwd_config = cwd / "config.ini";
    if (fs::exists(cwd_config, ec) && !ec)
    {
        return cwd_config;
    }

    ec.clear();
    const fs::path source_config = source_dir / "config.ini";
    if (fs::exists(source_config, ec) && !ec)
    {
        return source_config;
    }

    throw std::runtime_error("config.ini not found in current working directory or source directory");
}

/** @brief 从显式 config.ini 路径加载并校验 SimulationConfig。 */
inline SimulationConfig loadSimulationConfig(const fs::path &config_path)
{
    const detail::IniConfig ini = detail::loadIniFile(config_path);

    SimulationConfig cfg;
    cfg.config_file = config_path.string();
    cfg.dl = detail::getReal(ini, "geometry", "dl");
    cfg.dh = detail::getReal(ini, "geometry", "dh");
    cfg.global_resolution = detail::getReal(ini, "geometry", "global_resolution");
    cfg.sponge_width_factor = detail::getReal(ini, "geometry", "sponge_width_factor");
    cfg.cylinder_center_x = detail::getReal(ini, "geometry", "cylinder_center_x");
    cfg.cylinder_center_y = detail::getReal(ini, "geometry", "cylinder_center_y");
    cfg.cylinder_radius = detail::getReal(ini, "geometry", "cylinder_radius");

    cfg.rho0_f = detail::getReal(ini, "fluid", "rho0_f");
    cfg.u_f = detail::getReal(ini, "fluid", "u_f");
    cfg.sound_speed_factor = detail::getReal(ini, "fluid", "sound_speed_factor");
    cfg.re = detail::getReal(ini, "fluid", "re");

    cfg.end_time = detail::getReal(ini, "simulation", "end_time");
    cfg.output_interval = detail::getReal(ini, "simulation", "output_interval");
    cfg.screen_output_interval = detail::getInt(ini, "simulation", "screen_output_interval");
    cfg.acoustic_cfl = detail::getRealOrDefault(ini, "simulation", "acoustic_cfl", 0.5);
    cfg.riemann_limiter_parameter = detail::getReal(ini, "simulation", "riemann_limiter_parameter");
    cfg.run_particle_relaxation = detail::getBool(ini, "simulation", "run_particle_relaxation");
    cfg.reload_particles = detail::getBool(ini, "simulation", "reload_particles");
    cfg.pressure_contact_only_normal = detail::getBool(ini, "simulation", "pressure_contact_only_normal");
    cfg.enable_restart = detail::getBoolOrDefault(ini, "restart", "enable_restart", false);
    cfg.restart_step = detail::getIntOrDefault(ini, "restart", "restart_step", 0);
    cfg.restart_output_factor = detail::getIntOrDefault(ini, "restart", "restart_output_factor", 10);
    cfg.restart_keep_last_n = detail::getIntOrDefault(ini, "restart", "restart_keep_last_n", 10);
    const Real cylinder_diameter = 2.0 * cfg.cylinder_radius;
    cfg.enable_staged_refinement = detail::getBoolOrDefault(
        ini, "localrefinement", "enable_staged_refinement", false);
    cfg.coarse_stage_fraction = detail::getRealOrDefault(
        ini, "localrefinement", "coarse_stage_fraction", 0.5);
    cfg.enable_local_refinement = detail::getBoolOrDefault(ini, "localrefinement", "enable_local_refinement", false);
    cfg.load_remapped_state = detail::getBoolOrDefault(
        ini, "localrefinement", "load_remapped_state", false);
    cfg.remap_physical_time = detail::getRealOrDefault(
        ini, "localrefinement", "remap_physical_time", 0.0);
    cfg.local_refinement_runtime_diagnostics = detail::getBoolOrDefault(
        ini, "localrefinement", "runtime_diagnostics", false);
    cfg.local_refinement_level = detail::getIntOrDefault(ini, "localrefinement", "refinement_level", 1);
    cfg.local_refinement_spacing_factor = detail::getRealOrDefault(
        ini, "localrefinement", "spacing_factor",
        std::pow(Real(0.5), static_cast<Real>(cfg.local_refinement_level)));
    cfg.local_refinement_level = std::max(
        cfg.local_refinement_level,
        detail::meshLevelFromSpacingFactor(cfg.local_refinement_spacing_factor));
    cfg.local_refinement_region_x_min = detail::getRealOrDefault(
        ini, "localrefinement", "region_x_min", cfg.cylinder_center_x);
    cfg.local_refinement_region_x_max = detail::getRealOrDefault(
        ini, "localrefinement", "region_x_max", cfg.cylinder_center_x + 5.0 * cylinder_diameter);
    cfg.local_refinement_region_y_min = detail::getRealOrDefault(
        ini, "localrefinement", "region_y_min", cfg.cylinder_center_y - 3.0 * cylinder_diameter);
    cfg.local_refinement_region_y_max = detail::getRealOrDefault(
        ini, "localrefinement", "region_y_max", cfg.cylinder_center_y + 3.0 * cylinder_diameter);

    cfg.end_time = detail::getRealEnvOrDefault("CYLINDER_LG_STAGE_END_TIME", cfg.end_time);
    cfg.run_particle_relaxation = detail::getBoolEnvOrDefault(
        "CYLINDER_LG_STAGE_RUN_PARTICLE_RELAXATION", cfg.run_particle_relaxation);
    cfg.reload_particles = detail::getBoolEnvOrDefault(
        "CYLINDER_LG_STAGE_RELOAD_PARTICLES", cfg.reload_particles);
    cfg.enable_restart = detail::getBoolEnvOrDefault(
        "CYLINDER_LG_STAGE_ENABLE_RESTART", cfg.enable_restart);
    cfg.restart_step = detail::getIntEnvOrDefault(
        "CYLINDER_LG_STAGE_RESTART_STEP", cfg.restart_step);
    cfg.enable_staged_refinement = detail::getBoolEnvOrDefault(
        "CYLINDER_LG_STAGE_ENABLE_STAGED_REFINEMENT", cfg.enable_staged_refinement);
    cfg.enable_local_refinement = detail::getBoolEnvOrDefault(
        "CYLINDER_LG_STAGE_ENABLE_LOCAL_REFINEMENT", cfg.enable_local_refinement);
    cfg.load_remapped_state = detail::getBoolEnvOrDefault(
        "CYLINDER_LG_STAGE_LOAD_REMAPPED_STATE", cfg.load_remapped_state);
    cfg.remap_physical_time = detail::getRealEnvOrDefault(
        "CYLINDER_LG_STAGE_REMAP_PHYSICAL_TIME", cfg.remap_physical_time);
    cfg.local_refinement_spacing_factor = detail::getRealEnvOrDefault(
        "CYLINDER_LG_STAGE_SPACING_FACTOR", cfg.local_refinement_spacing_factor);
    cfg.local_refinement_level = std::max(
        cfg.local_refinement_level,
        detail::meshLevelFromSpacingFactor(cfg.local_refinement_spacing_factor));

    detail::validateConfig(cfg);
    return cfg;
}

/** @brief 加载当前运行上下文解析到的默认 config.ini。 */
inline SimulationConfig loadSimulationConfig()
{
    return loadSimulationConfig(resolveDefaultConfigPath());
}

/** @brief 本算例加载运行时配置的公开入口。 */
inline SimulationConfig loadCaseConfig()
{
    return loadSimulationConfig();
}

/** @brief 描述自动两段计算中某个子阶段的内存覆盖项。 */
struct StageOverrideSpec
{
    Real end_time = 0.0;
    bool run_particle_relaxation = false;
    bool reload_particles = false;
    bool enable_restart = false;
    int restart_step = 0;
    bool enable_local_refinement = false;
    bool load_remapped_state = false;
    Real remap_physical_time = 0.0;
    Real spacing_factor = 1.0;
    bool force_single_thread = false;
};

/** @brief 记录自动两段计算选中的 coarse restart 文件和物理时间。 */
struct RestartSnapshot
{
    fs::path path;
    Real physical_time = 0.0;
};

/** @brief 用紧凑但足够精确的格式写出 Real 值，供日志和内部阶段环境变量使用。 */
inline std::string formatReal(Real value)
{
    std::ostringstream stream;
    stream << std::setprecision(17) << value;
    return stream.str();
}

/** @brief 将布尔值写成 config.ini 里使用的小写文本。 */
inline std::string formatBool(bool value)
{
    return value ? "true" : "false";
}

/** @brief 判断命令行是否请求只生成阶段计划，不实际启动计算。 */
inline bool hasPlanOnlyArgument(int ac, char *av[])
{
    for (int i = 1; i < ac; ++i)
    {
        const std::string arg = av[i] == nullptr ? "" : av[i];
        if (arg == "-PlanOnly" || arg == "--plan-only")
        {
            return true;
        }
    }
    return false;
}

/** @brief 收集需要转发给子阶段主程序的命令行参数，并剔除父层专用的计划参数。 */
inline std::vector<std::string> collectForwardedSimulationArgs(int ac, char *av[])
{
    std::vector<std::string> args;
    for (int i = 1; i < ac; ++i)
    {
        const std::string arg = av[i] == nullptr ? "" : av[i];
        if (arg == "-PlanOnly" || arg == "--plan-only")
        {
            continue;
        }
        args.push_back(arg);
    }
    return args;
}

/** @brief 以跨平台方式返回绝对路径；失败时保留原始路径以便错误信息可读。 */
inline fs::path absolutePathOrSelf(const fs::path &path)
{
    std::error_code ec;
    const fs::path absolute_path = fs::absolute(path, ec);
    return ec ? path : absolute_path;
}

/** @brief 从 Restart_*.xml 文本中解析 restart_time 属性。 */
inline bool tryReadRestartPhysicalTime(const fs::path &path, Real &physical_time)
{
    std::ifstream input(path);
    if (!input.is_open())
    {
        return false;
    }

    std::string line;
    const std::string marker = "restart_time=\"";
    while (std::getline(input, line))
    {
        const size_t begin = line.find(marker);
        if (begin == std::string::npos)
        {
            continue;
        }
        const size_t value_begin = begin + marker.size();
        const size_t value_end = line.find('"', value_begin);
        if (value_end == std::string::npos)
        {
            return false;
        }
        return tryParseLeadingRealField(line.substr(value_begin, value_end - value_begin),
                                        physical_time);
    }
    return false;
}

/** @brief 在 restart 目录中按物理时间选择最新 coarse checkpoint。 */
inline RestartSnapshot findLatestRestartSnapshot(const fs::path &restart_dir)
{
    std::error_code ec;
    if (!fs::exists(restart_dir, ec) || !fs::is_directory(restart_dir, ec))
    {
        throw std::runtime_error("Restart folder not found: " + restart_dir.string());
    }

    bool found = false;
    RestartSnapshot latest;
    for (const auto &entry : fs::directory_iterator(restart_dir, ec))
    {
        if (ec || !entry.is_regular_file(ec))
        {
            continue;
        }

        const fs::path path = entry.path();
        const std::string filename = path.filename().string();
        int step = 0;
        if (!tryExtractFixedWidthRestartStep(filename, "Restart_", ".xml", step))
        {
            continue;
        }

        Real physical_time = 0.0;
        if (!tryReadRestartPhysicalTime(path, physical_time))
        {
            continue;
        }
        if (!found || physical_time > latest.physical_time)
        {
            found = true;
            latest = {absolutePathOrSelf(path), physical_time};
        }
    }

    if (!found)
    {
        throw std::runtime_error("No valid Restart_*.xml with restart_time found in " +
                                 restart_dir.string());
    }
    return latest;
}

/** @brief 生成用于自动两段计算的唯一内部工作目录名。 */
inline fs::path makeUniqueStagedRunRoot(const fs::path &case_dir)
{
    const fs::path staged_root = case_dir / ".staged_refinement_work";
    fs::create_directories(staged_root);

    const auto now = std::chrono::system_clock::now();
    const std::time_t raw_time = std::chrono::system_clock::to_time_t(now);
    std::tm local_time{};
#ifdef _WIN32
    localtime_s(&local_time, &raw_time);
#else
    localtime_r(&raw_time, &local_time);
#endif

    std::ostringstream stamp_stream;
    stamp_stream << std::put_time(&local_time, "%Y%m%d_%H%M%S");
    const std::string stamp = stamp_stream.str();

    for (int suffix = 0; suffix < 1000; ++suffix)
    {
        std::ostringstream name_stream;
        name_stream << "run_" << stamp;
        if (suffix > 0)
        {
            name_stream << "_" << suffix;
        }
        const fs::path candidate = staged_root / name_stream.str();
        std::error_code ec;
        if (!fs::exists(candidate, ec))
        {
            fs::create_directories(candidate);
            return candidate;
        }
    }
    throw std::runtime_error("Failed to create a unique staged refinement work directory.");
}

/** @brief 临时切换环境变量，保证子进程加载同一 config 并接收阶段覆盖。 */
class ScopedEnvironmentVariable
{
  public:
    ScopedEnvironmentVariable(std::string name, const std::string &value)
        : name_(std::move(name))
    {
        if (const char *old_value = std::getenv(name_.c_str()))
        {
            had_old_value_ = true;
            old_value_ = old_value;
        }
        set(value);
    }

    ~ScopedEnvironmentVariable()
    {
        try
        {
            if (had_old_value_)
            {
                set(old_value_);
            }
            else
            {
                unset();
            }
        }
        catch (const std::exception &)
        {
        }
    }

  private:
    void set(const std::string &value) const
    {
#ifdef _WIN32
        if (_putenv_s(name_.c_str(), value.c_str()) != 0)
        {
            throw std::runtime_error("Failed to set environment variable: " + name_);
        }
#else
        if (setenv(name_.c_str(), value.c_str(), 1) != 0)
        {
            throw std::runtime_error("Failed to set environment variable: " + name_);
        }
#endif
    }

    void unset() const
    {
#ifdef _WIN32
        _putenv_s(name_.c_str(), "");
#else
        unsetenv(name_.c_str());
#endif
    }

    std::string name_;
    bool had_old_value_ = false;
    std::string old_value_;
};

/** @brief 批量临时设置环境变量，析构时自动恢复。 */
class ScopedEnvironmentVariables
{
  public:
    explicit ScopedEnvironmentVariables(const std::vector<std::pair<std::string, std::string>> &values)
    {
        variables_.reserve(values.size());
        for (const auto &item : values)
        {
            variables_.push_back(std::make_unique<ScopedEnvironmentVariable>(item.first, item.second));
        }
    }

  private:
    std::vector<std::unique_ptr<ScopedEnvironmentVariable>> variables_;
};

/** @brief 临时切换当前工作目录，让正式两段计算保持在 case 根目录输出。 */
class ScopedCurrentPath
{
  public:
    explicit ScopedCurrentPath(const fs::path &path)
        : previous_path_(fs::current_path())
    {
        fs::current_path(path);
    }

    ~ScopedCurrentPath()
    {
        try
        {
            fs::current_path(previous_path_);
        }
        catch (const fs::filesystem_error &)
        {
        }
    }

  private:
    fs::path previous_path_;
};

/** @brief 对命令行参数加引号，避免路径空格影响 std::system 调用。 */
inline std::string quoteCommandArgument(const std::string &arg)
{
#ifdef _WIN32
    std::string quoted = "\"";
    size_t backslashes = 0;
    for (char c : arg)
    {
        if (c == '\\')
        {
            ++backslashes;
        }
        else if (c == '"')
        {
            quoted.append(backslashes * 2 + 1, '\\');
            quoted.push_back('"');
            backslashes = 0;
        }
        else
        {
            quoted.append(backslashes, '\\');
            quoted.push_back(c);
            backslashes = 0;
        }
    }
    quoted.append(backslashes * 2, '\\');
    quoted.push_back('"');
    return quoted;
#else
    std::string quoted = "'";
    for (char c : arg)
    {
        if (c == '\'')
        {
            quoted += "'\\''";
        }
        else
        {
            quoted.push_back(c);
        }
    }
    quoted.push_back('\'');
    return quoted;
#endif
}

/** @brief 拼出一个完整命令行。 */
inline std::string buildCommandLine(const fs::path &executable,
                                    const std::vector<std::string> &args)
{
    std::string command = quoteCommandArgument(executable.string());
    for (const auto &arg : args)
    {
        command += " ";
        command += quoteCommandArgument(arg);
    }
    return command;
}

/** @brief 启动外部进程；Windows 下直接 spawn，避免 cmd.exe 重新解释路径和引号。 */
inline int runExternalProcess(const fs::path &executable,
                              const std::vector<std::string> &args)
{
#ifdef _WIN32
    std::vector<std::string> argv_storage;
    argv_storage.reserve(args.size() + 1);
    argv_storage.push_back(executable.string());
    argv_storage.insert(argv_storage.end(), args.begin(), args.end());

    std::vector<const char *> argv_ptrs;
    argv_ptrs.reserve(argv_storage.size() + 1);
    for (const auto &arg : argv_storage)
    {
        argv_ptrs.push_back(arg.c_str());
    }
    argv_ptrs.push_back(nullptr);

    const intptr_t result = _spawnv(_P_WAIT, executable.string().c_str(), argv_ptrs.data());
    if (result == -1)
    {
        throw std::runtime_error("Failed to start process: " + executable.string() +
                                 ", errno=" + std::to_string(errno));
    }
    return static_cast<int>(result);
#else
    return std::system(buildCommandLine(executable, args).c_str());
#endif
}

/** @brief 生成内部阶段覆盖环境变量，子进程仍读取同一个用户 config.ini。 */
inline std::vector<std::pair<std::string, std::string>>
makeStageEnvironmentOverrides(const StageOverrideSpec &stage)
{
    std::vector<std::pair<std::string, std::string>> env_values = {
        {"CYLINDER_LG_STAGE_END_TIME", formatReal(stage.end_time)},
        {"CYLINDER_LG_STAGE_RUN_PARTICLE_RELAXATION", formatBool(stage.run_particle_relaxation)},
        {"CYLINDER_LG_STAGE_RELOAD_PARTICLES", formatBool(stage.reload_particles)},
        {"CYLINDER_LG_STAGE_ENABLE_RESTART", formatBool(stage.enable_restart)},
        {"CYLINDER_LG_STAGE_RESTART_STEP", std::to_string(stage.restart_step)},
        {"CYLINDER_LG_STAGE_ENABLE_STAGED_REFINEMENT", "false"},
        {"CYLINDER_LG_STAGE_ENABLE_LOCAL_REFINEMENT", formatBool(stage.enable_local_refinement)},
        {"CYLINDER_LG_STAGE_LOAD_REMAPPED_STATE", formatBool(stage.load_remapped_state)},
        {"CYLINDER_LG_STAGE_REMAP_PHYSICAL_TIME", formatReal(stage.remap_physical_time)},
        {"CYLINDER_LG_STAGE_SPACING_FACTOR", formatReal(stage.spacing_factor)}};
    return env_values;
}

/** @brief 在指定工作目录和阶段覆盖下运行一个子进程，并检查退出码。 */
inline void invokeCheckedStageProcess(const fs::path &executable,
                                      const std::vector<std::string> &args,
                                      const fs::path &working_dir,
                                      const std::string &description,
                                      const fs::path &config_path,
                                      const StageOverrideSpec &stage)
{
    std::cout << "=== " << description << " ===" << std::endl;
    std::cout << "[StagedRefinement] Using config: "
              << absolutePathOrSelf(config_path).string() << std::endl;
    std::cout << "[StagedRefinement] Stage overrides: end_time="
              << formatReal(stage.end_time)
              << ", local_refinement=" << formatBool(stage.enable_local_refinement)
              << ", remapped_reload=" << formatBool(stage.load_remapped_state)
              << std::endl;

    std::vector<std::pair<std::string, std::string>> env_values =
        makeStageEnvironmentOverrides(stage);
    env_values.push_back({"CYLINDER_LG_CONFIG", absolutePathOrSelf(config_path).string()});
    ScopedEnvironmentVariables scoped_env(env_values);
    ScopedCurrentPath scoped_path(working_dir);
    const int exit_code = runExternalProcess(executable, args);
    if (exit_code != 0)
    {
        throw std::runtime_error(description + " failed with exit code " +
                                 std::to_string(exit_code));
    }
}

/** @brief 运行非主计算辅助进程，并检查退出码。 */
inline void invokeCheckedUtilityProcess(const fs::path &executable,
                                        const std::vector<std::string> &args,
                                        const fs::path &working_dir,
                                        const std::string &description)
{
    std::cout << "=== " << description << " ===" << std::endl;
    ScopedCurrentPath scoped_path(working_dir);
    const int exit_code = runExternalProcess(executable, args);
    if (exit_code != 0)
    {
        throw std::runtime_error(description + " failed with exit code " +
                                 std::to_string(exit_code));
    }
}

/** @brief 在同一 case 根目录内执行 coarse->remap->refined 自动两段计算。 */
inline int runStagedRefinementWorkflow(const SimulationConfig &base_config, int ac, char *av[])
{
    const fs::path case_dir = absolutePathOrSelf(fs::current_path());
    const bool plan_only = hasPlanOnlyArgument(ac, av);
    const std::vector<std::string> simulation_args = collectForwardedSimulationArgs(ac, av);
    const fs::path main_executable = absolutePathOrSelf(fs::path(av[0]));
    const fs::path remap_executable =
        main_executable.parent_path() /
        (main_executable.stem().string() + "_remap" + main_executable.extension().string());

    if (!plan_only && !fs::exists(main_executable))
    {
        throw std::runtime_error("Main executable not found: " + main_executable.string());
    }
    if (!plan_only && !fs::exists(remap_executable))
    {
        throw std::runtime_error("Remap executable not found: " + remap_executable.string());
    }

    const Real final_time = base_config.end_time;
    const Real coarse_end_time = final_time * base_config.coarse_stage_fraction;
    const fs::path output_dir = case_dir / "output";
    const fs::path user_config_path = absolutePathOrSelf(base_config.config_file);

    StageOverrideSpec coarse_stage;
    coarse_stage.end_time = coarse_end_time;
    coarse_stage.run_particle_relaxation = false;
    coarse_stage.reload_particles = false;
    coarse_stage.enable_restart = true;
    coarse_stage.restart_step = 0;
    coarse_stage.enable_local_refinement = false;
    coarse_stage.load_remapped_state = false;
    coarse_stage.remap_physical_time = 0.0;
    coarse_stage.spacing_factor = base_config.local_refinement_spacing_factor;

    StageOverrideSpec geometry_stage;
    geometry_stage.end_time = 0.1;
    geometry_stage.run_particle_relaxation = true;
    geometry_stage.reload_particles = false;
    geometry_stage.enable_restart = false;
    geometry_stage.restart_step = 0;
    geometry_stage.enable_local_refinement = true;
    geometry_stage.load_remapped_state = false;
    geometry_stage.remap_physical_time = 0.0;
    geometry_stage.spacing_factor = base_config.local_refinement_spacing_factor;
    geometry_stage.force_single_thread = true;

    StageOverrideSpec refined_stage;
    refined_stage.end_time = final_time;
    refined_stage.run_particle_relaxation = false;
    refined_stage.reload_particles = true;
    refined_stage.enable_restart = true;
    refined_stage.restart_step = 0;
    refined_stage.enable_local_refinement = true;
    refined_stage.load_remapped_state = true;
    refined_stage.remap_physical_time = coarse_end_time;
    refined_stage.spacing_factor = base_config.local_refinement_spacing_factor;

    if (plan_only)
    {
        std::cout << "[StagedRefinement] Plan only; no stage files were created." << std::endl;
        std::cout << "[StagedRefinement] simulation output folder: "
                  << output_dir.string() << std::endl;
        std::cout << "[StagedRefinement] temporary stage files are created under .staged_refinement_work during a real run and removed after success."
                  << std::endl;
        std::cout << "[StagedRefinement] coarse fraction="
                  << formatReal(base_config.coarse_stage_fraction)
                  << ", coarse end_time=" << formatReal(coarse_end_time)
                  << ", refined final end_time=" << formatReal(final_time)
                  << ", spacing_factor="
                  << formatReal(base_config.local_refinement_spacing_factor)
                  << std::endl;
        return 0;
    }

    const fs::path run_root = makeUniqueStagedRunRoot(case_dir);
    const fs::path geometry_dir = run_root / "relaxed_geometry";

    invokeCheckedStageProcess(main_executable, simulation_args, case_dir,
                              "Stage 1: coarse run to " + formatReal(coarse_end_time) + "s",
                              user_config_path, coarse_stage);

    const RestartSnapshot source_restart = findLatestRestartSnapshot(case_dir / "restart");
    std::cout << "[StagedRefinement] Source restart: "
              << source_restart.path.string() << std::endl;
    std::cout << "[StagedRefinement] Source restart physical time: "
              << formatReal(source_restart.physical_time) << "s" << std::endl;

    fs::create_directories(geometry_dir);
    {
        std::error_code ec;
        const fs::path preserved_restart_dir = run_root / "coarse_restart";
        fs::create_directories(preserved_restart_dir, ec);
        ec.clear();
        const fs::path preserved_restart = preserved_restart_dir / source_restart.path.filename();
        fs::copy_file(source_restart.path, preserved_restart,
                      fs::copy_options::overwrite_existing, ec);
        if (ec)
        {
            std::cerr << "[StagedRefinement][WARNING] Failed to preserve source restart: "
                      << ec.message() << std::endl;
        }
        else
        {
            std::cout << "[StagedRefinement] Preserved source restart copy: "
                      << preserved_restart.string() << std::endl;
        }
    }
    refined_stage.remap_physical_time = source_restart.physical_time;

    invokeCheckedStageProcess(main_executable, {}, geometry_dir,
                              "Preparing relaxed refined geometry",
                              user_config_path, geometry_stage);

    const fs::path target_water_reload = geometry_dir / "reload" / "WaterBlock_rld.xml";
    const fs::path target_cylinder_reload = geometry_dir / "reload" / "Cylinder_rld.xml";
    if (!fs::exists(target_water_reload))
    {
        throw std::runtime_error("Missing relaxed target WaterBlock reload: " +
                                 target_water_reload.string());
    }
    if (!fs::exists(target_cylinder_reload))
    {
        throw std::runtime_error("Missing relaxed target Cylinder reload: " +
                                 target_cylinder_reload.string());
    }

    const std::vector<std::string> remap_args = {
        "--source-restart", source_restart.path.string(),
        "--target-water-reload", absolutePathOrSelf(target_water_reload).string(),
        "--target-cylinder-reload", absolutePathOrSelf(target_cylinder_reload).string(),
        "--output-dir", "reload",
        "--report", absolutePathOrSelf(output_dir / "remap_report.csv").string(),
        "--rho0", formatReal(base_config.rho0_f),
        "--sound-speed", formatReal(base_config.soundSpeed()),
        "--velocity-scale", formatReal(base_config.u_f),
        "--region-x-min", formatReal(base_config.local_refinement_region_x_min),
        "--region-x-max", formatReal(base_config.local_refinement_region_x_max),
        "--region-y-min", formatReal(base_config.local_refinement_region_y_min),
        "--region-y-max", formatReal(base_config.local_refinement_region_y_max)};
    invokeCheckedUtilityProcess(remap_executable, remap_args, case_dir,
                                "Remapping coarse restart to refined reload");

    {
        std::error_code ec;
        const fs::path geometry_output_dir = geometry_dir / "output";
        if (fs::exists(geometry_output_dir, ec) && !ec)
        {
            const auto removed_count = fs::remove_all(geometry_output_dir, ec);
            if (!ec)
            {
                std::cout << "[StagedRefinement] Removed temporary geometry output folder: "
                          << geometry_output_dir.string() << " (" << removed_count
                          << " item(s))" << std::endl;
            }
        }
    }

    invokeCheckedStageProcess(main_executable, simulation_args, case_dir,
                              "Stage 2: refined run to " + formatReal(final_time) + "s",
                              user_config_path, refined_stage);

    {
        std::error_code ec;
        if (fs::exists(run_root, ec) && !ec)
        {
            const auto removed_count = fs::remove_all(run_root, ec);
            if (!ec)
            {
                std::cout << "[StagedRefinement] Removed temporary stage work: "
                          << run_root.string() << " (" << removed_count
                          << " item(s))" << std::endl;
            }
        }
        fs::remove(run_root.parent_path(), ec);
    }

    std::cout << "[StagedRefinement] Completed." << std::endl;
    std::cout << "[StagedRefinement] Simulation output: "
              << output_dir.string() << std::endl;
    return 0;
}
} // namespace cylinder_lg

//----------------------------------------------------------------------
//	Case parameters populated from config.ini before SPHSystem construction.
//----------------------------------------------------------------------
inline Real DL;
inline Real DH;
inline Real global_resolution;
inline Real DL_sponge;
inline Real DH_sponge;
inline Vec2d cylinder_center = Vec2d::Zero();
inline Real cylinder_radius;
inline Real rho0_f;
inline Real U_f;
inline Real c_f;
inline Real Re;
inline Real mu_f;
inline bool enable_local_refinement;
inline bool load_remapped_state;
inline Real remap_physical_time;
inline bool local_refinement_runtime_diagnostics;
inline int local_refinement_level;
inline Real local_refinement_spacing_factor;
inline Real local_refinement_region_x_min;
inline Real local_refinement_region_x_max;
inline Real local_refinement_region_y_min;
inline Real local_refinement_region_y_max;
inline Real acoustic_cfl;

/** @brief 将配置值发布到 SPHinXsys 初始化仍使用的 case-level 全局参数。 */
inline void applyCaseConfig(const cylinder_lg::SimulationConfig &config)
{
    DL = config.dl;
    DH = config.dh;
    global_resolution = config.global_resolution;
    DL_sponge = config.sponge_width_factor * global_resolution;
    DH_sponge = config.sponge_width_factor * global_resolution;
    cylinder_center = Vec2d(config.cylinder_center_x, config.cylinder_center_y);
    cylinder_radius = config.cylinder_radius;

    rho0_f = config.rho0_f;
    U_f = config.u_f;
    c_f = config.soundSpeed();
    Re = config.re;
    mu_f = config.dynamicViscosity();
    acoustic_cfl = config.acoustic_cfl;

    enable_local_refinement = config.enable_local_refinement;
    load_remapped_state = config.load_remapped_state;
    remap_physical_time = config.remap_physical_time;
    local_refinement_runtime_diagnostics = config.local_refinement_runtime_diagnostics;
    local_refinement_level = config.local_refinement_level;
    local_refinement_spacing_factor = config.local_refinement_spacing_factor;
    local_refinement_region_x_min = config.local_refinement_region_x_min;
    local_refinement_region_x_max = config.local_refinement_region_x_max;
    local_refinement_region_y_min = config.local_refinement_region_y_min;
    local_refinement_region_y_max = config.local_refinement_region_y_max;
}
} // namespace SPH

#endif // CYLINDER_LG_DATA_HPP
