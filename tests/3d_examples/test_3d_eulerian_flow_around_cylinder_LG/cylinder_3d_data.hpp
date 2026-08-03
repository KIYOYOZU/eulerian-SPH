#ifndef CYLINDER_3D_DATA_HPP
#define CYLINDER_3D_DATA_HPP

#include "sphinxsys.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace SPH
{
namespace cylinder_3d
{

using IniConfig = std::map<std::string, std::map<std::string, std::string>>;

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

inline std::string toLowerCopy(const std::string &value)
{
    std::string out = value;
    std::transform(out.begin(), out.end(), out.begin(),
                   [](unsigned char ch)
                   { return static_cast<char>(std::tolower(ch)); });
    return out;
}

inline std::string stripComment(const std::string &value)
{
    const size_t hash_pos = value.find('#');
    const size_t semi_pos = value.find(';');
    size_t cut = std::string::npos;
    if (hash_pos != std::string::npos)
    {
        cut = hash_pos;
    }
    if (semi_pos != std::string::npos)
    {
        cut = cut == std::string::npos ? semi_pos : std::min(cut, semi_pos);
    }
    return cut == std::string::npos ? value : value.substr(0, cut);
}

inline IniConfig loadIniFile(const std::string &path)
{
    std::ifstream in(path, std::ios::in | std::ios::binary);
    if (!in.is_open())
    {
        throw std::runtime_error("Cannot open config file: " + path);
    }

    char bom[3] = {};
    in.read(bom, 3);
    if (!(static_cast<unsigned char>(bom[0]) == 0xEF &&
          static_cast<unsigned char>(bom[1]) == 0xBB &&
          static_cast<unsigned char>(bom[2]) == 0xBF))
    {
        in.clear();
        in.seekg(0);
    }

    IniConfig ini;
    std::string section;
    std::string line;
    while (std::getline(in, line))
    {
        std::string cleaned = trimCopy(stripComment(line));
        if (cleaned.empty())
        {
            continue;
        }
        if (cleaned.front() == '[' && cleaned.back() == ']')
        {
            section = toLowerCopy(trimCopy(cleaned.substr(1, cleaned.size() - 2)));
            continue;
        }
        const size_t eq_pos = cleaned.find('=');
        if (eq_pos == std::string::npos || section.empty())
        {
            continue;
        }
        const std::string key = toLowerCopy(trimCopy(cleaned.substr(0, eq_pos)));
        const std::string value = trimCopy(cleaned.substr(eq_pos + 1));
        if (!key.empty())
        {
            ini[section][key] = value;
        }
    }
    return ini;
}

inline Real parseReal(const std::string &value)
{
    try
    {
        return static_cast<Real>(std::stod(trimCopy(value)));
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error("Failed to parse Real value '" + value + "': " + e.what());
    }
}

inline int parseInt(const std::string &value)
{
    try
    {
        return std::stoi(trimCopy(value));
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error("Failed to parse int value '" + value + "': " + e.what());
    }
}

inline bool parseBool(const std::string &value)
{
    const std::string normalized = toLowerCopy(trimCopy(value));
    if (normalized == "true" || normalized == "1" || normalized == "yes" || normalized == "on")
    {
        return true;
    }
    if (normalized == "false" || normalized == "0" || normalized == "no" || normalized == "off")
    {
        return false;
    }
    throw std::runtime_error("Failed to parse bool value '" + value + "'.");
}

inline std::string requireValue(const IniConfig &ini, const std::string &section, const std::string &key)
{
    const std::string section_key = toLowerCopy(section);
    const std::string item_key = toLowerCopy(key);
    auto section_it = ini.find(section_key);
    if (section_it == ini.end())
    {
        throw std::runtime_error("Config missing section: [" + section + "]");
    }
    auto item_it = section_it->second.find(item_key);
    if (item_it == section_it->second.end())
    {
        throw std::runtime_error("Config missing key '" + key + "' in section [" + section + "]");
    }
    return item_it->second;
}

inline Real optionalReal(const IniConfig &ini, const std::string &section, const std::string &key,
                         Real default_value)
{
    auto section_it = ini.find(toLowerCopy(section));
    if (section_it == ini.end())
    {
        return default_value;
    }
    auto item_it = section_it->second.find(toLowerCopy(key));
    return item_it == section_it->second.end() ? default_value : parseReal(item_it->second);
}

inline int optionalInt(const IniConfig &ini, const std::string &section, const std::string &key,
                       int default_value)
{
    auto section_it = ini.find(toLowerCopy(section));
    if (section_it == ini.end())
    {
        return default_value;
    }
    auto item_it = section_it->second.find(toLowerCopy(key));
    return item_it == section_it->second.end() ? default_value : parseInt(item_it->second);
}

inline bool optionalBool(const IniConfig &ini, const std::string &section, const std::string &key,
                         bool default_value)
{
    auto section_it = ini.find(toLowerCopy(section));
    if (section_it == ini.end())
    {
        return default_value;
    }
    auto item_it = section_it->second.find(toLowerCopy(key));
    return item_it == section_it->second.end() ? default_value : parseBool(item_it->second);
}

struct Cylinder3DConfig
{
    Real DL = 0.0;
    Real DH = 0.0;
    Real DW = 0.0;
    Real global_resolution = 0.0;
    Real sponge_width_factor = 3.0;
    Real cylinder_center_x = 0.0;
    Real cylinder_center_y = 0.0;
    Real cylinder_radius = 0.0;
    int max_local_particles = 1000000;

    Real rho0_f = 0.0;
    Real u_f = 0.0;
    Real sound_speed_factor = 0.0;
    Real re = 0.0;
    Real outlet_pressure = 0.0;

    Real end_time = 0.0;
    int output_interval = 1;
    int screen_output_interval = 10;
    Real acoustic_cfl = 0.25;
    int boundary_n_layers = 3;
    bool enable_restart = false;
    int restart_step = 0;
    int restart_output_factor = 10;
    int restart_keep_last_n = 3;

    int smoke_min_steps = 50;
    Real cl_cd_epsilon = 1.0e-12;

    Real c_f = 0.0;
    Real mu_f = 0.0;
    Real dp = 0.0;
    Real D = 0.0;
    Real sponge_width = 0.0;
};

inline void deriveAndValidate(Cylinder3DConfig &cfg)
{
    cfg.dp = cfg.global_resolution;
    cfg.D = 2.0 * cfg.cylinder_radius;
    cfg.c_f = cfg.sound_speed_factor * cfg.u_f;
    cfg.mu_f = cfg.rho0_f * cfg.u_f * cfg.D / cfg.re;
    cfg.sponge_width = cfg.sponge_width_factor * cfg.dp;

    if (cfg.DL <= 0.0 || cfg.DH <= 0.0 || cfg.DW <= 0.0 || cfg.dp <= 0.0)
    {
        throw std::runtime_error("Invalid geometry: dl, dh, dw and global_resolution must be positive.");
    }
    if (cfg.sponge_width_factor <= 0.0)
    {
        throw std::runtime_error("Invalid geometry: sponge_width_factor must be positive.");
    }
    if (cfg.cylinder_radius <= 0.0)
    {
        throw std::runtime_error("Invalid geometry: cylinder_radius must be positive.");
    }
    const Real clearance_x = std::min(cfg.cylinder_center_x, cfg.DL - cfg.cylinder_center_x) - cfg.cylinder_radius;
    const Real clearance_y = std::min(cfg.cylinder_center_y, cfg.DH - cfg.cylinder_center_y) - cfg.cylinder_radius;
    if (clearance_x <= cfg.dp || clearance_y <= cfg.dp)
    {
        throw std::runtime_error("Invalid geometry: cylinder must stay inside the physical box with at least one dp clearance.");
    }
    if (cfg.rho0_f <= 0.0 || cfg.u_f <= 0.0 || cfg.sound_speed_factor <= 0.0 || cfg.re <= 0.0)
    {
        throw std::runtime_error("Invalid physical parameters: rho0_f, u_f, sound_speed_factor and re must be positive.");
    }
    if (std::fabs(cfg.outlet_pressure) > 1.0e-14)
    {
        throw std::runtime_error("Unsupported outlet_pressure: first 3D LG cylinder version only supports p_out = 0 gauge.");
    }
    if (cfg.c_f < 10.0 * cfg.u_f)
    {
        throw std::runtime_error("Low-Mach guard violated: sound_speed_factor * u_f must be at least 10 * u_f.");
    }
    if (cfg.end_time <= 0.0 || cfg.output_interval <= 0 || cfg.screen_output_interval <= 0 ||
        cfg.acoustic_cfl <= 0.0 || cfg.boundary_n_layers <= 0 || cfg.smoke_min_steps <= 0)
    {
        throw std::runtime_error("Invalid simulation controls: time, intervals, CFL and layer counts must be positive.");
    }
    if (!cfg.enable_restart && cfg.restart_step != 0)
    {
        throw std::runtime_error("restart_step must be 0 when enable_restart is false.");
    }
    if (cfg.restart_step < -1 || cfg.restart_output_factor <= 0)
    {
        throw std::runtime_error("Invalid restart controls: restart_step must be -1/0/>0 and restart_output_factor must be positive.");
    }
    if (cfg.restart_keep_last_n == 0)
    {
        cfg.restart_keep_last_n = -1;
    }

    const Real nx = (cfg.DL + 2.0 * cfg.sponge_width) / cfg.dp;
    const Real ny = (cfg.DH + 2.0 * cfg.sponge_width) / cfg.dp;
    const Real nz = cfg.DW / cfg.dp;
    const Real estimated_particles = nx * ny * nz;
    if (estimated_particles > static_cast<Real>(cfg.max_local_particles))
    {
        throw std::runtime_error("Estimated fluid particle count exceeds max_local_particles; refusing local large case.");
    }
}

inline Cylinder3DConfig loadConfig(const std::string &config_path)
{
    IniConfig ini = loadIniFile(config_path);
    Cylinder3DConfig cfg;

    cfg.DL = parseReal(requireValue(ini, "geometry", "dl"));
    cfg.DH = parseReal(requireValue(ini, "geometry", "dh"));
    cfg.DW = parseReal(requireValue(ini, "geometry", "dw"));
    cfg.global_resolution = parseReal(requireValue(ini, "geometry", "global_resolution"));
    cfg.sponge_width_factor = optionalReal(ini, "geometry", "sponge_width_factor", 3.0);
    cfg.cylinder_center_x = parseReal(requireValue(ini, "geometry", "cylinder_center_x"));
    cfg.cylinder_center_y = parseReal(requireValue(ini, "geometry", "cylinder_center_y"));
    cfg.cylinder_radius = parseReal(requireValue(ini, "geometry", "cylinder_radius"));
    cfg.max_local_particles = optionalInt(ini, "geometry", "max_local_particles", 1000000);

    cfg.rho0_f = parseReal(requireValue(ini, "physical", "rho0_f"));
    cfg.u_f = parseReal(requireValue(ini, "physical", "u_f"));
    cfg.sound_speed_factor = parseReal(requireValue(ini, "physical", "sound_speed_factor"));
    cfg.re = parseReal(requireValue(ini, "physical", "re"));
    cfg.outlet_pressure = optionalReal(ini, "physical", "outlet_pressure", 0.0);

    cfg.end_time = parseReal(requireValue(ini, "simulation", "end_time"));
    cfg.output_interval = parseInt(requireValue(ini, "simulation", "output_interval"));
    cfg.screen_output_interval = parseInt(requireValue(ini, "simulation", "screen_output_interval"));
    cfg.acoustic_cfl = parseReal(requireValue(ini, "simulation", "acoustic_cfl"));
    cfg.boundary_n_layers = parseInt(requireValue(ini, "simulation", "boundary_n_layers"));
    cfg.enable_restart = optionalBool(ini, "simulation", "enable_restart", false);
    cfg.restart_step = optionalInt(ini, "simulation", "restart_step", 0);
    cfg.restart_output_factor = optionalInt(ini, "simulation", "restart_output_factor", 10);
    cfg.restart_keep_last_n = optionalInt(ini, "simulation", "restart_keep_last_n", 3);

    cfg.smoke_min_steps = optionalInt(ini, "validation", "smoke_min_steps", 50);
    cfg.cl_cd_epsilon = optionalReal(ini, "validation", "cl_cd_epsilon", 1.0e-12);

    deriveAndValidate(cfg);
    return cfg;
}

inline bool tryParseRestartStepToken(const std::string &step_token, int &step)
{
    if (step_token.size() != 10 ||
        !std::all_of(step_token.begin(), step_token.end(),
                     [](unsigned char ch)
                     { return std::isdigit(ch) != 0; }))
    {
        return false;
    }
    try
    {
        step = std::stoi(step_token);
        return true;
    }
    catch (...)
    {
        return false;
    }
}

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

inline bool tryExtractRestartAnchorStep(const std::string &filename, int &step)
{
    return tryExtractFixedWidthRestartStep(filename, "Restart_time_", ".dat", step) ||
           tryExtractFixedWidthRestartStep(filename, "Restart_time_", ".xml", step) ||
           tryExtractFixedWidthRestartStep(filename, "Restart_", ".xml", step);
}

inline bool tryExtractRestartRelatedStep(const std::string &filename, int &step)
{
    if (tryExtractRestartAnchorStep(filename, step))
    {
        return true;
    }

    const std::string marker = "_rst_";
    const std::string suffix = ".xml";
    if (filename.size() <= marker.size() + 10 + suffix.size() ||
        filename.rfind(suffix) != filename.size() - suffix.size())
    {
        return false;
    }

    const size_t marker_pos = filename.rfind(marker);
    if (marker_pos == std::string::npos)
    {
        return false;
    }
    const size_t step_start = marker_pos + marker.size();
    const size_t step_length = filename.size() - suffix.size() - step_start;
    return step_length == 10 &&
           tryParseRestartStepToken(filename.substr(step_start, step_length), step);
}

inline int detectLatestRestartStep(const std::string &restart_dir = "restart")
{
    std::error_code ec;
    if (!std::filesystem::exists(restart_dir, ec) ||
        !std::filesystem::is_directory(restart_dir, ec))
    {
        std::cout << "Restart folder not found: " << restart_dir << std::endl;
        return 0;
    }

    int max_step = 0;
    std::string max_step_filename;
    for (const auto &entry : std::filesystem::directory_iterator(restart_dir, ec))
    {
        if (ec || !entry.is_regular_file(ec))
        {
            continue;
        }
        int step = 0;
        const std::string filename = entry.path().filename().string();
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
        std::cout << "No valid restart file found in " << restart_dir << std::endl;
    }
    return max_step;
}

inline void cleanupOldRestartCheckpoints(const std::string &restart_dir, int keep_last_n)
{
    if (keep_last_n < 0)
    {
        return;
    }

    std::error_code ec;
    if (!std::filesystem::exists(restart_dir, ec) ||
        !std::filesystem::is_directory(restart_dir, ec))
    {
        return;
    }

    std::vector<int> checkpoint_steps;
    std::unordered_set<int> checkpoint_step_set;
    std::unordered_map<int, std::vector<std::filesystem::path>> step_to_files;
    for (const auto &entry : std::filesystem::directory_iterator(restart_dir, ec))
    {
        if (ec || !entry.is_regular_file(ec))
        {
            continue;
        }
        int step = 0;
        const std::string filename = entry.path().filename().string();
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
        auto it = step_to_files.find(step);
        if (it == step_to_files.end())
        {
            continue;
        }
        for (const auto &path : it->second)
        {
            if (std::filesystem::remove(path, ec))
            {
                ++deleted_count;
            }
        }
    }
    if (deleted_count > 0)
    {
        std::cout << "[Cylinder3D][Restart] removed " << deleted_count
                  << " old checkpoint file(s); kept latest "
                  << keep_last_n << " checkpoint group(s)." << std::endl;
    }
}

inline std::filesystem::path resolveDefaultConfigPath(
    const std::filesystem::path &source_dir = std::filesystem::path(__FILE__).parent_path(),
    const std::filesystem::path &cwd = std::filesystem::current_path())
{
    std::error_code ec;
    const std::filesystem::path cwd_config = cwd / "config.ini";
    if (std::filesystem::exists(cwd_config, ec) && !ec)
    {
        return cwd_config;
    }
    ec.clear();
    const std::filesystem::path source_config = source_dir / "config.ini";
    if (std::filesystem::exists(source_config, ec) && !ec)
    {
        return source_config;
    }
    throw std::runtime_error("config.ini not found in current working directory or source directory.");
}

inline Vecd freestreamVelocity(const Cylinder3DConfig &cfg)
{
    return Vecd(cfg.u_f, 0.0, 0.0);
}

inline Real forceScale(const Cylinder3DConfig &cfg)
{
    return 0.5 * cfg.rho0_f * cfg.u_f * cfg.u_f * cfg.D * cfg.DW;
}

} // namespace cylinder_3d
} // namespace SPH

#endif // CYLINDER_3D_DATA_HPP
