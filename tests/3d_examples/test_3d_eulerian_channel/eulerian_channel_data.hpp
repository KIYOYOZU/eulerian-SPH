/**
 * @file 	eulerian_channel_data.hpp
 * @brief 	Configuration layer for the 3D Eulerian weakly-compressible laminar channel flow.
 * @details Weakly-compressible Eulerian SPH, 3D laminar channel. Flow direction x,
 *          wall-normal y, spanwise z (periodic). Inlet parabolic profile
 *          u_x(y) = U_max * 4 * eta * (1 - eta), eta = y / DH, U_max = 1.5 * U_bulk.
 *          INI parsing style mirrors wmles_channel_data.hpp (trim_copy / to_lower_copy /
 *          strip_comment / section map) but is simplified for this smoke-test case.
 */
#ifndef EULERIAN_CHANNEL_DATA_H
#define EULERIAN_CHANNEL_DATA_H

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>

namespace SPH
{

// INI section/key -> value table, both lower-cased for case-insensitive lookup.
using IniConfig = std::map<std::string, std::map<std::string, std::string>>;

// Copy a string with leading/trailing whitespace removed.
inline std::string trim_copy(const std::string &value)
{
    const auto begin = value.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos)
    {
        return "";
    }
    const auto end = value.find_last_not_of(" \t\r\n");
    return value.substr(begin, end - begin + 1);
}

// Copy a string and convert to lower case (section/key/enum parsing).
inline std::string to_lower_copy(const std::string &value)
{
    std::string out = value;
    std::transform(out.begin(), out.end(), out.begin(),
                   [](unsigned char ch)
                   { return static_cast<char>(std::tolower(ch)); });
    return out;
}

// Strip an inline '#' or ';' comment, keeping the original value before it.
inline std::string strip_comment(const std::string &value)
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
        cut = (cut == std::string::npos) ? semi_pos : std::min(cut, semi_pos);
    }
    if (cut == std::string::npos)
    {
        return value;
    }
    return value.substr(0, cut);
}

/**
 * @brief Simulation configuration for the 3D Eulerian channel smoke-test case.
 *        All fields are public with defaults so missing optional keys fall back gracefully.
 */
struct SimulationConfig
{
    // --- Geometry (m) ---
    Real DL = 0.0;               ///< Streamwise length.
    Real DH = 0.0;               ///< Wall-normal height.
    Real DW = 0.0;               ///< Spanwise width.
    Real global_resolution = 0.0; ///< Particle spacing dp.
    Real sponge_width_factor = 5.0; ///< Sponge layer thickness in dp units (x-direction inlet/outlet), aligned with 2D LG.
    Real z_buffer_factor = 4.0; ///< z-direction wall body buffer extension in dp units (周期方向，wall z 外延消除端面退化).

    // --- Physical parameters ---
    Real rho0_f = 0.0; ///< Reference density.
    Real u_bulk = 0.0; ///< Bulk (volume-averaged) velocity.
    Real c_f = 0.0;    ///< Artificial sound speed.
    Real nu = 0.0;     ///< Kinematic viscosity.

    // --- Inlet ---
    std::string inlet_profile = "parabolic"; ///< Only "parabolic" is supported.
    Real inlet_relaxation_rate = 1.0;        ///< Reserved for Task 2/4 blending.

    // --- Simulation controls ---
    Real end_time = 0.0;          ///< End simulation time (s).
    Real acoustic_cfl = 0.25;     ///< Acoustic CFL.
    int output_interval = 20;     ///< Reload-step interval for VTP output.
    int screen_output_interval = 10; ///< Screen print interval (iterations).
    bool allow_high_mach = false; ///< If true, bypass the low-Mach guard.

    // --- Boundary ---
    std::string outlet_pressure_mode = "zero_gradient"; ///< Only "zero_gradient" is supported.
    int boundary_n_layers = 3; ///< Ghost-layer thickness (in dp units) used for face classification.

    // --- Validation tolerances ---
    Real inlet_profile_l2_tol = 0.05; ///< L2 tolerance for inlet profile check.
    Real wall_slip_tol = 0.01;        ///< Tolerance for wall no-slip check.
};

// Read a UTF-8 INI file (BOM-safe) into a section/key/value lower-case index.
inline IniConfig load_ini_file(const std::string &path)
{
    std::ifstream in(path, std::ios::in);
    if (!in.is_open())
    {
        throw std::runtime_error("Cannot open config file: " + path);
    }

    // Skip a leading UTF-8 BOM (EF BB BF) if present.
    int first = in.peek();
    if (first == 0xEF)
    {
        char bom[3];
        in.read(bom, 3);
        if (static_cast<unsigned char>(bom[0]) != 0xEF ||
            static_cast<unsigned char>(bom[1]) != 0xBB ||
            static_cast<unsigned char>(bom[2]) != 0xBF)
        {
            // Not a BOM: rewind so the bytes are re-read as content.
            in.clear();
            in.seekg(0);
        }
    }

    IniConfig ini;
    std::string current_section;
    std::string line;
    while (std::getline(in, line))
    {
        std::string cleaned = trim_copy(strip_comment(line));
        if (cleaned.empty())
        {
            continue;
        }
        if (cleaned.front() == '[' && cleaned.back() == ']')
        {
            current_section = to_lower_copy(trim_copy(cleaned.substr(1, cleaned.size() - 2)));
            continue;
        }
        const size_t eq_pos = cleaned.find('=');
        if (eq_pos == std::string::npos)
        {
            continue;
        }
        std::string key = to_lower_copy(trim_copy(cleaned.substr(0, eq_pos)));
        std::string value = trim_copy(cleaned.substr(eq_pos + 1));
        if (!current_section.empty() && !key.empty())
        {
            ini[current_section][key] = value;
        }
    }
    return ini;
}

// Parse a Real config value (plain decimal, no pi-expression support needed here).
inline Real parse_real_value(const std::string &value)
{
    try
    {
        return static_cast<Real>(std::stod(trim_copy(value)));
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error("Failed to parse Real value '" + value + "': " + e.what());
    }
}

// Parse an integer config value.
inline int parse_int_value(const std::string &value)
{
    try
    {
        return std::stoi(trim_copy(value));
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error("Failed to parse int value '" + value + "': " + e.what());
    }
}

// Parse a boolean config value: true/false, 1/0, yes/no, on/off.
inline bool parse_bool_value(const std::string &value)
{
    const std::string s = to_lower_copy(trim_copy(value));
    if (s == "true" || s == "1" || s == "yes" || s == "on")
    {
        return true;
    }
    if (s == "false" || s == "0" || s == "no" || s == "off")
    {
        return false;
    }
    throw std::runtime_error("Cannot parse bool value '" + value + "'");
}

// Read a required INI key; throws with section/key context if missing.
inline std::string get_ini_value(const IniConfig &ini,
                                 const std::string &section,
                                 const std::string &key)
{
    const std::string section_key = to_lower_copy(section);
    const std::string item_key = to_lower_copy(key);
    auto sec_it = ini.find(section_key);
    if (sec_it == ini.end())
    {
        throw std::runtime_error("Config missing section: [" + section + "]");
    }
    auto item_it = sec_it->second.find(item_key);
    if (item_it == sec_it->second.end())
    {
        throw std::runtime_error("Config missing key '" + key + "' in section [" + section + "]");
    }
    return item_it->second;
}

// Read an optional Real key with a default fallback.
inline Real get_ini_real_or_default(const IniConfig &ini, const std::string &section,
                                    const std::string &key, Real default_value)
{
    const std::string sec_lower = to_lower_copy(section);
    const std::string key_lower = to_lower_copy(key);
    auto sec_it = ini.find(sec_lower);
    if (sec_it == ini.end())
        return default_value;
    auto key_it = sec_it->second.find(key_lower);
    if (key_it == sec_it->second.end())
        return default_value;
    try
    {
        return parse_real_value(key_it->second);
    }
    catch (...)
    {
        return default_value;
    }
}

// Read an optional int key with a default fallback.
inline int get_ini_int_or_default(const IniConfig &ini, const std::string &section,
                                  const std::string &key, int default_value)
{
    const std::string sec_lower = to_lower_copy(section);
    const std::string key_lower = to_lower_copy(key);
    auto sec_it = ini.find(sec_lower);
    if (sec_it == ini.end())
        return default_value;
    auto key_it = sec_it->second.find(key_lower);
    if (key_it == sec_it->second.end())
        return default_value;
    try
    {
        return parse_int_value(key_it->second);
    }
    catch (...)
    {
        return default_value;
    }
}

// Read an optional bool key with a default fallback.
inline bool get_ini_bool_or_default(const IniConfig &ini, const std::string &section,
                                    const std::string &key, bool default_value)
{
    const std::string sec_lower = to_lower_copy(section);
    const std::string key_lower = to_lower_copy(key);
    auto sec_it = ini.find(sec_lower);
    if (sec_it == ini.end())
        return default_value;
    auto key_it = sec_it->second.find(key_lower);
    if (key_it == sec_it->second.end())
        return default_value;
    try
    {
        return parse_bool_value(key_it->second);
    }
    catch (...)
    {
        return default_value;
    }
}

// Read an optional string key with a default fallback.
inline std::string get_ini_string_or_default(const IniConfig &ini, const std::string &section,
                                             const std::string &key, const std::string &default_value)
{
    const std::string sec_lower = to_lower_copy(section);
    const std::string key_lower = to_lower_copy(key);
    auto sec_it = ini.find(sec_lower);
    if (sec_it == ini.end())
        return default_value;
    auto key_it = sec_it->second.find(key_lower);
    if (key_it == sec_it->second.end())
        return default_value;
    return key_it->second;
}

/**
 * @brief Resolve config.ini path: current working directory first, source directory fallback.
 */
inline std::filesystem::path resolve_default_config_path(
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
    throw std::runtime_error("config.ini not found in current working directory or source directory");
}

/**
 * @brief Parse config.ini into SimulationConfig and run consistency validation.
 *
 * Validation:
 *   - Geometry positivity (DL, DH, DW, global_resolution > 0).
 *   - Acoustic CFL positivity.
 *   - inlet_profile must be "parabolic".
 *   - outlet_pressure_mode must be "zero_gradient".
 *   - Low-Mach guard: c_f >= 10 * max(u_bulk, 1.5 * u_bulk) unless allow_high_mach.
 */
inline SimulationConfig load_config(const std::string &config_path)
{
    IniConfig ini = load_ini_file(config_path);

    SimulationConfig cfg;

    // --- [geometry] (required) ---
    cfg.DL = parse_real_value(get_ini_value(ini, "geometry", "dl"));
    cfg.DH = parse_real_value(get_ini_value(ini, "geometry", "dh"));
    cfg.DW = parse_real_value(get_ini_value(ini, "geometry", "dw"));
    cfg.global_resolution = parse_real_value(get_ini_value(ini, "geometry", "global_resolution"));
    cfg.sponge_width_factor =
        get_ini_real_or_default(ini, "geometry", "sponge_width_factor", 5.0);
    cfg.z_buffer_factor =
        get_ini_real_or_default(ini, "geometry", "z_buffer_factor", 4.0);

    // --- [physical] (required) ---
    cfg.rho0_f = parse_real_value(get_ini_value(ini, "physical", "rho0_f"));
    cfg.u_bulk = parse_real_value(get_ini_value(ini, "physical", "u_bulk"));
    cfg.c_f = parse_real_value(get_ini_value(ini, "physical", "c_f"));
    cfg.nu = parse_real_value(get_ini_value(ini, "physical", "nu"));

    // --- [inlet] (optional with defaults) ---
    cfg.inlet_profile = to_lower_copy(
        get_ini_string_or_default(ini, "inlet", "profile", "parabolic"));
    cfg.inlet_relaxation_rate =
        get_ini_real_or_default(ini, "inlet", "relaxation_rate", 1.0);

    // --- [simulation] (required + optional) ---
    cfg.end_time = parse_real_value(get_ini_value(ini, "simulation", "end_time"));
    cfg.acoustic_cfl = parse_real_value(get_ini_value(ini, "simulation", "acoustic_cfl"));
    cfg.output_interval = parse_int_value(get_ini_value(ini, "simulation", "output_interval"));
    cfg.screen_output_interval = parse_int_value(
        get_ini_value(ini, "simulation", "screen_output_interval"));
    cfg.allow_high_mach =
        get_ini_bool_or_default(ini, "simulation", "allow_high_mach", false);

    // --- [boundary] (optional with defaults) ---
    cfg.outlet_pressure_mode = to_lower_copy(
        get_ini_string_or_default(ini, "boundary", "outlet_pressure_mode", "zero_gradient"));
    cfg.boundary_n_layers =
        get_ini_int_or_default(ini, "boundary", "boundary_n_layers", 3);

    // --- [validation] (optional with defaults) ---
    cfg.inlet_profile_l2_tol =
        get_ini_real_or_default(ini, "validation", "inlet_profile_l2_tol", 0.05);
    cfg.wall_slip_tol =
        get_ini_real_or_default(ini, "validation", "wall_slip_tol", 0.01);

    // --- Validation: geometry positivity ---
    if (cfg.DL <= 0.0 || cfg.DH <= 0.0 || cfg.DW <= 0.0)
    {
        throw std::runtime_error(
            "Invalid geometry: DL, DH, DW must be positive (got DL=" +
            std::to_string(cfg.DL) + ", DH=" + std::to_string(cfg.DH) +
            ", DW=" + std::to_string(cfg.DW) + ")");
    }
    if (cfg.global_resolution <= 0.0)
    {
        throw std::runtime_error(
            "Invalid geometry: global_resolution must be positive (got " +
            std::to_string(cfg.global_resolution) + ")");
    }
    if (cfg.sponge_width_factor <= 0.0)
    {
        throw std::runtime_error(
            "Invalid geometry: sponge_width_factor must be positive (got " +
            std::to_string(cfg.sponge_width_factor) + ")");
    }
    if (cfg.z_buffer_factor <= 0.0)
    {
        throw std::runtime_error(
            "Invalid geometry: z_buffer_factor must be positive (got " +
            std::to_string(cfg.z_buffer_factor) + ")");
    }
    if (cfg.rho0_f <= 0.0 || cfg.u_bulk <= 0.0 || cfg.c_f <= 0.0 || cfg.nu <= 0.0)
    {
        throw std::runtime_error(
            "Invalid physical parameters: rho0_f, u_bulk, c_f, nu must be positive");
    }

    // --- Validation: acoustic CFL positivity ---
    if (cfg.acoustic_cfl <= 0.0)
    {
        throw std::runtime_error(
            "Invalid simulation.acoustic_cfl: must be positive (got " +
            std::to_string(cfg.acoustic_cfl) + ")");
    }

    // --- Validation: inlet profile enum ---
    if (cfg.inlet_profile != "parabolic")
    {
        throw std::runtime_error(
            "Unsupported inlet.profile '" + cfg.inlet_profile +
            "', only 'parabolic' is supported");
    }

    // --- Validation: outlet pressure mode enum ---
    if (cfg.outlet_pressure_mode != "zero_gradient")
    {
        throw std::runtime_error(
            "Unsupported boundary.outlet_pressure_mode '" + cfg.outlet_pressure_mode +
            "', only 'zero_gradient' is supported");
    }

    // --- Validation: low-Mach guard ---
    // Weakly-compressible SPH requires c_f >= 10 * max velocity scale to keep
    // density fluctuations below ~1%. The peak parabolic velocity is 1.5 * u_bulk.
    const Real u_max_scale = std::max(cfg.u_bulk, 1.5 * cfg.u_bulk);
    if (!cfg.allow_high_mach && cfg.c_f < 10.0 * u_max_scale)
    {
        throw std::runtime_error(
            "Low-Mach guard violated: c_f=" + std::to_string(cfg.c_f) +
            " < 10 * max(u_bulk, 1.5*u_bulk)=" + std::to_string(10.0 * u_max_scale) +
            ". Raise c_f or set simulation.allow_high_mach=true.");
    }

    return cfg;
}

/**
 * @brief Inlet peak velocity for the parabolic profile: U_max = 1.5 * u_bulk.
 *        With u_x(y) = U_max * 4 * eta * (1 - eta), the bulk average is 2/3 * U_max,
 *        so 1.5 * u_bulk gives u_bulk as the bulk-averaged velocity.
 */
inline Real inletUMax(const SimulationConfig &cfg)
{
    return 1.5 * cfg.u_bulk;
}

/**
 * @brief Parabolic inlet/streamwise velocity at wall-normal coordinate y.
 *        eta = y / DH, u_x = U_max * 4 * eta * (1 - eta).
 *        Vanishes at y=0 and y=DH (no-slip walls), peaks at y=DH/2.
 */
inline Real inletProfileVelocity(Real y, Real DH, Real U_max)
{
    if (DH <= 0.0)
    {
        return 0.0;
    }
    const Real eta = y / DH;
    return U_max * 4.0 * eta * (1.0 - eta);
}

/**
 * @brief Analytical bulk-averaged velocity of the parabolic profile.
 *        u_bulk = (1/DH) * integral_0^DH u_x(y) dy = 2/3 * U_max.
 *        Used to verify the U_max = 1.5 * u_bulk relationship.
 */
inline Real inletProfileBulk(Real DH, Real U_max)
{
    (void)DH; // The bulk average is independent of DH for this profile.
    return (2.0 / 3.0) * U_max;
}

} // namespace SPH

#endif // EULERIAN_CHANNEL_DATA_H
