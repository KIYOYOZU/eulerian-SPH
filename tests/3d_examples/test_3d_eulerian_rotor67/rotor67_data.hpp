/**
 * @file 	rotor67_data.hpp
 * @brief 	Configuration layer for the Rotor 67 Phase A smoke test.
 * @details Fully-compressible (Ideal-gas, HLLC-Riemann) Eulerian SPH in a
 *          rotating reference frame. Single-passage sector theta in
 *          [theta0, theta0 + delta_theta] about the Z axis, simplified
 *          constant-radius hub/casing walls (Phase A), circumferential
 *          rotating-periodic ghost boundary, annular-swirl inlet
 *          (relative velocity w = (Omega*y, -Omega*x, v_z)) and
 *          back-pressure outlet.
 *          INI parsing style mirrors eulerian_channel_data.hpp (trim_copy /
 *          to_lower_copy / strip_comment / section map).
 *
 *          Coordinate convention: particles are stored in the standard
 *          Cartesian frame (x, y, z); the rotation axis is Z; the sector
 *          spans [theta0, theta0 + delta_theta] in the x-y plane. The code
 *          vel_ field carries the rotating-frame relative velocity w.
 */
#ifndef ROTOR67_DATA_H
#define ROTOR67_DATA_H

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>

namespace SPH
{

// INI section/key -> value table, both lower-cased for case-insensitive lookup.
using RotorIniConfig = std::map<std::string, std::map<std::string, std::string>>;

inline constexpr int rotor67_literature_blade_count = 22;
inline constexpr Real rotor67_literature_design_rpm = 16043.0;
inline constexpr Real rotor67_literature_design_mass_flow = 33.25;
inline constexpr Real rotor67_literature_design_pressure_ratio = 1.63;
inline constexpr Real rotor67_literature_tip_clearance = 0.00101;

// Copy a string with leading/trailing whitespace removed.
inline std::string rotor_trim_copy(const std::string &value)
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
inline std::string rotor_to_lower_copy(const std::string &value)
{
    std::string out = value;
    std::transform(out.begin(), out.end(), out.begin(),
                   [](unsigned char ch)
                   { return static_cast<char>(std::tolower(ch)); });
    return out;
}

// Strip an inline '#' or ';' comment, keeping the original value before it.
inline std::string rotor_strip_comment(const std::string &value)
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
 * @brief Simulation configuration for the Rotor 67 Phase A smoke test.
 *        All fields are public with defaults so missing optional keys fall
 *        back gracefully.
 */
struct Rotor67Config
{
    // --- Case / claim profile ---
    std::string case_profile = "smoke"; ///< smoke or literature_clean.

    // --- Physical (Ideal gas) ---
    Real rho0_f = 0.0;               ///< Reference (inlet) density (kg/m^3).
    Real gamma = 1.4;                ///< Heat capacity ratio.
    Real inlet_axial_velocity = 0.0; ///< Inlet axial relative velocity v_x (m/s).
    Real inlet_total_pressure = 0.0; ///< Inlet total pressure (Pa).
    Real inlet_total_temperature = 0.0; ///< Inlet total temperature (K).
    Real gas_constant = 287.05;      ///< Specific gas constant R (J/(kg.K)).

    // --- Rotor geometry & operating point ---
    int n_blades = 22;               ///< Blade count (Rot.inf).
    Real rpm = 0.0;                  ///< Rotational speed (rev/min).
    Real omega = 0.0;                ///< Angular velocity Omega (rad/s) = 2*pi*rpm/60.
    Real hub_radius = 0.0;           ///< Hub radius (m), simplified constant.
    Real shroud_radius = 0.0;        ///< Shroud radius (m), simplified constant.
    Real blade_chord_max = 0.0;      ///< Max blade chord (m), Phase A placeholder.
    Real theta0 = 0.0;               ///< Sector start angle (rad).
    Real delta_theta = 0.0;          ///< Blade pitch angle = 2*pi/n_blades (rad).
    Real geometry_scale = 1.0;       ///< Geometry scale factor (record only).

    // --- Geometry (domain) ---
    Real z_in = 0.0;                 ///< Inlet axial coordinate (m).
    Real z_out = 0.0;                ///< Outlet axial coordinate (m).
    Real global_resolution = 0.0;    ///< Particle spacing dp (m).
    Real sponge_width_factor = 5.0;  ///< Sponge thickness (dp units).
    Real z_buffer_factor = 4.0;      ///< Hub/casing z buffer extension (dp units).

    // --- Inlet ---
    std::string inlet_profile = "annular_swirl"; ///< Only "annular_swirl" supported.
    std::string inlet_case_type = "clean"; ///< clean/pressure_distortion/swirl_distortion/combined.
    std::string swirl_axis = "Z";   ///< Rotation axis identifier.
    std::string inlet_flow_direction = "axial"; ///< Clean inlet baseline direction.
    std::string inlet_swirl_type = "none"; ///< none/co/counter absolute-frame pre-swirl.
    Real inlet_swirl_angle_deg = 0.0; ///< Absolute-frame pre-swirl angle.
    std::string inlet_total_pressure_profile = "uniform"; ///< Clean total-pressure profile.
    Real inlet_distortion_intensity = 0.0; ///< Total-pressure distortion intensity DI.

    // --- Outlet ---
    Real back_pressure = 0.0;        ///< Outlet back pressure (Pa).

    // --- Wall ---
    std::string casing_motion = "stationary_absolute"; ///< Casing wall motion contract.
    int wall_layers = 3; ///< Wall shell thickness in dp layers.

    // --- Simulation controls ---
    Real end_time = 0.0;             ///< End simulation time (s).
    int output_interval = 20;        ///< VTP output interval (acoustic steps).
    int screen_output_interval = 20; ///< Screen print interval.
    bool allow_high_mach = true;     ///< Bypass low-Mach guard (transonic).
    Real c_f = 0.0;                  ///< Unused for compressible path (record only).

    // --- Viscosity (Phase B) ---
    Real dynamic_viscosity = 0.0;    ///< Dynamic viscosity mu (Pa·s). 0 = inviscid.

    // --- Domain geometry (named-solid STL) ---
    std::string geometry_source = "domain_stl"; ///< Active geometry source label.
    std::string domain_stl_path; ///< ASCII multi-solid STL with named Rotor67 domain surfaces.
    Real particle_boundary_offset = 0.5; ///< First fluid/wall layer offset in dp units.

    // --- Meridional r(z) curves (Phase B+) ---
    /// Hub CFX .curve file (full-scale metres; loaded directly, no scaling).
    std::string meridional_hub_curve_path;
    /// Shroud CFX .curve file (full-scale metres; loaded directly, no scaling).
    std::string meridional_shroud_curve_path;

    // --- Boundary ---
    int boundary_n_layers = 3;       ///< Ghost-layer thickness (dp units).
    Real inlet_outlet_tolerance_factor = 3.0; ///< Axial open-boundary classifier band in dp units.
    Real periodic_tolerance_factor = 1.5;     ///< Circumferential classifier band in dp units.
    Real wall_near_tolerance_factor = 1.0;    ///< Hub/casing/blade near-wall classifier band in dp units.

    // --- Validation tolerances (trend indicators) ---
    Real mass_flux_imbalance_tol = 0.05;
    Real pressure_ratio_min = 1.2;
    Real pressure_ratio_max = 2.0;
};

// Read a UTF-8 INI file (BOM-safe) into a section/key/value lower-case index.
inline RotorIniConfig rotor_load_ini_file(const std::string &path)
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
            in.clear();
            in.seekg(0);
        }
    }

    RotorIniConfig ini;
    std::string current_section;
    std::string line;
    while (std::getline(in, line))
    {
        std::string cleaned = rotor_trim_copy(rotor_strip_comment(line));
        if (cleaned.empty())
        {
            continue;
        }
        if (cleaned.front() == '[' && cleaned.back() == ']')
        {
            current_section = rotor_to_lower_copy(rotor_trim_copy(cleaned.substr(1, cleaned.size() - 2)));
            continue;
        }
        const size_t eq_pos = cleaned.find('=');
        if (eq_pos == std::string::npos)
        {
            continue;
        }
        std::string key = rotor_to_lower_copy(rotor_trim_copy(cleaned.substr(0, eq_pos)));
        std::string value = rotor_trim_copy(cleaned.substr(eq_pos + 1));
        if (!current_section.empty() && !key.empty())
        {
            ini[current_section][key] = value;
        }
    }
    return ini;
}

inline Real rotor_parse_real(const std::string &value)
{
    try
    {
        return static_cast<Real>(std::stod(rotor_trim_copy(value)));
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error("Failed to parse Real value '" + value + "': " + e.what());
    }
}

inline int rotor_parse_int(const std::string &value)
{
    try
    {
        return std::stoi(rotor_trim_copy(value));
    }
    catch (const std::exception &e)
    {
        throw std::runtime_error("Failed to parse int value '" + value + "': " + e.what());
    }
}

inline bool rotor_parse_bool(const std::string &value)
{
    const std::string s = rotor_to_lower_copy(rotor_trim_copy(value));
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

inline std::string rotor_get_ini_value(const RotorIniConfig &ini,
                                       const std::string &section,
                                       const std::string &key)
{
    const std::string section_key = rotor_to_lower_copy(section);
    const std::string item_key = rotor_to_lower_copy(key);
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

inline Real rotor_get_real_or_default(const RotorIniConfig &ini, const std::string &section,
                                      const std::string &key, Real default_value)
{
    const std::string sec_lower = rotor_to_lower_copy(section);
    const std::string key_lower = rotor_to_lower_copy(key);
    auto sec_it = ini.find(sec_lower);
    if (sec_it == ini.end())
        return default_value;
    auto key_it = sec_it->second.find(key_lower);
    if (key_it == sec_it->second.end())
        return default_value;
    try
    {
        return rotor_parse_real(key_it->second);
    }
    catch (...)
    {
        return default_value;
    }
}

inline int rotor_get_int_or_default(const RotorIniConfig &ini, const std::string &section,
                                    const std::string &key, int default_value)
{
    const std::string sec_lower = rotor_to_lower_copy(section);
    const std::string key_lower = rotor_to_lower_copy(key);
    auto sec_it = ini.find(sec_lower);
    if (sec_it == ini.end())
        return default_value;
    auto key_it = sec_it->second.find(key_lower);
    if (key_it == sec_it->second.end())
        return default_value;
    try
    {
        return rotor_parse_int(key_it->second);
    }
    catch (...)
    {
        return default_value;
    }
}

inline bool rotor_get_bool_or_default(const RotorIniConfig &ini, const std::string &section,
                                      const std::string &key, bool default_value)
{
    const std::string sec_lower = rotor_to_lower_copy(section);
    const std::string key_lower = rotor_to_lower_copy(key);
    auto sec_it = ini.find(sec_lower);
    if (sec_it == ini.end())
        return default_value;
    auto key_it = sec_it->second.find(key_lower);
    if (key_it == sec_it->second.end())
        return default_value;
    try
    {
        return rotor_parse_bool(key_it->second);
    }
    catch (...)
    {
        return default_value;
    }
}

inline std::string rotor_get_string_or_default(const RotorIniConfig &ini, const std::string &section,
                                               const std::string &key, const std::string &default_value)
{
    const std::string sec_lower = rotor_to_lower_copy(section);
    const std::string key_lower = rotor_to_lower_copy(key);
    auto sec_it = ini.find(sec_lower);
    if (sec_it == ini.end())
        return default_value;
    auto key_it = sec_it->second.find(key_lower);
    if (key_it == sec_it->second.end())
        return default_value;
    return key_it->second;
}

inline bool rotor_ini_has_key(const RotorIniConfig &ini, const std::string &section,
                              const std::string &key)
{
    const std::string sec_lower = rotor_to_lower_copy(section);
    const std::string key_lower = rotor_to_lower_copy(key);
    auto sec_it = ini.find(sec_lower);
    return sec_it != ini.end() && sec_it->second.find(key_lower) != sec_it->second.end();
}

inline bool rotor67_is_allowed_distortion_intensity(Real value)
{
    return std::fabs(value) <= Real(1.0e-12) ||
           std::fabs(value - Real(0.05)) <= Real(1.0e-12) ||
           std::fabs(value - Real(0.10)) <= Real(1.0e-12);
}

inline std::filesystem::path rotor_resolve_case_path(const std::filesystem::path &case_dir,
                                                     const std::string &raw_path)
{
    std::filesystem::path resolved = raw_path;
    if (resolved.empty())
    {
        return resolved;
    }
    if (!resolved.is_absolute())
    {
        resolved = case_dir / resolved;
    }
    return resolved.lexically_normal();
}

inline void rotor_reject_active_old_geometry_keys(const RotorIniConfig &ini)
{
    static const std::set<std::string> forbidden_keys{
        "geometry_source", "step_mesh", "step_blade_wall_mesh", "step_mesh_report",
        "step_source", "step_solid_selector", "profile_curve", "profile_nsamples",
        "blade_stl_path", "stl_path", "blade_scale", "tip_clearance",
        "reload_normals", "blade_wall_band_half_width",
        "fluid_blade_clearance_half_width", "theta_blade"};
    auto sec_it = ini.find("blade");
    if (sec_it != ini.end())
    {
        for (const auto &entry : sec_it->second)
        {
            const std::string &key = entry.first;
            if (forbidden_keys.find(key) != forbidden_keys.end())
            {
                throw std::runtime_error(
                    "Forbidden active [blade] key '" + key +
                    "': Rotor67 initialization now uses [geometry] domain_stl only");
            }
        }
    }
    if (ini.find("step_domain") != ini.end())
    {
        throw std::runtime_error(
            "Forbidden active [step_domain] section: Rotor67 boundary geometry now derives from [geometry] domain_stl");
    }
    if (ini.find("reload") != ini.end())
    {
        throw std::runtime_error(
            "Forbidden active [reload] section: Rotor67 initialization no longer loads reload particles");
    }
    static const std::set<std::string> forbidden_reload_keys{
        "reload", "reload_xml", "reload_particles"};
    for (const auto &section : ini)
    {
        for (const auto &key : forbidden_reload_keys)
        {
            if (section.second.find(key) != section.second.end())
            {
                throw std::runtime_error(
                    "Forbidden active key '" + section.first + "." + key +
                    "': Rotor67 initialization no longer loads reload particles");
            }
        }
    }
}

inline std::filesystem::path rotor_resolve_default_config_path(
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
 * @brief Parse config.ini into Rotor67Config and run consistency validation.
 *
 * Validation:
 *   - n_blades >= 3, omega > 0 (derived from rpm if omega_rad_s absent).
 *   - hub/shroud radii positive, hub < shroud.
 *   - delta_theta derived as 2*pi/n_blades unless overridden.
 *   - back_pressure > inlet static pressure (isentropic from total P/T).
 *   - geometry positivity (z_in < z_out, global_resolution > 0).
 *   - inlet_profile == "annular_swirl".
 *   - allow_high_mach required true for the transonic rotor; if false and
 *     the relative tip Mach exceeds the low-Mach bound, throw.
 */
inline Rotor67Config load_rotor67_config(const std::string &config_path)
{
    RotorIniConfig ini = rotor_load_ini_file(config_path);

    Rotor67Config cfg;

    // --- [case] ---
    cfg.case_profile = rotor_to_lower_copy(
        rotor_get_string_or_default(ini, "case", "profile", "smoke"));

    // --- [physical] (required) ---
    cfg.rho0_f = rotor_parse_real(rotor_get_ini_value(ini, "physical", "rho0_f"));
    cfg.gamma = rotor_get_real_or_default(ini, "physical", "gamma", 1.4);
    cfg.inlet_axial_velocity =
        rotor_parse_real(rotor_get_ini_value(ini, "physical", "inlet_axial_velocity"));
    cfg.inlet_total_pressure =
        rotor_parse_real(rotor_get_ini_value(ini, "physical", "inlet_total_pressure"));
    cfg.inlet_total_temperature =
        rotor_parse_real(rotor_get_ini_value(ini, "physical", "inlet_total_temperature"));
    cfg.gas_constant =
        rotor_get_real_or_default(ini, "physical", "gas_constant", 287.05);

    // --- [rotor] (required) ---
    cfg.n_blades = rotor_parse_int(rotor_get_ini_value(ini, "rotor", "n_blades"));
    cfg.rpm = rotor_parse_real(rotor_get_ini_value(ini, "rotor", "rpm"));
    cfg.omega = rotor_get_real_or_default(ini, "rotor", "omega_rad_s", -1.0);
    cfg.hub_radius = rotor_parse_real(rotor_get_ini_value(ini, "rotor", "hub_radius"));
    cfg.shroud_radius = rotor_parse_real(rotor_get_ini_value(ini, "rotor", "shroud_radius"));
    cfg.blade_chord_max =
        rotor_get_real_or_default(ini, "rotor", "blade_chord_max", 0.0);
    cfg.theta0 = rotor_get_real_or_default(ini, "rotor", "theta0", 0.0);
    cfg.delta_theta = rotor_get_real_or_default(ini, "rotor", "delta_theta", -1.0);
    cfg.geometry_scale =
        rotor_get_real_or_default(ini, "rotor", "geometry_scale", 1.0);

    // --- [geometry] (required) ---
    cfg.z_in = rotor_parse_real(rotor_get_ini_value(ini, "geometry", "z_in"));
    cfg.z_out = rotor_parse_real(rotor_get_ini_value(ini, "geometry", "z_out"));
    cfg.global_resolution =
        rotor_parse_real(rotor_get_ini_value(ini, "geometry", "global_resolution"));
    {
        const std::filesystem::path case_dir =
            std::filesystem::path(config_path).parent_path();
        rotor_reject_active_old_geometry_keys(ini);
        cfg.domain_stl_path = rotor_resolve_case_path(
                                  case_dir,
                                  rotor_get_ini_value(ini, "geometry", "domain_stl"))
                                  .string();
    }
    cfg.particle_boundary_offset =
        rotor_get_real_or_default(ini, "geometry", "particle_boundary_offset", 0.5);
    cfg.sponge_width_factor =
        rotor_get_real_or_default(ini, "geometry", "sponge_width_factor", 5.0);
    cfg.z_buffer_factor =
        rotor_get_real_or_default(ini, "geometry", "z_buffer_factor", 4.0);

    // --- [inlet] ---
    cfg.inlet_profile = rotor_to_lower_copy(
        rotor_get_string_or_default(ini, "inlet", "profile", "annular_swirl"));
    cfg.inlet_case_type = rotor_to_lower_copy(
        rotor_get_string_or_default(ini, "inlet", "case_type", "clean"));
    cfg.swirl_axis = rotor_to_lower_copy(
        rotor_get_string_or_default(ini, "inlet", "swirl_axis", "Z"));
    cfg.inlet_flow_direction = rotor_to_lower_copy(
        rotor_get_string_or_default(ini, "inlet", "flow_direction", "axial"));
    cfg.inlet_swirl_type = rotor_to_lower_copy(
        rotor_get_string_or_default(ini, "inlet", "swirl_type", "none"));
    cfg.inlet_swirl_angle_deg =
        rotor_get_real_or_default(ini, "inlet", "swirl_angle_deg", 0.0);
    cfg.inlet_total_pressure_profile = rotor_to_lower_copy(
        rotor_get_string_or_default(ini, "inlet", "total_pressure_profile", "uniform"));
    cfg.inlet_distortion_intensity =
        rotor_get_real_or_default(ini, "inlet", "distortion_intensity", 0.0);

    // --- [outlet] ---
    cfg.back_pressure = rotor_parse_real(rotor_get_ini_value(ini, "outlet", "back_pressure"));

    // --- [wall] ---
    cfg.casing_motion = rotor_to_lower_copy(
        rotor_get_string_or_default(ini, "wall", "casing_motion", "stationary_absolute"));
    cfg.wall_layers =
        rotor_get_int_or_default(ini, "wall", "wall_layers", 3);

    // --- [simulation] ---
    cfg.end_time = rotor_parse_real(rotor_get_ini_value(ini, "simulation", "end_time"));
    cfg.output_interval =
        rotor_parse_int(rotor_get_ini_value(ini, "simulation", "output_interval"));
    cfg.screen_output_interval =
        rotor_parse_int(rotor_get_ini_value(ini, "simulation", "screen_output_interval"));
    cfg.allow_high_mach =
        rotor_get_bool_or_default(ini, "simulation", "allow_high_mach", true);
    cfg.c_f = rotor_get_real_or_default(ini, "simulation", "c_f", 0.0);

    // --- [boundary] ---
    cfg.boundary_n_layers =
        rotor_get_int_or_default(ini, "boundary", "boundary_n_layers", 3);
    cfg.inlet_outlet_tolerance_factor =
        rotor_get_real_or_default(ini, "boundary", "inlet_outlet_tolerance_factor",
                                  static_cast<Real>(cfg.boundary_n_layers));
    cfg.periodic_tolerance_factor =
        rotor_get_real_or_default(ini, "boundary", "periodic_tolerance_factor", 1.5);
    cfg.wall_near_tolerance_factor =
        rotor_get_real_or_default(ini, "boundary", "wall_near_tolerance_factor", 1.0);

    // --- [viscosity] (Phase B, optional) ---
    cfg.dynamic_viscosity =
        rotor_get_real_or_default(ini, "viscosity", "dynamic_viscosity", 0.0);

    // --- [meridional] (Phase B+, required for true hub/shroud r(z) profiles) ---
    {
        const std::filesystem::path case_dir =
            std::filesystem::path(config_path).parent_path();
        auto resolve_curve = [&](const std::string &raw)
        {
            std::filesystem::path resolved = raw;
            if (resolved.empty())
            {
                throw std::runtime_error(
                    "[meridional] curve path is empty; specify hub_curve and shroud_curve in config.ini");
            }
            if (!resolved.is_absolute())
            {
                resolved = case_dir / resolved;
            }
            return resolved.string();
        };
        const std::string hub_raw = rotor_get_string_or_default(
            ini, "meridional", "hub_curve", "docu/R67_CFX/Rot_Hub.curve");
        const std::string shd_raw = rotor_get_string_or_default(
            ini, "meridional", "shroud_curve", "docu/R67_CFX/Rot_Shd.curve");
        cfg.meridional_hub_curve_path = resolve_curve(hub_raw);
        cfg.meridional_shroud_curve_path = resolve_curve(shd_raw);
    }

    // --- [validation] ---
    cfg.mass_flux_imbalance_tol =
        rotor_get_real_or_default(ini, "validation", "mass_flux_imbalance_tol", 0.05);
    cfg.pressure_ratio_min =
        rotor_get_real_or_default(ini, "validation", "pressure_ratio_min", 1.2);
    cfg.pressure_ratio_max =
        rotor_get_real_or_default(ini, "validation", "pressure_ratio_max", 2.0);

    // --- Derived: omega and delta_theta ---
    if (cfg.omega <= 0.0)
    {
        cfg.omega = Real(2.0) * Pi * cfg.rpm / Real(60.0);
    }
    if (cfg.delta_theta <= 0.0)
    {
        cfg.delta_theta = Real(2.0) * Pi / static_cast<Real>(cfg.n_blades);
    }

    // --- Validation ---
    if (cfg.n_blades < 3)
    {
        throw std::runtime_error("Invalid rotor.n_blades: must be >= 3 (got " +
                                 std::to_string(cfg.n_blades) + ")");
    }
    if (cfg.case_profile != "smoke" && cfg.case_profile != "literature_clean")
    {
        throw std::runtime_error("Unsupported case.profile '" + cfg.case_profile +
                                 "': use 'smoke' or 'literature_clean'");
    }
    if (cfg.omega <= 0.0)
    {
        throw std::runtime_error("Invalid rotor omega: must be positive (got " +
                                 std::to_string(cfg.omega) + ")");
    }
    if (cfg.hub_radius <= 0.0 || cfg.shroud_radius <= 0.0 ||
        cfg.hub_radius >= cfg.shroud_radius)
    {
        throw std::runtime_error("Invalid rotor radii: require 0 < hub < shroud (hub=" +
                                 std::to_string(cfg.hub_radius) + ", shroud=" +
                                 std::to_string(cfg.shroud_radius) + ")");
    }
    if (cfg.global_resolution <= 0.0)
    {
        throw std::runtime_error("Invalid geometry.global_resolution: must be positive (got " +
                                 std::to_string(cfg.global_resolution) + ")");
    }
    if (cfg.z_in >= cfg.z_out)
    {
        throw std::runtime_error("Invalid geometry: require z_in < z_out (z_in=" +
                                 std::to_string(cfg.z_in) + ", z_out=" +
                                 std::to_string(cfg.z_out) + ")");
    }
    if (cfg.rho0_f <= 0.0 || cfg.gamma <= 1.0)
    {
        throw std::runtime_error("Invalid physical: rho0_f > 0 and gamma > 1 required");
    }
    if (cfg.inlet_axial_velocity < 0.0)
    {
        throw std::runtime_error("Invalid physical.inlet_axial_velocity: must be >= 0");
    }
    if (cfg.inlet_profile != "annular_swirl")
    {
        throw std::runtime_error("Unsupported inlet.profile '" + cfg.inlet_profile +
                                 "', only 'annular_swirl' is supported");
    }
    if (cfg.inlet_case_type != "clean" &&
        cfg.inlet_case_type != "pressure_distortion" &&
        cfg.inlet_case_type != "swirl_distortion" &&
        cfg.inlet_case_type != "combined")
    {
        throw std::runtime_error("Unsupported inlet.case_type '" + cfg.inlet_case_type +
                                 "': use clean, pressure_distortion, swirl_distortion, or combined");
    }
    if (cfg.inlet_flow_direction != "axial")
    {
        throw std::runtime_error("Unsupported inlet.flow_direction '" +
                                 cfg.inlet_flow_direction + "', only 'axial' is supported");
    }
    if (cfg.swirl_axis != "z")
    {
        throw std::runtime_error("Unsupported inlet.swirl_axis '" + cfg.swirl_axis +
                                 "': Rotor67 runtime uses the Z rotation axis");
    }
    if (cfg.inlet_total_pressure_profile != "uniform" &&
        cfg.inlet_total_pressure_profile != "tip_radial" &&
        cfg.inlet_total_pressure_profile != "hub_radial")
    {
        throw std::runtime_error("Unsupported inlet.total_pressure_profile '" +
                                 cfg.inlet_total_pressure_profile +
                                 "': use uniform, tip_radial, or hub_radial");
    }
    if (!rotor67_is_allowed_distortion_intensity(cfg.inlet_distortion_intensity))
    {
        throw std::runtime_error("Unsupported inlet.distortion_intensity: use 0.0, 0.05, or 0.10");
    }
    if (cfg.inlet_swirl_type != "none" &&
        cfg.inlet_swirl_type != "co" &&
        cfg.inlet_swirl_type != "counter")
    {
        throw std::runtime_error("Unsupported inlet.swirl_type '" + cfg.inlet_swirl_type +
                                 "': use none, co, or counter");
    }
    const bool pressure_profile_requested = cfg.inlet_total_pressure_profile != "uniform";
    const bool pressure_intensity_positive =
        cfg.inlet_distortion_intensity > Real(1.0e-12);
    const bool swirl_type_requested = cfg.inlet_swirl_type != "none";
    const bool swirl_angle_positive =
        std::fabs(cfg.inlet_swirl_angle_deg) > Real(1.0e-12);
    if (cfg.inlet_case_type == "clean" &&
        (pressure_profile_requested || pressure_intensity_positive ||
         swirl_type_requested || swirl_angle_positive))
    {
        throw std::runtime_error(
            "Invalid clean inlet contract: pressure distortion and swirl settings must be zero/none");
    }
    if (cfg.inlet_case_type == "pressure_distortion" &&
        (!pressure_profile_requested || !pressure_intensity_positive ||
         swirl_type_requested || swirl_angle_positive))
    {
        throw std::runtime_error(
            "Invalid pressure_distortion inlet contract: require tip_radial/hub_radial with DI > 0 and no swirl");
    }
    if (cfg.inlet_case_type == "swirl_distortion" &&
        (pressure_profile_requested || pressure_intensity_positive ||
         !swirl_type_requested || !swirl_angle_positive))
    {
        throw std::runtime_error(
            "Invalid swirl_distortion inlet contract: require co/counter swirl with nonzero angle and uniform total pressure");
    }
    if (cfg.inlet_case_type == "combined" &&
        (!pressure_profile_requested || !pressure_intensity_positive ||
         !swirl_type_requested || !swirl_angle_positive))
    {
        throw std::runtime_error(
            "Invalid combined inlet contract: require radial total-pressure distortion with DI > 0 and co/counter swirl with nonzero angle");
    }
    if (cfg.casing_motion != "stationary_absolute" && cfg.casing_motion != "rotating_frame_static")
    {
        throw std::runtime_error("Unsupported wall.casing_motion '" + cfg.casing_motion +
                                 "': use 'stationary_absolute' or 'rotating_frame_static'");
    }
    if (cfg.wall_layers < 3)
    {
        throw std::runtime_error("Invalid wall.wall_layers: Rotor67 baseline requires at least 3");
    }
    // Back pressure must exceed the inlet static pressure (isentropic from total).
    // p_static = p_total / (1 + (gamma-1)/2 * M^2)^(gamma/(gamma-1)), with the
    // inlet Mach based on the axial relative velocity and the sound speed at T_total.
    const Real cp_inlet = cfg.gamma * cfg.gas_constant / (cfg.gamma - Real(1.0));
    const Real T_static_inlet =
        cfg.inlet_total_temperature -
        cfg.inlet_axial_velocity * cfg.inlet_axial_velocity / (Real(2.0) * cp_inlet);
    if (T_static_inlet <= TinyReal)
    {
        throw std::runtime_error(
            "Invalid inlet total state: axial velocity exceeds total-temperature enthalpy");
    }
    const Real p_static_inlet = cfg.inlet_total_pressure *
        std::pow(T_static_inlet / cfg.inlet_total_temperature,
                 cfg.gamma / (cfg.gamma - Real(1.0)));
    if (cfg.back_pressure <= p_static_inlet)
    {
        throw std::runtime_error(
            "Invalid outlet.back_pressure: must exceed inlet static pressure " +
            std::to_string(p_static_inlet) + " (got " + std::to_string(cfg.back_pressure) + ")");
    }
    if (cfg.geometry_scale <= 0.0)
    {
        throw std::runtime_error("Invalid rotor.geometry_scale: must be positive");
    }
    if (cfg.geometry_source != "domain_stl")
    {
        throw std::runtime_error("Internal Rotor67 geometry_source mismatch: only 'domain_stl' is active");
    }
    if (cfg.domain_stl_path.empty())
    {
        throw std::runtime_error("Invalid [geometry] domain_stl: path must not be empty");
    }
    if (cfg.particle_boundary_offset <= 0.0)
    {
        throw std::runtime_error("Invalid geometry.particle_boundary_offset: must be positive");
    }
    if (cfg.inlet_outlet_tolerance_factor <= 0.0 ||
        cfg.periodic_tolerance_factor <= 0.0 ||
        cfg.wall_near_tolerance_factor <= 0.0)
    {
        throw std::runtime_error("Invalid [boundary] classifier tolerances: all must be positive");
    }

    return cfg;
}

/**
 * @brief Inlet static pressure (isentropic from total P/T and axial Mach).
 */
inline Real rotorInletStaticPressure(const Rotor67Config &cfg)
{
    const Real cp = cfg.gamma * cfg.gas_constant / (cfg.gamma - Real(1.0));
    const Real T_static =
        cfg.inlet_total_temperature -
        cfg.inlet_axial_velocity * cfg.inlet_axial_velocity / (Real(2.0) * cp);
    if (T_static <= TinyReal)
    {
        throw std::runtime_error("Invalid inlet total state: axial velocity exceeds total-temperature enthalpy");
    }
    return cfg.inlet_total_pressure *
           std::pow(T_static / cfg.inlet_total_temperature,
                    cfg.gamma / (cfg.gamma - Real(1.0)));
}

/**
 * @brief Inlet static density (isentropic from total P/T and axial Mach).
 */
inline Real rotorInletStaticDensity(const Rotor67Config &cfg)
{
    const Real p_static = rotorInletStaticPressure(cfg);
    const Real cp = cfg.gamma * cfg.gas_constant / (cfg.gamma - Real(1.0));
    const Real T_static =
        cfg.inlet_total_temperature -
        cfg.inlet_axial_velocity * cfg.inlet_axial_velocity / (Real(2.0) * cp);
    if (T_static <= TinyReal)
    {
        throw std::runtime_error("Invalid inlet total state: axial velocity exceeds total-temperature enthalpy");
    }
    return p_static / (cfg.gas_constant * T_static);
}

/**
 * @brief Inlet static temperature (isentropic from total T and axial Mach).
 */
inline Real rotorInletStaticTemperature(const Rotor67Config &cfg)
{
    const Real cp = cfg.gamma * cfg.gas_constant / (cfg.gamma - Real(1.0));
    const Real T_static =
        cfg.inlet_total_temperature -
        cfg.inlet_axial_velocity * cfg.inlet_axial_velocity / (Real(2.0) * cp);
    if (T_static <= TinyReal)
    {
        throw std::runtime_error("Invalid inlet total state: axial velocity exceeds total-temperature enthalpy");
    }
    return T_static;
}

/**
 * @brief Annular-swirl relative velocity at a Cartesian position (no pre-swirl).
 *        Rotation axis is Z with angular velocity Omega along +Z. The absolute
 *        circumferential velocity of the rotating frame is Omega x r =
 *        (-Omega*y, Omega*x, 0). With no inlet pre-swirl the absolute tangential
 *        velocity is zero, so the relative (rotating-frame) tangential velocity is
 *        the negative of the frame velocity: (Omega*y, -Omega*x, 0). The axial
 *        relative velocity equals the absolute axial velocity v_x (along z).
 *        Hence w(pos) = (Omega*y, -Omega*x, v_x) (m/s).
 * @param pos  Particle position in Cartesian (x, y, z).
 * @param cfg  Configuration carrying v_x and Omega.
 * @return Relative velocity vector w (m/s).
 */
inline Vecd rotorAnnularSwirlVelocity(const Vecd &pos, const Rotor67Config &cfg)
{
    const Real x = pos[0];
    const Real y = pos[1];
    // w = (Omega*y, -Omega*x, v_x); the axial component carries v_x along z.
    return Vecd(cfg.omega * y, -cfg.omega * x, cfg.inlet_axial_velocity);
}

} // namespace SPH

#endif // ROTOR67_DATA_H
