#ifndef CYLINDER_3D_COMPRESSIBLE_DATA_HPP
#define CYLINDER_3D_COMPRESSIBLE_DATA_HPP

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
namespace cylinder_3d_compressible
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

inline bool hasKey(const IniConfig &ini, const std::string &section, const std::string &key)
{
    auto section_it = ini.find(toLowerCopy(section));
    if (section_it == ini.end())
    {
        return false;
    }
    return section_it->second.find(toLowerCopy(key)) != section_it->second.end();
}

/** X faces either use the Mach-aware far field or become periodic. */
enum class XBoundaryMode
{
    FarField,
    Periodic
};

inline const char *xBoundaryModeName(XBoundaryMode mode)
{
    return mode == XBoundaryMode::Periodic ? "periodic" : "farfield";
}

inline XBoundaryMode parseXBoundaryMode(const std::string &value)
{
    const std::string normalized = toLowerCopy(trimCopy(value));
    if (normalized == "farfield")
    {
        return XBoundaryMode::FarField;
    }
    if (normalized == "periodic")
    {
        return XBoundaryMode::Periodic;
    }
    throw std::runtime_error("Invalid x_boundary: expected 'farfield' or 'periodic', got '" + value + "'.");
}

/** Y faces either use the Mach-aware far field or become static walls. */
enum class YBoundaryMode
{
    FarField,
    Wall
};

inline const char *yBoundaryModeName(YBoundaryMode mode)
{
    return mode == YBoundaryMode::Wall ? "wall" : "farfield";
}

inline YBoundaryMode parseYBoundaryMode(const std::string &value)
{
    const std::string normalized = toLowerCopy(trimCopy(value));
    if (normalized == "farfield")
    {
        return YBoundaryMode::FarField;
    }
    if (normalized == "wall")
    {
        return YBoundaryMode::Wall;
    }
    throw std::runtime_error("Invalid y_boundary: expected 'farfield' or 'wall', got '" + value + "'.");
}

/**
 * Single source of truth for the configurable-Mach ideal-gas freestream.
 *
 * Inputs are gamma, rho_inf, c_inf, mach_inf and re; u_inf, p_inf, the total
 * energy density and mu are derived, never configured independently.
 */
struct CompressibleCylinderConfig
{
    Real DL = 0.0;
    Real DH = 0.0;
    Real DW = 0.0;
    Real global_resolution = 0.0;
    Real cylinder_center_x = 0.0;
    Real cylinder_center_y = 0.0;
    Real cylinder_radius = 0.0;
    int max_local_particles = 1000000;

    Real gamma = 0.0;
    Real rho_inf = 0.0;
    Real c_inf = 0.0;
    Real mach_inf = 0.0;
    Real re = 0.0;

    Real end_time = 0.0;
    int output_interval = 1;
    int screen_output_interval = 10;
    Real acoustic_cfl = 0.1;
    int boundary_n_layers = 3;
    Real ghost_reserve_factor = 0.5;
    bool enable_restart = false;
    int restart_step = 0;
    int restart_output_factor = 10;
    int restart_keep_last_n = 3;
    XBoundaryMode x_boundary_mode = XBoundaryMode::FarField;
    YBoundaryMode y_boundary_mode = YBoundaryMode::FarField;

    int smoke_min_steps = 50;
    /**
     * Relative amplitude of the streamwise velocity perturbation in the initial
     * condition.
     *
     * A non-zero value breaks the y-symmetry so the Re=100 wake can go unstable:
     * at 0 the dp=0.04 / t=49.3 run stayed on the symmetric steady branch with
     * Cl ~ 5e-7 and never shed.
     *
     * One zero-residual measurement is NOT exercised when it is non-zero: the
     * wall-flux absolute floor relies on the mirrored ghost state cancelling on a
     * uniform field, so with a scatter the gate falls through to its relative
     * branch. The t=0 pressure-closure baseline is unaffected, because this case's
     * DissipativePJump is identically zero (see the baseline comment in the main
     * translation unit).
     */
    Real velocity_noise_amplitude = 0.0;
    Real mass_drift_tolerance = 1.0e-4;
    Real energy_drift_tolerance = 1.0e-3;
    Real wall_flux_relative_tolerance = 1.0e-2;
    /**
     * Absolute floor for the wall net flux, as a fraction of the reference mass
     * flow rho_inf*u_inf*D*DW. Below it the wall is provably tight in absolute
     * terms and the relative criterion (which degenerates when net and gross are
     * both at machine precision) is bypassed.
     */
    Real wall_flux_absolute_floor_factor = 1.0e-10;
    /**
     * Diagnostic spatial A/B: when true, force first-order reconstruction in
     * the kernel-support band of the physical y walls. This is intentionally
     * opt-in: it is an isolated test of the wall-adjacent MUSCL path, not a
     * pressure or conservative-state correction.
     */
    bool y_wall_first_order_reconstruction = false;
    Real y_wall_first_order_band_dp = 2.6;
    /**
     * Diagnostic spatial A/B for the mixed x-outlet/y-wall corners only.
     * This keeps the inlet corners and the y-wall mid-span at full MUSCL order.
     */
    bool outlet_y_wall_corner_first_order_reconstruction = false;
    Real outlet_corner_first_order_band_dp = 2.6;
    /**
     * Diagnostic spatial A/B selected by the actual wall-contact relation.
     * Every fluid owner with at least one neighbour from any configured wall
     * body uses cell-centred first-order reconstruction.
     */
    bool wall_contact_first_order_reconstruction = false;
    int wall_flux_probe_interval = 50;
    /**
     * Stability A/B: route the MUSCL-HLLC bridge through the
     * dissipation-limited HLLC solver (Roe-averaged widened wave speeds)
     * instead of the plain Davis-speed HLLC. The plain path upwinds
     * supersonic normal flow at mirrored wall ghosts (s_l = u_n - c > 0),
     * so the rear cylinder face leaks mass/energy through the wall from
     * t = 0 at Ma = 2 -- invisible at Ma = 0.3 where |u_n| < c always
     * holds. See findings.md "Mach 2 vs 2D" for the evidence.
     */
    bool muscl_hllc_dissipation_limiter = false;
    Real muscl_hllc_limiter_parameter = 5.0;
    /**
     * Positivity floor (Zhang--Shu style): after each conservative update,
     * cells whose EOS-recovered pressure or density falls below
     * factor x freestream are reset to the floor with a consistent
     * TotalEnergy rewrite (Cylinder3DPositivityFloor). The A3 limited HLLC
     * leaves isolated bow-shock MUSCL undershoot cells at p ~ -0.1*p_inf;
     * the floor turns those into a bounded, visible correction so the
     * strict positivity gates stay green. 0 disables the correction.
     */
    bool positivity_floor = false;
    Real positivity_floor_factor = 1.0e-2;

    // Derived quantities -- written only by deriveAndValidate().
    Real u_inf = 0.0;
    Real p_inf = 0.0;
    Real E_inf_per_volume = 0.0;
    Real mu_inf = 0.0;
    Real dp = 0.0;
    Real D = 0.0;
};

/** Primitive plus conservative freestream state derived from the config. */
struct CompressibleFreestreamState
{
    Real gamma = 0.0;
    Real rho = 0.0;
    Real p = 0.0;
    Real sound_speed = 0.0;
    Real mach = 0.0;
    Vecd vel = Vecd::Zero();
    /** Total energy per unit volume: p/(gamma-1) + 0.5*rho*|u|^2. */
    Real E_per_volume = 0.0;
};

inline void deriveAndValidate(CompressibleCylinderConfig &cfg)
{
    // ---- geometry ----
    cfg.dp = cfg.global_resolution;
    cfg.D = 2.0 * cfg.cylinder_radius;

    if (cfg.DL <= 0.0 || cfg.DH <= 0.0 || cfg.DW <= 0.0 || cfg.dp <= 0.0)
    {
        throw std::runtime_error("Invalid geometry: dl, dh, dw and global_resolution must be positive.");
    }
    if (cfg.cylinder_radius <= 0.0)
    {
        throw std::runtime_error("Invalid geometry: cylinder_radius must be positive.");
    }
    if (cfg.y_wall_first_order_reconstruction && cfg.y_boundary_mode != YBoundaryMode::Wall)
    {
        throw std::runtime_error("y_wall_first_order_reconstruction requires y_boundary=wall.");
    }
    if (cfg.y_wall_first_order_reconstruction && cfg.y_wall_first_order_band_dp <= 0.0)
    {
        throw std::runtime_error("y_wall_first_order_band_dp must be positive when the y-wall A/B is enabled.");
    }
    const int first_order_reconstruction_modes =
        static_cast<int>(cfg.y_wall_first_order_reconstruction) +
        static_cast<int>(cfg.outlet_y_wall_corner_first_order_reconstruction) +
        static_cast<int>(cfg.wall_contact_first_order_reconstruction);
    if (first_order_reconstruction_modes > 1)
    {
        throw std::runtime_error(
            "The full y-wall, outlet-corner and wall-contact first-order A/B modes are mutually exclusive.");
    }
    if (cfg.outlet_y_wall_corner_first_order_reconstruction &&
        cfg.y_boundary_mode != YBoundaryMode::Wall)
    {
        throw std::runtime_error(
            "outlet_y_wall_corner_first_order_reconstruction requires y_boundary=wall.");
    }
    if (cfg.outlet_y_wall_corner_first_order_reconstruction &&
        cfg.x_boundary_mode != XBoundaryMode::FarField)
    {
        throw std::runtime_error(
            "outlet_y_wall_corner_first_order_reconstruction requires x_boundary=farfield.");
    }
    if (cfg.outlet_y_wall_corner_first_order_reconstruction &&
        cfg.outlet_corner_first_order_band_dp <= 0.0)
    {
        throw std::runtime_error(
            "outlet_corner_first_order_band_dp must be positive when the outlet-corner A/B is enabled.");
    }
    if (cfg.outlet_y_wall_corner_first_order_reconstruction &&
        cfg.y_wall_first_order_band_dp <= 0.0)
    {
        throw std::runtime_error(
            "y_wall_first_order_band_dp must be positive when the outlet-corner A/B is enabled.");
    }
    if (cfg.outlet_y_wall_corner_first_order_reconstruction &&
        cfg.outlet_corner_first_order_band_dp * cfg.dp >= 0.5 * cfg.DL)
    {
        throw std::runtime_error(
            "The outlet-corner first-order x band must be narrower than half the domain length.");
    }
    if (cfg.outlet_y_wall_corner_first_order_reconstruction &&
        cfg.y_wall_first_order_band_dp * cfg.dp >= 0.5 * cfg.DH)
    {
        throw std::runtime_error(
            "The outlet-corner first-order y band must be narrower than half the channel height.");
    }
    if (cfg.wall_flux_probe_interval <= 0)
    {
        throw std::runtime_error("wall_flux_probe_interval must be positive.");
    }
    if (cfg.muscl_hllc_dissipation_limiter &&
        (!std::isfinite(cfg.muscl_hllc_limiter_parameter) || cfg.muscl_hllc_limiter_parameter <= 0.0))
    {
        throw std::runtime_error(
            "muscl_hllc_limiter_parameter must be finite and positive when the dissipation limiter is enabled.");
    }
    if (cfg.positivity_floor &&
        (!std::isfinite(cfg.positivity_floor_factor) || cfg.positivity_floor_factor <= 0.0 ||
         cfg.positivity_floor_factor >= 1.0))
    {
        throw std::runtime_error(
            "positivity_floor_factor must be finite and lie in (0, 1) when the positivity floor is enabled.");
    }
    if (std::round(cfg.DL / cfg.dp) < 1.0 || std::round(cfg.DH / cfg.dp) < 1.0 ||
        std::round(cfg.DW / cfg.dp) < 1.0)
    {
        throw std::runtime_error(
            "Invalid geometry: each periodic-lattice direction needs at least one cell.");
    }
    const Real clearance_x = std::min(cfg.cylinder_center_x, cfg.DL - cfg.cylinder_center_x) - cfg.cylinder_radius;
    const Real clearance_y = std::min(cfg.cylinder_center_y, cfg.DH - cfg.cylinder_center_y) - cfg.cylinder_radius;
    if (clearance_x <= cfg.dp || clearance_y <= cfg.dp)
    {
        throw std::runtime_error("Invalid geometry: cylinder must stay inside the physical box with at least one dp clearance.");
    }

    // ---- freestream: gamma / rho_inf / c_inf are locked by shared code ----
    if (cfg.gamma <= 1.0)
    {
        throw std::runtime_error("Invalid physical parameters: gamma must be greater than 1.");
    }
    if (std::fabs(cfg.gamma - 1.4) > 1.0e-12)
    {
        throw std::runtime_error("Unsupported gamma: the shared compressible integration hardcodes CompressibleFluid(1.4), so gamma must be 1.4.");
    }
    if (cfg.rho_inf <= 0.0 || cfg.c_inf <= 0.0)
    {
        throw std::runtime_error("Invalid physical parameters: rho_inf and c_inf must be positive.");
    }
    if (std::fabs(cfg.rho_inf - 1.0) > 1.0e-12 || std::fabs(cfg.c_inf - 1.0) > 1.0e-12)
    {
        throw std::runtime_error("Unsupported reference state: CompressibleFluid::ReferenceDensity()/ReferenceSoundSpeed() return 1.0, so rho_inf and c_inf must both be 1.0.");
    }
    if (!std::isfinite(cfg.mach_inf) || cfg.mach_inf <= 0.0)
    {
        throw std::runtime_error("Invalid mach_inf: Mach number must be finite and positive.");
    }
    if (cfg.re <= 0.0)
    {
        throw std::runtime_error("Invalid physical parameters: re must be positive.");
    }

    cfg.u_inf = cfg.mach_inf * cfg.c_inf;
    cfg.p_inf = cfg.rho_inf * cfg.c_inf * cfg.c_inf / cfg.gamma;
    cfg.E_inf_per_volume = cfg.p_inf / (cfg.gamma - 1.0) + 0.5 * cfg.rho_inf * cfg.u_inf * cfg.u_inf;
    cfg.mu_inf = cfg.rho_inf * cfg.u_inf * cfg.D / cfg.re;

    if (!std::isfinite(cfg.p_inf) || !std::isfinite(cfg.E_inf_per_volume) ||
        !std::isfinite(cfg.mu_inf) || cfg.p_inf <= 0.0 ||
        cfg.E_inf_per_volume <= 0.0 || cfg.mu_inf <= 0.0)
    {
        throw std::runtime_error("Derived freestream is non-physical: p_inf, total energy density and mu must all be positive.");
    }

    // ---- simulation controls ----
    if (cfg.end_time <= 0.0 || cfg.output_interval <= 0 || cfg.screen_output_interval <= 0 ||
        cfg.acoustic_cfl <= 0.0 || cfg.boundary_n_layers <= 0 || cfg.smoke_min_steps <= 0)
    {
        throw std::runtime_error("Invalid simulation controls: time, intervals, CFL and layer counts must be positive.");
    }
    if (cfg.ghost_reserve_factor <= 0.0)
    {
        throw std::runtime_error("Invalid ghost_reserve_factor: ghost reserve must be positive.");
    }
    if (cfg.velocity_noise_amplitude < 0.0 || cfg.velocity_noise_amplitude >= 1.0)
    {
        throw std::runtime_error("Invalid velocity_noise_amplitude: must lie in [0, 1).");
    }
    if (cfg.mass_drift_tolerance <= 0.0 || cfg.energy_drift_tolerance <= 0.0 ||
        cfg.wall_flux_relative_tolerance <= 0.0 || cfg.wall_flux_relative_tolerance >= 1.0 ||
        cfg.wall_flux_absolute_floor_factor <= 0.0 || cfg.wall_flux_absolute_floor_factor >= 1.0)
    {
        throw std::runtime_error("Invalid validation tolerances: drift tolerances must be positive and the wall-flux relative tolerance / absolute floor factor must lie in (0, 1).");
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
    if (cfg.x_boundary_mode == XBoundaryMode::Periodic &&
        cfg.y_boundary_mode != YBoundaryMode::Wall)
    {
        throw std::runtime_error(
            "x_boundary=periodic currently requires y_boundary=wall: this case only validates the fully closed x/z-periodic channel topology.");
    }

    // The fluid fills the physical box only -- no sponge / padding layers.
    const Real nx = cfg.DL / cfg.dp;
    const Real ny = cfg.DH / cfg.dp;
    const Real nz = cfg.DW / cfg.dp;
    const Real estimated_particles = nx * ny * nz;
    if (estimated_particles > static_cast<Real>(cfg.max_local_particles))
    {
        throw std::runtime_error("Estimated fluid particle count exceeds max_local_particles; refusing local large case.");
    }
}

inline CompressibleFreestreamState makeFreestreamState(const CompressibleCylinderConfig &cfg)
{
    CompressibleFreestreamState state;
    state.gamma = cfg.gamma;
    state.rho = cfg.rho_inf;
    state.p = cfg.p_inf;
    state.sound_speed = cfg.c_inf;
    state.mach = cfg.mach_inf;
    state.vel = Vecd(cfg.u_inf, 0.0, 0.0);
    state.E_per_volume = cfg.E_inf_per_volume;
    return state;
}

inline CompressibleCylinderConfig loadConfig(const std::string &config_path)
{
    IniConfig ini = loadIniFile(config_path);
    CompressibleCylinderConfig cfg;

    // Reject weakly-compressible legacy keys instead of silently ignoring them:
    // a stale config must never be reinterpreted as a valid compressible setup.
    static const std::pair<const char *, const char *> rejected_keys[] = {
        {"physical", "rho0_f"},
        {"physical", "u_f"},
        {"physical", "sound_speed_factor"},
        {"physical", "outlet_pressure"},
        {"geometry", "sponge_width_factor"},
    };
    for (const auto &entry : rejected_keys)
    {
        if (hasKey(ini, entry.first, entry.second))
        {
            throw std::runtime_error(
                std::string("Unsupported legacy key '") + entry.second + "' in section [" + entry.first +
                "]: the compressible cylinder case is configured by gamma/rho_inf/c_inf/mach_inf/re only.");
        }
    }

    cfg.DL = parseReal(requireValue(ini, "geometry", "dl"));
    cfg.DH = parseReal(requireValue(ini, "geometry", "dh"));
    cfg.DW = parseReal(requireValue(ini, "geometry", "dw"));
    cfg.global_resolution = parseReal(requireValue(ini, "geometry", "global_resolution"));
    cfg.cylinder_center_x = parseReal(requireValue(ini, "geometry", "cylinder_center_x"));
    cfg.cylinder_center_y = parseReal(requireValue(ini, "geometry", "cylinder_center_y"));
    cfg.cylinder_radius = parseReal(requireValue(ini, "geometry", "cylinder_radius"));
    cfg.max_local_particles = optionalInt(ini, "geometry", "max_local_particles", 1000000);

    cfg.gamma = parseReal(requireValue(ini, "physical", "gamma"));
    cfg.rho_inf = parseReal(requireValue(ini, "physical", "rho_inf"));
    cfg.c_inf = parseReal(requireValue(ini, "physical", "c_inf"));
    cfg.mach_inf = parseReal(requireValue(ini, "physical", "mach_inf"));
    cfg.re = parseReal(requireValue(ini, "physical", "re"));

    cfg.end_time = parseReal(requireValue(ini, "simulation", "end_time"));
    cfg.output_interval = parseInt(requireValue(ini, "simulation", "output_interval"));
    cfg.screen_output_interval = parseInt(requireValue(ini, "simulation", "screen_output_interval"));
    cfg.acoustic_cfl = parseReal(requireValue(ini, "simulation", "acoustic_cfl"));
    cfg.boundary_n_layers = parseInt(requireValue(ini, "simulation", "boundary_n_layers"));
    cfg.ghost_reserve_factor = optionalReal(ini, "simulation", "ghost_reserve_factor", 0.5);
    cfg.enable_restart = optionalBool(ini, "simulation", "enable_restart", false);
    cfg.restart_step = optionalInt(ini, "simulation", "restart_step", 0);
    cfg.restart_output_factor = optionalInt(ini, "simulation", "restart_output_factor", 10);
    cfg.restart_keep_last_n = optionalInt(ini, "simulation", "restart_keep_last_n", 3);
    if (hasKey(ini, "boundary", "x_boundary"))
    {
        cfg.x_boundary_mode = parseXBoundaryMode(requireValue(ini, "boundary", "x_boundary"));
    }
    if (hasKey(ini, "boundary", "y_boundary"))
    {
        cfg.y_boundary_mode = parseYBoundaryMode(requireValue(ini, "boundary", "y_boundary"));
    }

    cfg.smoke_min_steps = optionalInt(ini, "validation", "smoke_min_steps", 50);
    cfg.velocity_noise_amplitude = optionalReal(ini, "validation", "velocity_noise_amplitude", 0.0);
    cfg.mass_drift_tolerance = optionalReal(ini, "validation", "mass_drift_tolerance", 1.0e-4);
    cfg.energy_drift_tolerance = optionalReal(ini, "validation", "energy_drift_tolerance", 1.0e-3);
    cfg.wall_flux_relative_tolerance =
        optionalReal(ini, "validation", "wall_flux_relative_tolerance", 1.0e-2);
    cfg.wall_flux_absolute_floor_factor =
        optionalReal(ini, "validation", "wall_flux_absolute_floor_factor", 1.0e-10);
    cfg.y_wall_first_order_reconstruction =
        optionalBool(ini, "validation", "y_wall_first_order_reconstruction", false);
    cfg.y_wall_first_order_band_dp =
        optionalReal(ini, "validation", "y_wall_first_order_band_dp", 2.6);
    cfg.outlet_y_wall_corner_first_order_reconstruction =
        optionalBool(ini, "validation", "outlet_y_wall_corner_first_order_reconstruction", false);
    cfg.outlet_corner_first_order_band_dp =
        optionalReal(ini, "validation", "outlet_corner_first_order_band_dp", 2.6);
    cfg.wall_contact_first_order_reconstruction =
        optionalBool(ini, "validation", "wall_contact_first_order_reconstruction", false);
    cfg.wall_flux_probe_interval =
        optionalInt(ini, "validation", "wall_flux_probe_interval", 50);
    cfg.muscl_hllc_dissipation_limiter =
        optionalBool(ini, "validation", "muscl_hllc_dissipation_limiter", false);
    cfg.muscl_hllc_limiter_parameter =
        optionalReal(ini, "validation", "muscl_hllc_limiter_parameter", 5.0);
    cfg.positivity_floor =
        optionalBool(ini, "validation", "positivity_floor", false);
    cfg.positivity_floor_factor =
        optionalReal(ini, "validation", "positivity_floor_factor", 1.0e-2);

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

inline Vecd freestreamVelocity(const CompressibleCylinderConfig &cfg)
{
    return Vecd(cfg.u_inf, 0.0, 0.0);
}

inline Real forceScale(const CompressibleCylinderConfig &cfg)
{
    return 0.5 * cfg.rho_inf * cfg.u_inf * cfg.u_inf * cfg.D * cfg.DW;
}

inline void printConfigSummary(const CompressibleCylinderConfig &cfg)
{
    const CompressibleFreestreamState freestream = makeFreestreamState(cfg);
    std::cout << "[Cylinder3DCompressible][Config]"
              << " gamma=" << cfg.gamma
              << " Ma=" << freestream.mach
              << " Re=" << cfg.re
              << " rho_inf=" << freestream.rho
              << " c_inf=" << freestream.sound_speed
              << " u_inf=" << cfg.u_inf
              << " p_inf=" << freestream.p
              << " mu=" << cfg.mu_inf
              << " E_inf/V=" << freestream.E_per_volume
              << " D=" << cfg.D
              << " dp=" << cfg.dp
              << " D/dp=" << cfg.D / cfg.dp
              << " DW=" << cfg.DW
              << " x_boundary=" << xBoundaryModeName(cfg.x_boundary_mode)
              << " y_boundary=" << yBoundaryModeName(cfg.y_boundary_mode)
              << " acoustic_cfl=" << cfg.acoustic_cfl
              << " y_wall_first_order_reconstruction=" << cfg.y_wall_first_order_reconstruction
              << " y_wall_first_order_band_dp=" << cfg.y_wall_first_order_band_dp
              << " outlet_y_wall_corner_first_order_reconstruction="
              << cfg.outlet_y_wall_corner_first_order_reconstruction
              << " outlet_corner_first_order_band_dp=" << cfg.outlet_corner_first_order_band_dp
              << " wall_contact_first_order_reconstruction="
              << cfg.wall_contact_first_order_reconstruction
              << " wall_flux_probe_interval=" << cfg.wall_flux_probe_interval
              << " muscl_hllc_dissipation_limiter=" << cfg.muscl_hllc_dissipation_limiter
              << " muscl_hllc_limiter_parameter=" << cfg.muscl_hllc_limiter_parameter
              << " positivity_floor=" << cfg.positivity_floor
              << " positivity_floor_factor=" << cfg.positivity_floor_factor << std::endl;
}

} // namespace cylinder_3d_compressible
} // namespace SPH

#endif // CYLINDER_3D_COMPRESSIBLE_DATA_HPP
