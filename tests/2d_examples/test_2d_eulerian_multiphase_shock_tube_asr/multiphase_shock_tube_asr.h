/**
 * @file 	multiphase_shock_tube_asr.h
 * @brief 	2D Eulerian SPH multiphase shock tube with particle-band adaptive
 *          spatial resolution (SPH-ASR, Yang/Kong/Liu PRE 104:055308).
 *          Kapila five-equation model, first-order Godunov, reflective wall
 *          ends, y-periodic. Configuration from an INI file selected on the
 *          command line (--config=name.ini, default config_twogas.ini).
 *
 *          Initial particle layouts:
 *          - graded: banded lattice, 5 columns per band with geometric
 *            spacing ds_k = ds_min C_r^k outward from the interface and
 *            row-stretched cells dy_k = H / round(H/ds_k) so that the y
 *            periodicity is preserved at every band;
 *          - uniform: single spacing lattice (USR baseline through the same
 *            code path and kernel).
 * @author 	KIYOYOZU
 */
#include "background_pressure_correction.h"
#include "eulerian_multiphase_riemann_solver.h"
#include "gradient_correction.h"
#include "particle_band_adaptation.h"
#include "particle_band_relation.h"
#include "particle_split_merge.h"
#include "shepard_density_filter.h"
#include "update_particle_bands.h"
#include "sphinxsys.h"
#include "time_step_local_h.h"
#include "update_smoothing_length_by_band.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
using namespace SPH;
//----------------------------------------------------------------------
//	Minimal INI parser (same convention as the validated shock tube case).
//----------------------------------------------------------------------
namespace asr_cfg_detail
{
using Ini = std::map<std::string, std::map<std::string, std::string>>;

inline std::string trim(const std::string &s)
{
    size_t b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos)
        return "";
    size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

inline Ini parse(const std::string &path)
{
    Ini ini;
    std::ifstream in(path);
    std::string section, line;
    while (std::getline(in, line))
    {
        std::string cleaned = trim(line);
        if (cleaned.empty() || cleaned[0] == '#' || cleaned[0] == ';')
            continue;
        if (cleaned.front() == '[' && cleaned.back() == ']')
        {
            section = trim(cleaned.substr(1, cleaned.size() - 2));
            std::transform(section.begin(), section.end(), section.begin(), ::tolower);
            continue;
        }
        size_t eq = cleaned.find('=');
        if (eq == std::string::npos || section.empty())
            continue;
        std::string key = trim(cleaned.substr(0, eq));
        std::string val = trim(cleaned.substr(eq + 1));
        size_t cmt = val.find_first_of("#;");
        if (cmt != std::string::npos)
            val = trim(val.substr(0, cmt));
        std::transform(key.begin(), key.end(), key.begin(), ::tolower);
        ini[section][key] = val;
    }
    return ini;
}

inline std::string lower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

inline Real real(const Ini &ini, const std::string &sec, const std::string &key, Real def)
{
    auto s = ini.find(lower(sec));
    if (s == ini.end())
        return def;
    auto k = s->second.find(lower(key));
    if (k == s->second.end())
        return def;
    try
    {
        size_t pos = 0;
        Real v = static_cast<Real>(std::stod(k->second, &pos));
        if (pos != k->second.size() || !std::isfinite(v))
            throw std::invalid_argument("bad number");
        return v;
    }
    catch (const std::exception &)
    {
        std::cerr << "[config] bad number for [" << sec << "]/" << key
                  << " = '" << k->second << "', using default " << def << "\n";
        return def;
    }
}

inline int integer(const Ini &ini, const std::string &sec, const std::string &key, int def)
{
    auto s = ini.find(lower(sec));
    if (s == ini.end())
        return def;
    auto k = s->second.find(lower(key));
    if (k == s->second.end())
        return def;
    try
    {
        size_t pos = 0;
        int v = std::stoi(k->second, &pos);
        if (pos != k->second.size())
            throw std::invalid_argument("bad integer");
        return v;
    }
    catch (const std::exception &)
    {
        std::cerr << "[config] bad integer for [" << sec << "]/" << key
                  << " = '" << k->second << "', using default " << def << "\n";
        return def;
    }
}

inline bool boolean(const Ini &ini, const std::string &sec, const std::string &key, bool def)
{
    auto s = ini.find(lower(sec));
    if (s == ini.end())
        return def;
    auto k = s->second.find(lower(key));
    if (k == s->second.end())
        return def;
    std::string v = lower(k->second);
    if (v == "on" || v == "true" || v == "1" || v == "yes")
        return true;
    if (v == "off" || v == "false" || v == "0" || v == "no")
        return false;
    std::cerr << "[config] bad boolean for [" << sec << "]/" << key
              << " = '" << k->second << "', using default " << def << "\n";
    return def;
}

inline std::string str(const Ini &ini, const std::string &sec, const std::string &key, const std::string &def)
{
    auto s = ini.find(lower(sec));
    if (s == ini.end())
        return def;
    auto k = s->second.find(lower(key));
    return (k == s->second.end()) ? def : k->second;
}

inline std::string resolvePath(const std::string &file_name)
{
    std::error_code ec;
#ifdef CASE_SOURCE_DIR
    std::filesystem::path src = std::filesystem::path(CASE_SOURCE_DIR) / file_name;
    if (std::filesystem::exists(src, ec) && !ec)
        return src.string();
#endif
    std::filesystem::path cwd = std::filesystem::current_path() / file_name;
    if (std::filesystem::exists(cwd, ec) && !ec)
        return cwd.string();
    return file_name;
}
} // namespace asr_cfg_detail
//----------------------------------------------------------------------
//	Configuration struct.
//----------------------------------------------------------------------
struct AsrConfig
{
    // geometry
    Real dp, L, x_interface;
    Real H;              /**< 20 dp, y-periodic slab height */
    Real wall_thickness; /**< 14 dp >= max kernel support 2*h_ref = 12 dp */
    // material
    Real gamma1, p_inf1, gamma2, p_inf2;
    // initial condition
    Real rho_L, p_L, u_L, alpha_L, rho_R, p_R, u_R, alpha_R;
    std::string preset; /**< "sod" or "quiescent" (uniform p = p_L, u = u_L) */
    // simulation
    Real end_time, output_interval, acoustic_cfl;
    int riemann_order;       /**< 1 = first-order, 2 = MUSCL second-order */
    std::string limiter;     /**< MUSCL slope limiter: minmod | mc | vanleer | none */
    std::string reconstruct; /**< MUSCL scope: "full" (rho/vel/p/alpha) | "vel_p" */
    // asr
    Real ds_max_factor;      /**< ds_max = factor * dp */
    Real band_coef;          /**< C_r, Eq. (34), 2^{1/2} in 2D */
    Real band_width_factor;  /**< dS_k = factor * ds_k, Eq. (38) */
    Real h_spacing_ratio;    /**< h_ref / ds_max, paper h_r = 1.5 V^{1/d} */
    Real gamma_split;        /**< split threshold gamma_s, Eq. (41) */
    Real gamma_merge;        /**< merge threshold gamma_m, Eq. (50) */
    Real split_lambda;       /**< split offset factor lambda, Eq. (45) */
    int adapt_interval;      /**< adaptation every N steps */
    bool initial_banding;    /**< graded banded lattice (false = uniform USR) */
    bool adapt_h;            /**< smoothing length evolution, Eqs. (6)-(9) */
    bool adapt_split_merge;  /**< particle split/merge, Eqs. (41)-(53) */
    bool adapt_bands;        /**< interface band tracking */
    Real buffer_growth;      /**< reserved buffer = growth * initial count */
    std::string shepard_filter; /**< "on" | "off" | "event_window" */
    int shepard_window_steps;   /**< window length after split/merge events */
    std::string kernel;      /**< "hyperbolic" (paper Eq. 5) | "wendland" */
    bool steady_correction;  /**< consistent-flux well-balanced correction on/off */
    bool corr_upwind_select; /**< correction: advective upwind endpoint at density jumps */
    Real corr_interface_density_ratio; /**< density jump ratio that triggers upwind selection */
    Real uniform_spacing_factor; /**< uniform-mode lattice ds = factor*dp (-1: follow ds_max) */
    int force_band;          /**< test hook: overwrite all bands (-1 = off) */
    std::string band_tracking; /**< "auto" (planar + Dijkstra fallback) | "dijkstra" */
    Real interface_alpha_tol;  /**< mixed-state window for Dijkstra band-0 seeds */
    Real band_hysteresis;      /**< band change margin in units of ds_target */
    bool shock_band;           /**< also track the strongest planar same-phase pressure jump */
    Real shock_rel_jump;       /**< min relative pressure jump accepted as a shock */
    bool check_symmetry;     /**< pair-list symmetry audit at t=0 */
    bool print_moment;       /**< first-moment diagnostic at t=0 */
    int screen_output_interval;
};

inline AsrConfig loadAsrConfig(const std::string &file_name)
{
    auto ini = asr_cfg_detail::parse(asr_cfg_detail::resolvePath(file_name));
    AsrConfig c;
    c.dp = asr_cfg_detail::real(ini, "geometry", "dp", 1.0 / 400.0);
    c.L = asr_cfg_detail::real(ini, "geometry", "L", 1.0);
    c.x_interface = asr_cfg_detail::real(ini, "geometry", "x_interface", 0.5);
    // H defaults to 20 dp; the USR-coarse baseline overrides it to keep the
    // same physical tube height when dp is coarsened.
    c.H = asr_cfg_detail::real(ini, "geometry", "H", 20.0 * c.dp);
    c.wall_thickness = asr_cfg_detail::real(ini, "geometry", "wall_thickness", 14.0 * c.dp);
    c.gamma1 = asr_cfg_detail::real(ini, "material", "gamma1", 1.4);
    c.p_inf1 = asr_cfg_detail::real(ini, "material", "p_inf1", 0.0);
    c.gamma2 = asr_cfg_detail::real(ini, "material", "gamma2", 1.6);
    c.p_inf2 = asr_cfg_detail::real(ini, "material", "p_inf2", 0.0);
    c.rho_L = asr_cfg_detail::real(ini, "ic", "rho_L", 1.0);
    c.p_L = asr_cfg_detail::real(ini, "ic", "p_L", 0.425);
    c.u_L = asr_cfg_detail::real(ini, "ic", "u_L", 0.0);
    c.alpha_L = asr_cfg_detail::real(ini, "ic", "alpha_L", 1.0);
    c.rho_R = asr_cfg_detail::real(ini, "ic", "rho_R", 0.125);
    c.p_R = asr_cfg_detail::real(ini, "ic", "p_R", 0.1);
    c.u_R = asr_cfg_detail::real(ini, "ic", "u_R", 0.0);
    c.alpha_R = asr_cfg_detail::real(ini, "ic", "alpha_R", 0.0);
    c.preset = asr_cfg_detail::lower(
        asr_cfg_detail::str(ini, "ic", "preset", "sod"));
    c.end_time = asr_cfg_detail::real(ini, "simulation", "end_time", 0.2);
    c.output_interval = asr_cfg_detail::real(ini, "simulation", "output_interval", 0.02);
    c.acoustic_cfl = asr_cfg_detail::real(ini, "simulation", "acoustic_cfl", 0.1);
    c.riemann_order = asr_cfg_detail::integer(ini, "simulation", "riemann_order", 1);
    c.limiter = asr_cfg_detail::lower(
        asr_cfg_detail::str(ini, "simulation", "muscl_limiter", "minmod"));
    c.reconstruct = asr_cfg_detail::lower(
        asr_cfg_detail::str(ini, "simulation", "muscl_reconstruct", "full"));
    c.ds_max_factor = asr_cfg_detail::real(ini, "asr", "ds_max_factor", 4.0);
    c.band_coef = asr_cfg_detail::real(ini, "asr", "band_coef", std::sqrt(2.0));
    c.band_width_factor = asr_cfg_detail::real(ini, "asr", "band_width_factor", 5.0);
    c.h_spacing_ratio = asr_cfg_detail::real(ini, "asr", "h_spacing_ratio", 1.5);
    c.gamma_split = asr_cfg_detail::real(ini, "asr", "gamma_split", 1.5);
    c.gamma_merge = asr_cfg_detail::real(ini, "asr", "gamma_merge", 0.7);
    c.split_lambda = asr_cfg_detail::real(ini, "asr", "split_lambda", 0.6);
    c.adapt_interval = asr_cfg_detail::integer(ini, "asr", "adapt_interval", 5);
    c.initial_banding = asr_cfg_detail::boolean(ini, "asr", "initial_banding", true);
    c.adapt_h = asr_cfg_detail::boolean(ini, "asr", "adapt_h", false);
    c.adapt_split_merge = asr_cfg_detail::boolean(ini, "asr", "adapt_split_merge", false);
    c.adapt_bands = asr_cfg_detail::boolean(ini, "asr", "adapt_bands", false);
    c.buffer_growth = asr_cfg_detail::real(ini, "asr", "buffer_growth", 1.0);
    c.shepard_filter = asr_cfg_detail::lower(
        asr_cfg_detail::str(ini, "asr", "shepard_filter", "off"));
    c.shepard_window_steps = asr_cfg_detail::integer(ini, "asr", "shepard_window_steps", 50);
    c.kernel = asr_cfg_detail::lower(
        asr_cfg_detail::str(ini, "asr", "kernel", "hyperbolic"));
    c.steady_correction = asr_cfg_detail::boolean(ini, "asr", "steady_correction", false);
    c.corr_upwind_select = asr_cfg_detail::boolean(ini, "asr", "corr_upwind_select", true);
    c.corr_interface_density_ratio =
        asr_cfg_detail::real(ini, "asr", "corr_interface_density_ratio", 2.0);
    c.uniform_spacing_factor = asr_cfg_detail::real(ini, "asr", "uniform_spacing_factor", -1.0);
    c.force_band = asr_cfg_detail::integer(ini, "asr", "force_band", -1);
    c.band_tracking = asr_cfg_detail::lower(
        asr_cfg_detail::str(ini, "asr", "band_tracking", "auto"));
    c.interface_alpha_tol = asr_cfg_detail::real(ini, "asr", "interface_alpha_tol", 0.1);
    c.band_hysteresis = asr_cfg_detail::real(ini, "asr", "band_hysteresis", 0.5);
    c.shock_band = asr_cfg_detail::boolean(ini, "asr", "shock_band", false);
    c.shock_rel_jump = asr_cfg_detail::real(ini, "asr", "shock_rel_jump", 0.01);
    c.check_symmetry = asr_cfg_detail::boolean(ini, "asr", "check_symmetry", false);
    c.print_moment = asr_cfg_detail::boolean(ini, "asr", "print_moment", false);
    c.screen_output_interval = asr_cfg_detail::integer(ini, "simulation", "screen_output_interval", 100);

    if (c.preset != "sod" && c.preset != "quiescent")
    {
        std::cerr << "[config] unknown ic/preset '" << c.preset << "', using sod.\n";
        c.preset = "sod";
    }
    if (c.shepard_filter != "on" && c.shepard_filter != "off" &&
        c.shepard_filter != "event_window")
    {
        std::cerr << "[config] unknown asr/shepard_filter '" << c.shepard_filter
                  << "', using off.\n";
        c.shepard_filter = "off";
    }
    if (c.kernel != "hyperbolic" && c.kernel != "wendland")
    {
        std::cerr << "[config] unknown asr/kernel '" << c.kernel << "', using hyperbolic.\n";
        c.kernel = "hyperbolic";
    }
    if (c.band_tracking != "auto" && c.band_tracking != "dijkstra")
    {
        std::cerr << "[config] unknown asr/band_tracking '" << c.band_tracking
                  << "', using auto.\n";
        c.band_tracking = "auto";
    }
    if (c.riemann_order != 1 && c.riemann_order != 2)
    {
        std::cerr << "[config] unknown simulation/riemann_order " << c.riemann_order
                  << ", using 1.\n";
        c.riemann_order = 1;
    }
    if (c.limiter != "minmod" && c.limiter != "mc" && c.limiter != "vanleer" &&
        c.limiter != "none")
    {
        std::cerr << "[config] unknown simulation/muscl_limiter '" << c.limiter
                  << "', using minmod.\n";
        c.limiter = "minmod";
    }
    if (c.reconstruct != "full" && c.reconstruct != "vel_p")
    {
        std::cerr << "[config] unknown simulation/muscl_reconstruct '" << c.reconstruct
                  << "', using full.\n";
        c.reconstruct = "full";
    }
    return c;
}
//----------------------------------------------------------------------
//	Fluid block and wall shapes (H and wall thickness follow the ASR
//	geometry decision: H = 20 dp, wall = 14 dp >= 2 h_ref = 12 dp).
//----------------------------------------------------------------------
class AsrFluidBlock : public ComplexShape
{
  public:
    explicit AsrFluidBlock(const std::string &shape_name, Real L, Real H) : ComplexShape(shape_name)
    {
        std::vector<Vecd> rect_shape;
        rect_shape.push_back(Vecd(0.0, 0.0));
        rect_shape.push_back(Vecd(L, 0.0));
        rect_shape.push_back(Vecd(L, H));
        rect_shape.push_back(Vecd(0.0, H));
        rect_shape.push_back(Vecd(0.0, 0.0));
        MultiPolygon polygon(rect_shape);
        add<MultiPolygonShape>(polygon, "OuterBoundary");
    }
};

class AsrWallBlock : public MultiPolygonShape
{
  public:
    AsrWallBlock(const std::string &shape_name, Real x0, Real x1, Real H)
        : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> rect_shape;
        rect_shape.push_back(Vecd(x0, 0.0));
        rect_shape.push_back(Vecd(x1, 0.0));
        rect_shape.push_back(Vecd(x1, H));
        rect_shape.push_back(Vecd(x0, H));
        rect_shape.push_back(Vecd(x0, 0.0));
        multi_polygon_.addPolygon(rect_shape, GeometricOps::add);
    }
};
//----------------------------------------------------------------------
//	Banded lattice generator: 5 columns per band (band width 5 ds_k,
//	Eq. 38) marching outward from the interface with geometric spacing
//	ds_k, rows stretched to dy_k = H / round(H / ds_k) for y periodicity,
//	Vol = ds_k * dy_k.
//----------------------------------------------------------------------
class BandedLattice;

template <>
class ParticleGenerator<BaseParticles, BandedLattice> : public ParticleGenerator<BaseParticles>
{
  public:
    ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles, const AsrConfig &cfg)
        : ParticleGenerator<BaseParticles>(sph_body, base_particles), cfg_(cfg) {};

  protected:
    AsrConfig cfg_;

    void prepareGeometricData() override
    {
        ParticleBandAdaptation &adaptation =
            DynamicCast<ParticleBandAdaptation>(this, sph_body_.getSPHAdaptation());
        const Real L = cfg_.L;
        const Real H = cfg_.H;
        const Real x_if = cfg_.x_interface;
        const int band_max = adaptation.BandCount();
        const int columns_per_band = (int)std::round(cfg_.band_width_factor);

        auto add_column = [&](Real x, int k, Real ds_override = -1.0)
        {
            Real ds = ds_override > 0.0 ? ds_override : adaptation.BandSpacing(k);
            int n_rows = std::max(1, (int)std::round(H / ds));
            Real dy = H / Real(n_rows);
            Real Vol = ds * dy;
            for (int r = 0; r != n_rows; ++r)
                addPositionAndVolumetricMeasure(Vecd(x, (Real(r) + 0.5) * dy), Vol);
        };

        for (Real dir : {1.0, -1.0})
        {
            const Real x_edge = dir > 0.0 ? L : 0.0;
            Real x_cursor = x_if;

            // stage 1: graded bands 0..band_max-1 with exact spacings ds_k
            if (cfg_.initial_banding)
            {
                int m = 0;
                for (int k = 0; k < band_max; ++k)
                {
                    Real ds = adaptation.BandSpacing(k);
                    for (m = 0; m != columns_per_band; ++m)
                    {
                        if (dir * (x_edge - x_cursor) < ds)
                            break;
                        add_column(x_cursor + dir * 0.5 * ds, k);
                        x_cursor += dir * ds;
                    }
                    if (dir * (x_edge - x_cursor) < adaptation.BandSpacing(k))
                        break;
                }
            }

            // stage 2: fill the remaining span to the wall with columns on a
            // wall-symmetric phase (centers at (m + 1/2) s from the edge,
            // spacing stretched to divide the span exactly). The mirror wall
            // then sits on the lattice continuation sites, so the combined
            // inner + wall stencil is centrosymmetric about every fluid
            // particle and sum(Vol) = L*H holds exactly.
            Real span = dir * (x_edge - x_cursor);
            Real ds_coarse = adaptation.BandSpacing(band_max);
            // test hook: uniform lattice spacing decoupled from the band
            // spacing so that split/merge events can be driven by force_band
            if (!cfg_.initial_banding && cfg_.uniform_spacing_factor > 0.0)
                ds_coarse = cfg_.uniform_spacing_factor * cfg_.dp;
            if (span > 0.25 * ds_coarse)
            {
                int n_cols = std::max(1, (int)std::round(span / ds_coarse));
                Real s = span / Real(n_cols);
                for (int m = 0; m != n_cols; ++m)
                    add_column(x_edge - dir * (Real(m) + 0.5) * s, band_max, s);
            }
        }
    }
};
//----------------------------------------------------------------------
//	Mirror-wall lattice: wall particles are the reflection of the fluid
//	lattice across the wall face, with identical volumes. This makes the
//	wall stencil the exact negative image of the missing fluid half-space,
//	so a quiescent uniform-pressure state stays at rest even for coarse
//	particles whose 4 dp spacing does not match a dp-spaced wall (the
//	quadrature mismatch otherwise drives a spurious wall force). The wall
//	slab only carries reflected copies of fluid columns whose mirror falls
//	inside it; support coverage holds because the slab (14 dp) is thicker
//	than the distance any in-slab mirror can matter (kernel support of a
//	particle at distance d from the face reaches only d + 12 dp into the
//	wall, and mirrors beyond the slab belong to fluid columns farther than
//	the support from the face).
//----------------------------------------------------------------------
class MirrorWallLattice;

template <>
class ParticleGenerator<BaseParticles, MirrorWallLattice> : public ParticleGenerator<BaseParticles>
{
  public:
    ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles,
                      BaseParticles &fluid_particles, Real wall_face, Real inward, Real thickness)
        : ParticleGenerator<BaseParticles>(sph_body, base_particles),
          fluid_particles_(fluid_particles), wall_face_(wall_face),
          inward_(inward), thickness_(thickness) {};

  protected:
    BaseParticles &fluid_particles_;
    Real wall_face_;     /**< x coordinate of the fluid-side wall face */
    Real inward_;        /**< +1 (left wall) or -1 (right wall), towards fluid */
    Real thickness_;     /**< wall slab thickness */

    void prepareGeometricData() override
    {
        Vecd *fluid_pos = fluid_particles_.getVariableDataByName<Vecd>("Position");
        Real *fluid_Vol = fluid_particles_.getVariableDataByName<Real>("VolumetricMeasure");
        for (size_t i = 0; i != fluid_particles_.TotalRealParticles(); ++i)
        {
            Real depth = inward_ * (fluid_pos[i][0] - wall_face_);
            if (depth <= 0.0 || depth > thickness_)
                continue;
            Vecd mirror_pos = fluid_pos[i];
            mirror_pos[0] = wall_face_ - inward_ * depth;
            addPositionAndVolumetricMeasure(mirror_pos, fluid_Vol[i]);
        }
    }
};
//----------------------------------------------------------------------
//	Constant analytic wall normal (points from the solid into the fluid).
//----------------------------------------------------------------------
class AsrWallNormal : public LocalDynamics
{
  public:
    AsrWallNormal(SPHBody &sph_body, const Vecd &normal)
        : LocalDynamics(sph_body), normal_(normal),
          n_(particles_->registerStateVariableData<Vecd>("NormalDirection")) {};

    void update(size_t index_i, Real dt = 0.0)
    {
        n_[index_i] = normal_;
    }

  protected:
    Vecd normal_;
    Vecd *n_;
};
//----------------------------------------------------------------------
//	Initial condition. "sod": two-state Riemann data. "quiescent":
//	uniform pressure p_L and velocity u_L with the two-phase density and
//	alpha jump at the interface (exact equilibrium of the five-equation
//	model, used by the static and translating-interface gates).
//----------------------------------------------------------------------
class ShockTubeAsrInitialCondition : public LocalDynamics
{
  public:
    explicit ShockTubeAsrInitialCondition(SPHBody &sph_body, MultiphaseMixture &mixture,
                                          const AsrConfig &cfg)
        : LocalDynamics(sph_body), mixture_(mixture), cfg_(cfg),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          rho_(particles_->getVariableDataByName<Real>("Density")),
          mass_(particles_->getVariableDataByName<Real>("Mass")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
          p_(particles_->registerStateVariableData<Real>("Pressure")),
          vel_(particles_->registerStateVariableData<Vecd>("Velocity")),
          mom_(particles_->registerStateVariableData<Vecd>("Momentum")),
          E_(particles_->registerStateVariableData<Real>("TotalEnergy")),
          alpha_(particles_->registerStateVariableData<Real>("VolumeFraction")) {};

    void update(size_t index_i, Real dt)
    {
        const Real x = pos_[index_i][0];
        const bool left = x < cfg_.x_interface;
        const bool quiescent = cfg_.preset == "quiescent";

        rho_[index_i] = left ? cfg_.rho_L : cfg_.rho_R;
        p_[index_i] = quiescent ? cfg_.p_L : (left ? cfg_.p_L : cfg_.p_R);
        alpha_[index_i] = left ? cfg_.alpha_L : cfg_.alpha_R;
        vel_[index_i][0] = quiescent ? cfg_.u_L : (left ? cfg_.u_L : cfg_.u_R);
        vel_[index_i][1] = 0.0;

        Real rho_e_int = mixture_.MixtureInternalEnergyPerVolume(alpha_[index_i], p_[index_i]);
        mass_[index_i] = rho_[index_i] * Vol_[index_i];
        mom_[index_i] = mass_[index_i] * vel_[index_i];
        E_[index_i] = rho_e_int * Vol_[index_i]
                    + 0.5 * mass_[index_i] * vel_[index_i].squaredNorm();
    }

  protected:
    MultiphaseMixture &mixture_;
    AsrConfig cfg_;
    Vecd *pos_;
    Real *rho_, *mass_, *Vol_, *p_;
    Vecd *vel_, *mom_;
    Real *E_, *alpha_;
};
