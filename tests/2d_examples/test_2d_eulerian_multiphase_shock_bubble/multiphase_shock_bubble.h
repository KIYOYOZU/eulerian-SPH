/**
 * @file 	multiphase_shock_bubble.h
 * @brief 	2D Eulerian SPH shock-bubble interaction: a left-running water
 *          shock hits a circular gas bubble in water. Kapila five-equation
 *          model with stiffened gas EOS; all four sides are real SolidBody
 *          reflective walls (ContactRelation + five-equation WithWall
 *          integrators).
 *
 *          Reference setup (see case task_plan.md): domain window
 *          0.027 x 0.012 m, bubble center (0.012, 0.006) R = 0.003 m, shock
 *          plane x = 0.0201 with the post-shock water on its right moving
 *          left at -681.58 m/s. The right reservoir is extended to x_right
 *          (config) so the right-wall rarefaction stays outside the window
 *          of interest; within [0, 0.027] x [0, end_time] the setup is
 *          equivalent to a transmissive right boundary.
 *
 *          All parameters come from config.ini in the case source directory
 *          (defaults match the shipped file).
 * @author 	KIYOYOZU
 */
#include "sphinxsys.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
using namespace SPH;
//----------------------------------------------------------------------
//	Minimal INI parser (sections + key=value, '#'/';' comments, both
//	full-line and inline), same convention as the multiphase shock tube.
//----------------------------------------------------------------------
namespace sb_cfg_detail
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
        size_t cmt = val.find_first_of("#;"); // strip inline comments
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
    auto k = s->second.find(lower(key)); // keys are stored lowercased by parse()
    if (k == s->second.end())
        return def;
    try
    {
        size_t pos = 0;
        Real v = static_cast<Real>(std::stod(k->second, &pos));
        // Reject trailing garbage ("1.5abc") and non-finite values ("nan"/"inf").
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
        if (pos != k->second.size()) // reject trailing garbage ("2.5", "1abc")
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

inline std::string str(const Ini &ini, const std::string &sec, const std::string &key, const std::string &def)
{
    auto s = ini.find(lower(sec));
    if (s == ini.end())
        return def;
    auto k = s->second.find(lower(key));
    return (k == s->second.end()) ? def : k->second;
}

inline std::string resolvePath()
{
    std::error_code ec;
#ifdef CASE_SOURCE_DIR
    std::filesystem::path src = std::filesystem::path(CASE_SOURCE_DIR) / "config.ini";
    if (std::filesystem::exists(src, ec) && !ec)
        return src.string();
#endif
    std::filesystem::path cwd = std::filesystem::current_path() / "config.ini";
    if (std::filesystem::exists(cwd, ec) && !ec)
        return cwd.string();
    return "config.ini";
}
} // namespace sb_cfg_detail

struct ShockBubbleConfig
{
    Real dp, L, H, x_right;
    Real bubble_x, bubble_y, bubble_r, x_shock;
    Real gamma1, p_inf1, gamma2, p_inf2;
    Real rho_ps, p_ps, u_ps;   /**< post-shock water. */
    Real rho_0, p_0;           /**< pre-shock water. */
    Real rho_b, p_b;           /**< bubble gas. */
    Real end_time, output_interval, acoustic_cfl;
    int riemann_order;         // 1 = first-order, 2 = MUSCL second-order
    std::string limiter;
    std::string reconstruct;   // "full" or "vel_p" (MUSCL scope)
};

inline ShockBubbleConfig loadShockBubbleConfig()
{
    auto ini = sb_cfg_detail::parse(sb_cfg_detail::resolvePath());
    ShockBubbleConfig c;
    c.dp = sb_cfg_detail::real(ini, "geometry", "dp", 1.0e-4);
    c.L = sb_cfg_detail::real(ini, "geometry", "L", 0.027);
    c.H = sb_cfg_detail::real(ini, "geometry", "H", 0.012);
    c.x_right = sb_cfg_detail::real(ini, "geometry", "x_right", 0.05);
    c.bubble_x = sb_cfg_detail::real(ini, "geometry", "bubble_x", 0.012);
    c.bubble_y = sb_cfg_detail::real(ini, "geometry", "bubble_y", 0.006);
    c.bubble_r = sb_cfg_detail::real(ini, "geometry", "bubble_r", 0.003);
    c.x_shock = sb_cfg_detail::real(ini, "geometry", "x_shock", 0.0201);
    c.gamma1 = sb_cfg_detail::real(ini, "material", "gamma1", 1.4);
    c.p_inf1 = sb_cfg_detail::real(ini, "material", "p_inf1", 0.0);
    c.gamma2 = sb_cfg_detail::real(ini, "material", "gamma2", 4.4);
    c.p_inf2 = sb_cfg_detail::real(ini, "material", "p_inf2", 6.0e8);
    c.rho_ps = sb_cfg_detail::real(ini, "ic", "rho_postshock", 1323.65);
    c.p_ps = sb_cfg_detail::real(ini, "ic", "p_postshock", 1.9e9);
    c.u_ps = sb_cfg_detail::real(ini, "ic", "u_postshock", -681.58);
    c.rho_0 = sb_cfg_detail::real(ini, "ic", "rho_preshock", 1000.0);
    c.p_0 = sb_cfg_detail::real(ini, "ic", "p_preshock", 1.0e5);
    c.rho_b = sb_cfg_detail::real(ini, "ic", "rho_bubble", 1.0);
    c.p_b = sb_cfg_detail::real(ini, "ic", "p_bubble", 1.0e5);
    c.end_time = sb_cfg_detail::real(ini, "simulation", "end_time", 6.0e-6);
    c.output_interval = sb_cfg_detail::real(ini, "simulation", "output_interval", 0.25e-6);
    c.acoustic_cfl = sb_cfg_detail::real(ini, "simulation", "acoustic_cfl", 0.1);
    c.riemann_order = sb_cfg_detail::integer(ini, "simulation", "riemann_order", 2);
    c.limiter = sb_cfg_detail::str(ini, "simulation", "muscl_limiter", "minmod");
    c.reconstruct = sb_cfg_detail::str(ini, "simulation", "muscl_reconstruct", "vel_p");
    // Sanity warnings: silently building a non-physical layout is worse than
    // a loud warning at startup.
    if (c.x_right < c.L)
        std::cerr << "[config] WARNING: x_right < L; the reference window is not fully simulated.\n";
    if (!(c.bubble_x - c.bubble_r > 0.0 && c.bubble_x + c.bubble_r < c.x_shock))
        std::cerr << "[config] WARNING: bubble not strictly between the left wall and the shock plane.\n";
    if (!(c.x_shock < c.x_right))
        std::cerr << "[config] WARNING: x_shock not inside [0, x_right].\n";
    return c;
}

const ShockBubbleConfig SB_CFG = loadShockBubbleConfig();
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup (from config.ini).
//----------------------------------------------------------------------
Real particle_spacing_ref = SB_CFG.dp;
Real BW = 4.0 * particle_spacing_ref;     /**< System-domain padding around the box. */
Real L = SB_CFG.L;                        /**< Reference window length in x (figure). */
Real H = SB_CFG.H;                        /**< Box height in y. */
Real x_right = SB_CFG.x_right;            /**< Extended right reservoir end. */
Real wall_thickness = 3.0 * particle_spacing_ref; /**< Solid wall slab thickness. */
Real bubble_x = SB_CFG.bubble_x;
Real bubble_y = SB_CFG.bubble_y;
Real bubble_r = SB_CFG.bubble_r;
Real x_shock = SB_CFG.x_shock;
BoundingBoxd system_domain_bounds(Vec2d(-wall_thickness - BW, -wall_thickness - BW),
                                  Vec2d(x_right + wall_thickness + BW, H + wall_thickness + BW));
//----------------------------------------------------------------------
//	Material properties (stiffened gas EOS), from config.ini.
//----------------------------------------------------------------------
Real gamma_gas = SB_CFG.gamma1;
Real p_inf_gas = SB_CFG.p_inf1;
Real gamma_water = SB_CFG.gamma2;
Real p_inf_water = SB_CFG.p_inf2;
//----------------------------------------------------------------------
//	Initial condition, from config.ini.
//----------------------------------------------------------------------
Real end_time = SB_CFG.end_time;
//----------------------------------------------------------------------
//	Fluid body shape: rectangle [0, x_right] x [0, H].
//----------------------------------------------------------------------
class FluidBlock : public ComplexShape
{
  public:
    explicit FluidBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        std::vector<Vecd> rect_shape;
        rect_shape.push_back(Vecd(0.0, 0.0));
        rect_shape.push_back(Vecd(x_right, 0.0));
        rect_shape.push_back(Vecd(x_right, H));
        rect_shape.push_back(Vecd(0.0, H));
        rect_shape.push_back(Vecd(0.0, 0.0));
        MultiPolygon polygon(rect_shape);
        add<MultiPolygonShape>(polygon, "OuterBoundary");
    }
};
//----------------------------------------------------------------------
//	Solid wall slab [x0, x1] x [y0, y1] (bottom / top / left / right walls).
//----------------------------------------------------------------------
class WallBlock : public MultiPolygonShape
{
  public:
    WallBlock(const std::string &shape_name, Real x0, Real y0, Real x1, Real y1)
        : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> rect_shape;
        rect_shape.push_back(Vecd(x0, y0));
        rect_shape.push_back(Vecd(x1, y0));
        rect_shape.push_back(Vecd(x1, y1));
        rect_shape.push_back(Vecd(x0, y1));
        rect_shape.push_back(Vecd(x0, y0));
        multi_polygon_.addPolygon(rect_shape, GeometricOps::add);
    }
};
//----------------------------------------------------------------------
//	Constant analytic wall normal (points from the solid into the fluid).
//----------------------------------------------------------------------
class WallNormal : public LocalDynamics
{
  public:
    WallNormal(SPHBody &sph_body, const Vecd &normal)
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
//	Case-dependent initial condition for the five-equation model:
//	bubble gas inside the circle, pre-shock water left of the shock plane,
//	post-shock water right of it. Sets density, pressure, velocity, volume
//	fraction, mass, momentum and total energy (using the mixture EOS).
//----------------------------------------------------------------------
class ShockBubbleInitialCondition : public LocalDynamics
{
  public:
    explicit ShockBubbleInitialCondition(SPHBody &sph_body, MultiphaseMixture &mixture)
        : LocalDynamics(sph_body), mixture_(mixture),
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
        const Real y = pos_[index_i][1];
        Real rho_e_int; /**< Internal energy per unit volume. */

        if ((x - bubble_x) * (x - bubble_x) + (y - bubble_y) * (y - bubble_y) < bubble_r * bubble_r)
        {
            rho_[index_i] = SB_CFG.rho_b;
            p_[index_i] = SB_CFG.p_b;
            vel_[index_i][0] = 0.0;
            alpha_[index_i] = 1.0;
        }
        else if (x < x_shock)
        {
            rho_[index_i] = SB_CFG.rho_0;
            p_[index_i] = SB_CFG.p_0;
            vel_[index_i][0] = 0.0;
            alpha_[index_i] = 0.0;
        }
        else
        {
            rho_[index_i] = SB_CFG.rho_ps;
            p_[index_i] = SB_CFG.p_ps;
            vel_[index_i][0] = SB_CFG.u_ps;
            alpha_[index_i] = 0.0;
        }
        vel_[index_i][1] = 0.0;

        // Mixture internal energy per unit volume from pressure.
        rho_e_int = mixture_.MixtureInternalEnergyPerVolume(alpha_[index_i], p_[index_i]);

        mass_[index_i] = rho_[index_i] * Vol_[index_i];
        mom_[index_i] = mass_[index_i] * vel_[index_i];
        // Total energy = internal + kinetic.
        E_[index_i] = rho_e_int * Vol_[index_i]
                    + 0.5 * mass_[index_i] * vel_[index_i].squaredNorm();
    }

  protected:
    MultiphaseMixture &mixture_;
    Vecd *pos_;
    Real *rho_, *mass_, *Vol_, *p_;
    Vecd *vel_, *mom_;
    Real *E_, *alpha_;
};
