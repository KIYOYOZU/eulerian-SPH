/**
 * @file 	multiphase_shock_tube.h
 * @brief 	2D Eulerian SPH multiphase shock tube, Kapila five-equation model
 *          with stiffened gas EOS. x ends are real SolidBody reflective walls;
 *          y direction is periodic.
 *
 *          All parameters come from config.ini in the case source directory
 *          (see loadMultiphaseConfig below; defaults match the shipped file).
 *          The default setup is a two-material Sod shock tube: two DIFFERENT
 *          ideal gases (gamma1=1.4, gamma2=1.6, both p_inf=0) with the Sod
 *          left/right states, alpha=1 | alpha=0 across the membrane at
 *          x_interface. The impedance ratio is moderate, so the solution shows
 *          the classic rarefaction-fan + contact + shock shape while the alpha
 *          interface rides on the contact. This is a genuine two-phase
 *          (two-material) Riemann problem; the reference is the two-phase
 *          exact solution. [simulation]/riemann_order selects the Riemann
 *          interface treatment: 1 = first-order Godunov (validated default),
 *          2 = MUSCL second-order (experimental).
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
//	full-line and inline), mirroring the 3D cylinder case's config.ini
//	convention. Malformed numeric values fall back to the default with a
//	warning instead of throwing during static initialization.
//----------------------------------------------------------------------
namespace mp_cfg_detail
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
    // cwd first: each case runs from its own folder (cases/<name>/), so the
    // config and outputs stay per-case; the source dir is the fallback for
    // running the executable from the build tree directly
    std::error_code ec;
    std::filesystem::path cwd = std::filesystem::current_path() / "config.ini";
    if (std::filesystem::exists(cwd, ec) && !ec)
        return cwd.string();
#ifdef CASE_SOURCE_DIR
    std::filesystem::path src = std::filesystem::path(CASE_SOURCE_DIR) / "config.ini";
    if (std::filesystem::exists(src, ec) && !ec)
        return src.string();
#endif
    return "config.ini";
}
} // namespace mp_cfg_detail

struct MultiphaseConfig
{
    Real dp, L, x_interface;
    Real gamma1, p_inf1, gamma2, p_inf2;
    Real rho_L, p_L, u_L, alpha_L, rho_R, p_R, u_R, alpha_R;
    Real end_time, output_interval, acoustic_cfl;
    int riemann_order; // 1 = first-order, 2 = MUSCL second-order
    std::string limiter;
    std::string reconstruct; // "full" or "vel_p" (MUSCL scope)
};

inline MultiphaseConfig loadMultiphaseConfig()
{
    auto ini = mp_cfg_detail::parse(mp_cfg_detail::resolvePath());
    MultiphaseConfig c;
    c.dp = mp_cfg_detail::real(ini, "geometry", "dp", 1.0 / 400.0);
    c.L = mp_cfg_detail::real(ini, "geometry", "L", 1.0);
    c.x_interface = mp_cfg_detail::real(ini, "geometry", "x_interface", 0.5);
    c.gamma1 = mp_cfg_detail::real(ini, "material", "gamma1", 1.4);
    c.p_inf1 = mp_cfg_detail::real(ini, "material", "p_inf1", 0.0);
    c.gamma2 = mp_cfg_detail::real(ini, "material", "gamma2", 1.6);
    c.p_inf2 = mp_cfg_detail::real(ini, "material", "p_inf2", 0.0);
    c.rho_L = mp_cfg_detail::real(ini, "ic", "rho_L", 1.0);
    c.p_L = mp_cfg_detail::real(ini, "ic", "p_L", 0.425);
    c.u_L = mp_cfg_detail::real(ini, "ic", "u_L", 0.0);
    c.alpha_L = mp_cfg_detail::real(ini, "ic", "alpha_L", 1.0);
    c.rho_R = mp_cfg_detail::real(ini, "ic", "rho_R", 0.125);
    c.p_R = mp_cfg_detail::real(ini, "ic", "p_R", 0.1);
    c.u_R = mp_cfg_detail::real(ini, "ic", "u_R", 0.0);
    c.alpha_R = mp_cfg_detail::real(ini, "ic", "alpha_R", 0.0);
    c.end_time = mp_cfg_detail::real(ini, "simulation", "end_time", 0.2);
    c.output_interval = mp_cfg_detail::real(ini, "simulation", "output_interval", 0.02);
    c.acoustic_cfl = mp_cfg_detail::real(ini, "simulation", "acoustic_cfl", 0.1);
    c.riemann_order = mp_cfg_detail::integer(ini, "simulation", "riemann_order", 1);
    c.limiter = mp_cfg_detail::str(ini, "simulation", "muscl_limiter", "minmod");
    c.reconstruct = mp_cfg_detail::str(ini, "simulation", "muscl_reconstruct", "full");
    return c;
}

const MultiphaseConfig MP_CFG = loadMultiphaseConfig();
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup (from config.ini).
//----------------------------------------------------------------------
Real particle_spacing_ref = MP_CFG.dp;
Real BW = 4.0 * particle_spacing_ref;    /**< System-domain padding around the slab. */
Real L = MP_CFG.L;                       /**< Tube length in x. */
Real H = 10.0 * particle_spacing_ref;    /**< Tube height in y (periodic slab). */
Real wall_thickness = 3.0 * particle_spacing_ref; /**< Solid wall slab thickness. */
Real x_interface = MP_CFG.x_interface;   /**< Membrane / interface position. */
BoundingBoxd system_domain_bounds(Vec2d(-wall_thickness - BW, -BW),
                                  Vec2d(L + wall_thickness + BW, H + BW));
//----------------------------------------------------------------------
//	Material properties (stiffened gas EOS), from config.ini.
//----------------------------------------------------------------------
Real gamma_gas = MP_CFG.gamma1;
Real p_inf_gas = MP_CFG.p_inf1;
Real gamma_water = MP_CFG.gamma2;
Real p_inf_water = MP_CFG.p_inf2;
//----------------------------------------------------------------------
//	Initial condition, from config.ini.
//----------------------------------------------------------------------
Real rho_L = MP_CFG.rho_L;
Real p_L = MP_CFG.p_L;
Real u_L = MP_CFG.u_L;
Real alpha_L = MP_CFG.alpha_L;
Real rho_R = MP_CFG.rho_R;
Real p_R = MP_CFG.p_R;
Real u_R = MP_CFG.u_R;
Real alpha_R = MP_CFG.alpha_R;
Real end_time = MP_CFG.end_time;
//----------------------------------------------------------------------
//	Fluid body shape: thin rectangular slab [0, L] x [0, H].
//----------------------------------------------------------------------
class FluidBlock : public ComplexShape
{
  public:
    explicit FluidBlock(const std::string &shape_name) : ComplexShape(shape_name)
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
//----------------------------------------------------------------------
//	Solid wall slab [x0, x1] x [0, H] (left / right end walls).
//----------------------------------------------------------------------
class WallBlock : public MultiPolygonShape
{
  public:
    WallBlock(const std::string &shape_name, Real x0, Real x1) : MultiPolygonShape(shape_name)
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
//	Case-dependent initial condition for the five-equation model.
//	Sets density, pressure, velocity, volume fraction, mass, momentum,
//	and total energy (using the mixture EOS).
//----------------------------------------------------------------------
class MultiPhaseShockTubeInitialCondition : public LocalDynamics
{
  public:
    explicit MultiPhaseShockTubeInitialCondition(SPHBody &sph_body, MultiphaseMixture &mixture)
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
        Real rho_e_int; /**< Internal energy per unit volume. */

        if (x < x_interface)
        {
            rho_[index_i] = rho_L;
            p_[index_i] = p_L;
            vel_[index_i][0] = u_L;
            alpha_[index_i] = alpha_L;
        }
        else
        {
            rho_[index_i] = rho_R;
            p_[index_i] = p_R;
            vel_[index_i][0] = u_R;
            alpha_[index_i] = alpha_R;
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
