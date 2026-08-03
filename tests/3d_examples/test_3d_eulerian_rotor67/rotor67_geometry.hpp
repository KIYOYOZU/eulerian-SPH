/**
 * @file 	rotor67_geometry.hpp
 * @brief 	Geometry, initial condition and diagnostics for the Rotor 67 Phase B
 *          smoke test (fully-compressible Eulerian SPH + rotating frame).
 * @details Single-passage wedge sector about the Z axis: theta in [theta0,
 *          theta0 + delta_theta], r in [r_hub(z), r_shd(z)],
 *          z in [z_in, z_out]. Particles stored in the standard Cartesian
 *          frame. Hub and casing are no longer simplified constant-radius
 *          cylinders — they are tabulated r(z) meridional curves loaded
 *          from CFX-BladeGen .curve files at full scale. This
 *          gives the fluid domain and the wall bodies the same curved hub
 *          and shroud profiles as the reference NASA Rotor 67 geometry.
 *          Circumferential periodicity uses the shared-layer
 *          RotatingPeriodicConditionUsingGhostParticles. The code vel_ field
 *          carries the rotating-frame relative velocity w.
 *
 *          Phase B: adds blade body via TriangleMeshShapeSTL.
 *          Phase C: tip clearance modifies STL path.
 *
 *          Viscosity is handled separately in the main cpp via
 *          fluid_dynamics::ViscousForceInner, not via the geometry layer.
 */
#ifndef ROTOR67_GEOMETRY_H
#define ROTOR67_GEOMETRY_H

#include "rotor67_data.hpp"
#include "sphinxsys.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace SPH
{

//----------------------------------------------------------------------
//  Rotor67MeridionalCurve: tabulated r(z) profiles for the hub and shroud.
//  Loads CFX-BladeGen .curve files (3 whitespace-separated columns:
//  index, r, z; 1200 points each) at FULL SCALE (metres, no scaling),
//  stores sorted (z, r) pairs in ascending-z order, and provides
//  linear-interpolated r(z) queries with z-clamping at curve endpoints.
//  Used by Rotor67DomainScope and Rotor67HubShell/CasingShell for the true
//  meridional profiles.
//
//  Curve file format (whitespace-separated, no header line):
//      <index>\t<r>\t<z>
//      <index>\t<r>\t<z>
//      ...
//  Both r and z are in the CFX units (metres), stored directly.
//----------------------------------------------------------------------
class Rotor67MeridionalCurve
{
  public:
    /// Construct from Rotor67Config. Reads cfg.meridional_hub_curve_path
    /// and cfg.meridional_shroud_curve_path (resolved to absolute paths by
    /// load_rotor67_config) at full scale. Throws std::runtime_error if
    /// either file is missing or malformed.
    explicit Rotor67MeridionalCurve(const Rotor67Config &cfg)
    {
        if (cfg.meridional_hub_curve_path.empty() ||
            cfg.meridional_shroud_curve_path.empty())
        {
            throw std::runtime_error(
                "Rotor67MeridionalCurve: meridional curve paths must be set in [meridional] config");
        }
        hub_ = loadCurve(cfg.meridional_hub_curve_path, "hub");
        shd_ = loadCurve(cfg.meridional_shroud_curve_path, "shroud");
        z_hub_lo_ = hub_.front().first;
        z_hub_hi_ = hub_.back().first;
        z_shd_lo_ = shd_.front().first;
        z_shd_hi_ = shd_.back().first;
        std::cout << "[Rotor67MeridionalCurve] hub: " << hub_.size()
                  << " pts, z in [" << z_hub_lo_ << ", " << z_hub_hi_
                  << "], r in [" << hub_.front().second << ", "
                  << hub_.back().second << "] (full scale, m)\n";
        std::cout << "[Rotor67MeridionalCurve] shroud: " << shd_.size()
                  << " pts, z in [" << z_shd_lo_ << ", " << z_shd_hi_
                  << "], r in [" << shd_.front().second << ", "
                  << shd_.back().second << "]\n";
    };

    /// Linear-interpolated r at given z. Clamps to the curve endpoints
    /// (returns the endpoint r) if z is outside the curve's z range.
    Real hubR(Real z) const { return interpolateR(z, hub_); }
    Real shroudR(Real z) const { return interpolateR(z, shd_); }

    /// dr/dz at given z (linear-regression slope across the bracketing
    /// segment). Returns the endpoint one-sided slope if z is outside
    /// the curve's z range.
    Real hubDrDz(Real z) const { return slopeDrDz(z, hub_); }
    Real shroudDrDz(Real z) const { return slopeDrDz(z, shd_); }

    Real hubZMin() const { return z_hub_lo_; }
    Real hubZMax() const { return z_hub_hi_; }
    Real shroudZMin() const { return z_shd_lo_; }
    Real shroudZMax() const { return z_shd_hi_; }

    /// Largest shroud r over the union of the SPH z range and the
    /// shroud-curve z range. Used to size the system bounding box.
    Real maxShroudR(Real z_lo, Real z_hi) const
    {
        Real r_max = Real(0.0);
        const Real z_scan_lo = std::max(z_lo, z_shd_lo_);
        const Real z_scan_hi = std::min(z_hi, z_shd_hi_);
        if (z_scan_lo < z_scan_hi)
        {
            for (const auto &zp : shd_)
            {
                if (zp.first >= z_scan_lo && zp.first <= z_scan_hi)
                {
                    r_max = std::max(r_max, zp.second);
                }
            }
        }
        // Account for clamp extension at z-scan endpoints.
        r_max = std::max(r_max, shroudR(std::max(z_lo, z_shd_lo_)));
        r_max = std::max(r_max, shroudR(std::min(z_hi, z_shd_hi_)));
        return r_max;
    }

  private:
    using CurvePoints = std::vector<std::pair<Real, Real>>; // (z, r) ascending

    static CurvePoints loadCurve(const std::string &path,
                                 const std::string &label)
    {
        std::ifstream in(path);
        if (!in.good())
        {
            throw std::runtime_error(
                "Rotor67MeridionalCurve: cannot open " + label + " curve at '" +
                path + "'. Set [meridional] in config.ini or check the file.");
        }
        CurvePoints pts;
        pts.reserve(2048);
        std::string line;
        size_t line_no = 0;
        while (std::getline(in, line))
        {
            ++line_no;
            // Skip blank lines and any line that starts with a non-digit
            // (e.g. '#' comments if any are added later).
            std::istringstream iss(line);
            Real col0, col1, col2;
            if (!(iss >> col0 >> col1 >> col2))
            {
                continue;
            }
            pts.emplace_back(col2, col1); // (z, r) full scale
        }
        if (pts.size() < 2)
        {
            throw std::runtime_error(
                "Rotor67MeridionalCurve: " + label + " curve '" + path +
                "' has fewer than 2 valid points (got " +
                std::to_string(pts.size()) + ")");
        }
        // Sort ascending by z (in case the source is not monotonic).
        std::sort(pts.begin(), pts.end(),
                  [](const std::pair<Real, Real> &a,
                     const std::pair<Real, Real> &b)
                  { return a.first < b.first; });
        return pts;
    }

    static Real interpolateR(Real z, const CurvePoints &pts)
    {
        // Clamp to endpoints.
        if (z <= pts.front().first)
        {
            return pts.front().second;
        }
        if (z >= pts.back().first)
        {
            return pts.back().second;
        }
        // Binary search for the upper bracket z_{i+1} with z_i <= z < z_{i+1}.
        auto lo = pts.begin();
        auto hi = pts.end() - 1;
        while (hi - lo > 1)
        {
            auto mid = lo + (hi - lo) / 2;
            if (mid->first <= z)
            {
                lo = mid;
            }
            else
            {
                hi = mid;
            }
        }
        const Real z0 = lo->first;
        const Real z1 = hi->first;
        const Real r0 = lo->second;
        const Real r1 = hi->second;
        const Real t = (z - z0) / (z1 - z0);
        return r0 + t * (r1 - r0);
    }

    static Real slopeDrDz(Real z, const CurvePoints &pts)
    {
        if (z <= pts.front().first)
        {
            return (pts[1].second - pts[0].second) /
                   (pts[1].first - pts[0].first);
        }
        if (z >= pts.back().first)
        {
            const size_t n = pts.size();
            return (pts[n - 1].second - pts[n - 2].second) /
                   (pts[n - 1].first - pts[n - 2].first);
        }
        auto lo = pts.begin();
        auto hi = pts.end() - 1;
        while (hi - lo > 1)
        {
            auto mid = lo + (hi - lo) / 2;
            if (mid->first <= z)
            {
                lo = mid;
            }
            else
            {
                hi = mid;
            }
        }
        return (hi->second - lo->second) / (hi->first - lo->first);
    }

    CurvePoints hub_;
    CurvePoints shd_;
    Real z_hub_lo_, z_hub_hi_;
    Real z_shd_lo_, z_shd_hi_;
};

inline Real rotor67NormalisedAngleDiff(Real a, Real b)
{
    Real d = a - b;
    while (d > Pi)
        d -= Real(2.0) * Pi;
    while (d < -Pi)
        d += Real(2.0) * Pi;
    return d;
}

struct Rotor67InletPrimitiveState
{
    Real p_total;
    Real T_total;
    Real p_static;
    Real T_static;
    Real rho_static;
    Vecd velocity_absolute;
    Vecd velocity_relative;
};

inline Real rotor67SpanFraction(const Vecd &pos, const Rotor67Config &cfg,
                                const Rotor67MeridionalCurve &curve)
{
    const Real r = std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
    const Real r_hub = curve.hubR(pos[2]);
    const Real r_shroud = curve.shroudR(pos[2]);
    const Real denom = std::max(r_shroud - r_hub, TinyReal);
    const Real s = (r - r_hub) / denom;
    return std::max(Real(0.0), std::min(Real(1.0), s));
}

inline Real rotor67InletTotalPressure(const Vecd &pos, const Rotor67Config &cfg,
                                      const Rotor67MeridionalCurve &curve)
{
    const Real s = rotor67SpanFraction(pos, cfg, curve);
    if (cfg.inlet_total_pressure_profile == "tip_radial")
    {
        return cfg.inlet_total_pressure * (Real(1.0) - cfg.inlet_distortion_intensity * s);
    }
    if (cfg.inlet_total_pressure_profile == "hub_radial")
    {
        return cfg.inlet_total_pressure *
               (Real(1.0) - cfg.inlet_distortion_intensity * (Real(1.0) - s));
    }
    return cfg.inlet_total_pressure;
}

inline Vecd rotor67InletAbsoluteVelocity(const Vecd &pos, const Rotor67Config &cfg)
{
    Vecd velocity = Vecd(Real(0.0), Real(0.0), cfg.inlet_axial_velocity);
    if (cfg.inlet_swirl_type == "none" ||
        std::fabs(cfg.inlet_swirl_angle_deg) <= Real(1.0e-12))
    {
        return velocity;
    }
    const Real r = std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
    if (r <= TinyReal)
    {
        return velocity;
    }
    const Real sign = cfg.inlet_swirl_type == "counter" ? Real(-1.0) : Real(1.0);
    const Real angle = cfg.inlet_swirl_angle_deg * Pi / Real(180.0);
    const Real u_theta = sign * cfg.inlet_axial_velocity * std::tan(angle);
    const Vecd theta_dir(-pos[1] / r, pos[0] / r, Real(0.0));
    return velocity + u_theta * theta_dir;
}

inline Vecd rotor67InletRelativeVelocity(const Vecd &pos, const Rotor67Config &cfg)
{
    const Vecd frame_velocity(-cfg.omega * pos[1], cfg.omega * pos[0], Real(0.0));
    return rotor67InletAbsoluteVelocity(pos, cfg) - frame_velocity;
}

inline Rotor67InletPrimitiveState rotor67InletStaticState(
    const Vecd &pos, const Rotor67Config &cfg, const Rotor67MeridionalCurve &curve)
{
    Rotor67InletPrimitiveState state{};
    state.p_total = rotor67InletTotalPressure(pos, cfg, curve);
    state.T_total = cfg.inlet_total_temperature;
    state.velocity_absolute = rotor67InletAbsoluteVelocity(pos, cfg);
    state.velocity_relative = rotor67InletRelativeVelocity(pos, cfg);

    // Rotor67 inlet total quantities are defined in the absolute-frame AIP.
    // Use |u_abs|, not |w|, for the isentropic total-to-static conversion.
    const Real cp = cfg.gamma * cfg.gas_constant / (cfg.gamma - Real(1.0));
    state.T_static =
        state.T_total -
        state.velocity_absolute.squaredNorm() / (Real(2.0) * cp);
    if (state.T_static <= TinyReal)
    {
        throw std::runtime_error(
            "Rotor67 inlet total state is invalid: local absolute velocity exceeds total-temperature enthalpy");
    }
    state.p_static =
        state.p_total *
        std::pow(state.T_static / state.T_total, cfg.gamma / (cfg.gamma - Real(1.0)));
    state.rho_static = state.p_static / (cfg.gas_constant * state.T_static);
    return state;
}

enum class Rotor67DomainSurface
{
    Hub = 0,
    Shroud,
    Inlet,
    Outlet,
    PeriodicMin,
    PeriodicMax,
    Blade,
    Count
};

inline constexpr size_t rotor67DomainSurfaceCount()
{
    return static_cast<size_t>(Rotor67DomainSurface::Count);
}

inline const char *rotor67DomainSurfaceName(Rotor67DomainSurface surface)
{
    switch (surface)
    {
    case Rotor67DomainSurface::Hub:
        return "rotor67_domain_hub";
    case Rotor67DomainSurface::Shroud:
        return "rotor67_domain_shroud";
    case Rotor67DomainSurface::Inlet:
        return "rotor67_domain_inlet";
    case Rotor67DomainSurface::Outlet:
        return "rotor67_domain_outlet";
    case Rotor67DomainSurface::PeriodicMin:
        return "rotor67_domain_periodic_min";
    case Rotor67DomainSurface::PeriodicMax:
        return "rotor67_domain_periodic_max";
    case Rotor67DomainSurface::Blade:
        return "rotor67_blade_profile";
    case Rotor67DomainSurface::Count:
        break;
    }
    return "unknown";
}

inline Rotor67DomainSurface rotor67DomainSurfaceFromName(const std::string &name)
{
    for (size_t i = 0; i != rotor67DomainSurfaceCount(); ++i)
    {
        const auto surface = static_cast<Rotor67DomainSurface>(i);
        if (name == rotor67DomainSurfaceName(surface))
        {
            return surface;
        }
    }
    throw std::runtime_error("Unexpected Rotor67 domain STL solid name: " + name);
}

struct Rotor67DomainSurfaceMesh
{
    std::vector<std::array<Real, 3>> vertices_;
    std::vector<std::array<int, 3>> faces_;
    tmd::TriangleMeshDistance distance_;
    BoundingBoxd bounds_;
    bool has_bounds_ = false;
    bool distance_ready_ = false;
    std::map<std::string, int> vertex_index_;

    int addVertex(const std::array<Real, 3> &v)
    {
        if (!std::isfinite(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2]))
        {
            throw std::runtime_error("Rotor67 domain STL contains non-finite vertex coordinate");
        }
        std::ostringstream key;
        key << std::setprecision(17) << v[0] << "," << v[1] << "," << v[2];
        const std::string key_string = key.str();
        auto it = vertex_index_.find(key_string);
        if (it != vertex_index_.end())
        {
            return it->second;
        }
        const int index = static_cast<int>(vertices_.size());
        vertices_.push_back(v);
        vertex_index_[key_string] = index;
        const Vecd p(v[0], v[1], v[2]);
        if (!has_bounds_)
        {
            bounds_ = BoundingBoxd(p, p);
            has_bounds_ = true;
        }
        else
        {
            for (int d = 0; d != 3; ++d)
            {
                bounds_.lower_[d] = std::min(bounds_.lower_[d], p[d]);
                bounds_.upper_[d] = std::max(bounds_.upper_[d], p[d]);
            }
        }
        return index;
    }

    void addFace(const std::array<std::array<Real, 3>, 3> &tri)
    {
        faces_.push_back({addVertex(tri[0]), addVertex(tri[1]), addVertex(tri[2])});
    }
};

class Rotor67DomainSTL
{
  public:
    explicit Rotor67DomainSTL(const Rotor67Config &cfg)
        : path_(cfg.domain_stl_path)
    {
        load();
        validateRequiredSurfaces();
        constructExteriorDomainMesh();
        constructDistances();
    }

    const Rotor67DomainSurfaceMesh &surface(Rotor67DomainSurface surface) const
    {
        return surfaces_[static_cast<size_t>(surface)];
    }

    size_t faceCount(Rotor67DomainSurface surface) const
    {
        return this->surface(surface).faces_.size();
    }

    BoundingBoxd bounds(Real padding = Real(0.0)) const
    {
        BoundingBoxd box = global_bounds_;
        box.lower_ -= Vecd(padding, padding, padding);
        box.upper_ += Vecd(padding, padding, padding);
        return box;
    }

    Real periodicThetaMin() const { return periodic_theta_min_; }
    Real periodicThetaMax() const { return periodic_theta_max_; }
    Real periodicDeltaTheta() const { return periodic_delta_theta_; }
    Real inletZ() const { return surface(Rotor67DomainSurface::Inlet).bounds_.lower_[2]; }
    Real outletZ() const { return surface(Rotor67DomainSurface::Outlet).bounds_.upper_[2]; }

    Real unsignedDistance(Rotor67DomainSurface surface, const Vecd &pnt) const
    {
        auto result = this->surface(surface).distance_.unsigned_distance(pnt);
        return static_cast<Real>(result.distance);
    }

    Real bladeSignedDistance(const Vecd &pnt) const
    {
        auto result = surface(Rotor67DomainSurface::Blade).distance_.signed_distance(pnt);
        return static_cast<Real>(result.distance);
    }

    Real bladeUnsignedDistance(const Vecd &pnt) const
    {
        return unsignedDistance(Rotor67DomainSurface::Blade, pnt);
    }

    bool bladeContainsPoint(const Vecd &pnt) const
    {
        const auto &mesh = surface(Rotor67DomainSurface::Blade);
        if (!pointInBounds(mesh.bounds_, pnt, Eps))
        {
            return false;
        }
        const Real signed_distance = bladeSignedDistance(pnt);
        if (signed_distance > Eps)
        {
            return false;
        }
        if (std::abs(signed_distance) <= Eps)
        {
            return true;
        }
        return pointInsideClosedMeshByRay(mesh, pnt);
    }

    bool domainContainsPoint(const Vecd &pnt) const
    {
        if (!pointInBounds(exterior_bounds_, pnt, Eps))
        {
            return false;
        }
        return pointInsideClosedMeshByRay(exterior_domain_mesh_, pnt);
    }

    Vecd closestPoint(Rotor67DomainSurface surface, const Vecd &pnt) const
    {
        auto result = this->surface(surface).distance_.unsigned_distance(pnt);
        return Vecd(result.nearest_point[0], result.nearest_point[1], result.nearest_point[2]);
    }

    Vecd closestDomainPoint(const Vecd &pnt) const
    {
        Real best = MaxReal;
        Vecd best_point = pnt;
        for (size_t i = 0; i != rotor67DomainSurfaceCount(); ++i)
        {
            const auto surface_id = static_cast<Rotor67DomainSurface>(i);
            auto result = surface(surface_id).distance_.unsigned_distance(pnt);
            if (result.distance < best)
            {
                best = static_cast<Real>(result.distance);
                best_point = Vecd(result.nearest_point[0], result.nearest_point[1], result.nearest_point[2]);
            }
        }
        return best_point;
    }

    void printSummary() const
    {
        std::cout << "[RotorDiag][DomainSTL] path=" << path_ << "\n";
        for (size_t i = 0; i != rotor67DomainSurfaceCount(); ++i)
        {
            const auto surface_id = static_cast<Rotor67DomainSurface>(i);
            const auto &mesh = surface(surface_id);
            std::cout << "[RotorDiag][DomainSTL] solid=" << rotor67DomainSurfaceName(surface_id)
                      << " faces=" << mesh.faces_.size()
                      << " bounds=[" << mesh.bounds_.lower_.transpose()
                      << "]..[" << mesh.bounds_.upper_.transpose() << "]\n";
        }
    }

  private:
    void load()
    {
        std::ifstream in(path_);
        if (!in.is_open())
        {
            throw std::runtime_error("Cannot open Rotor67 domain STL: " + path_);
        }
        bool in_solid = false;
        Rotor67DomainSurface current = Rotor67DomainSurface::Count;
        std::array<std::array<Real, 3>, 3> tri{};
        int tri_vertices = 0;
        std::string line;
        while (std::getline(in, line))
        {
            std::istringstream iss(line);
            std::string token;
            iss >> token;
            if (token.empty())
            {
                continue;
            }
            if (token == "solid")
            {
                if (in_solid)
                {
                    throw std::runtime_error("Malformed Rotor67 domain STL: nested solid block");
                }
                std::string name;
                iss >> name;
                current = rotor67DomainSurfaceFromName(name);
                auto &block_count = solid_block_counts_[static_cast<size_t>(current)];
                block_count++;
                if (block_count > 1)
                {
                    throw std::runtime_error("Duplicate Rotor67 domain STL solid: " + name);
                }
                in_solid = true;
                tri_vertices = 0;
                continue;
            }
            if (token == "endsolid")
            {
                if (!in_solid || tri_vertices != 0)
                {
                    throw std::runtime_error("Malformed Rotor67 domain STL: invalid endsolid state");
                }
                in_solid = false;
                current = Rotor67DomainSurface::Count;
                tri_vertices = 0;
                continue;
            }
            if (token == "vertex" && in_solid)
            {
                if (tri_vertices >= 3)
                {
                    throw std::runtime_error("Malformed Rotor67 domain STL: more than 3 vertices in facet");
                }
                Real x, y, z;
                if (!(iss >> x >> y >> z))
                {
                    throw std::runtime_error("Malformed Rotor67 domain STL vertex line: " + line);
                }
                tri[tri_vertices++] = {x, y, z};
                continue;
            }
            if (token == "endfacet" && in_solid)
            {
                if (tri_vertices != 3)
                {
                    throw std::runtime_error("Malformed Rotor67 domain STL: facet without exactly 3 vertices");
                }
                surfaces_[static_cast<size_t>(current)].addFace(tri);
                tri_vertices = 0;
            }
        }
        if (in_solid || tri_vertices != 0)
        {
            throw std::runtime_error("Malformed Rotor67 domain STL: unterminated solid or facet");
        }
    }

    void validateRequiredSurfaces()
    {
        bool have_global = false;
        for (size_t i = 0; i != rotor67DomainSurfaceCount(); ++i)
        {
            const auto surface_id = static_cast<Rotor67DomainSurface>(i);
            auto &mesh = surfaces_[i];
            if (mesh.faces_.empty())
            {
                throw std::runtime_error(
                    "Rotor67 domain STL missing required non-empty solid: " +
                    std::string(rotor67DomainSurfaceName(surface_id)));
            }
            if (!have_global)
            {
                global_bounds_ = mesh.bounds_;
                have_global = true;
            }
            else
            {
                for (int d = 0; d != 3; ++d)
                {
                    global_bounds_.lower_[d] = std::min(global_bounds_.lower_[d], mesh.bounds_.lower_[d]);
                    global_bounds_.upper_[d] = std::max(global_bounds_.upper_[d], mesh.bounds_.upper_[d]);
                }
            }
        }
        periodic_theta_min_ = surfaceMeanTheta(surface(Rotor67DomainSurface::PeriodicMin));
        periodic_theta_max_ = surfaceMeanTheta(surface(Rotor67DomainSurface::PeriodicMax));
        periodic_delta_theta_ = positiveAngleDiff(periodic_theta_max_, periodic_theta_min_);
    }

    void constructDistances()
    {
        for (auto &mesh : surfaces_)
        {
            mesh.distance_.construct(mesh.vertices_, mesh.faces_);
            mesh.distance_ready_ = true;
        }
        if (!surface(Rotor67DomainSurface::Blade).distance_.is_mesh_manifold())
        {
            throw std::runtime_error(
                "Rotor67 domain STL blade solid is not manifold; signed distance is unsafe");
        }
        exterior_domain_mesh_.distance_.construct(exterior_domain_mesh_.vertices_,
                                                  exterior_domain_mesh_.faces_);
        exterior_domain_mesh_.distance_ready_ = true;
        if (!exterior_domain_mesh_.distance_.is_mesh_manifold())
        {
            throw std::runtime_error(
                "Rotor67 exterior domain STL surfaces do not form a closed manifold volume");
        }
        const Vecd inside_probe = findExteriorInsideProbe();
        auto result = exterior_domain_mesh_.distance_.signed_distance(inside_probe);
        exterior_inside_sign_ = result.distance >= 0.0 ? Real(1.0) : Real(-1.0);
    }

    void constructExteriorDomainMesh()
    {
        bool have_bounds = false;
        for (size_t i = 0; i != rotor67DomainSurfaceCount(); ++i)
        {
            const auto surface_id = static_cast<Rotor67DomainSurface>(i);
            if (surface_id == Rotor67DomainSurface::Blade)
            {
                continue;
            }
            const auto &mesh = surface(surface_id);
            for (const auto &face : mesh.faces_)
            {
                const auto &v0 = mesh.vertices_[static_cast<size_t>(face[0])];
                const auto &v1 = mesh.vertices_[static_cast<size_t>(face[1])];
                const auto &v2 = mesh.vertices_[static_cast<size_t>(face[2])];
                exterior_domain_mesh_.addFace({v0, v1, v2});
            }
            if (!have_bounds)
            {
                exterior_bounds_ = mesh.bounds_;
                have_bounds = true;
            }
            else
            {
                for (int d = 0; d != 3; ++d)
                {
                    exterior_bounds_.lower_[d] = std::min(exterior_bounds_.lower_[d], mesh.bounds_.lower_[d]);
                    exterior_bounds_.upper_[d] = std::max(exterior_bounds_.upper_[d], mesh.bounds_.upper_[d]);
                }
            }
        }
    }

    Vecd findExteriorInsideProbe() const
    {
        const Vecd lower = exterior_bounds_.lower_;
        const Vecd upper = exterior_bounds_.upper_;
        const Vecd center = Real(0.5) * (lower + upper);
        if (pointInsideClosedMeshByRay(exterior_domain_mesh_, center))
        {
            return center;
        }
        for (int ix = 1; ix != 8; ++ix)
        {
            for (int iy = 1; iy != 8; ++iy)
            {
                for (int iz = 1; iz != 8; ++iz)
                {
                    const Vecd pnt(lower[0] + (upper[0] - lower[0]) * Real(ix) / Real(8.0),
                                   lower[1] + (upper[1] - lower[1]) * Real(iy) / Real(8.0),
                                   lower[2] + (upper[2] - lower[2]) * Real(iz) / Real(8.0));
                    if (pointInsideClosedMeshByRay(exterior_domain_mesh_, pnt))
                    {
                        return pnt;
                    }
                }
            }
        }
        throw std::runtime_error("Cannot find an interior probe in Rotor67 exterior domain STL");
    }

    static bool pointInBounds(const BoundingBoxd &bounds, const Vecd &pnt, Real tol)
    {
        for (int d = 0; d != 3; ++d)
        {
            if (pnt[d] < bounds.lower_[d] - tol || pnt[d] > bounds.upper_[d] + tol)
            {
                return false;
            }
        }
        return true;
    }

    static Real surfaceMeanTheta(const Rotor67DomainSurfaceMesh &mesh)
    {
        Real c = 0.0;
        Real s = 0.0;
        for (const auto &v : mesh.vertices_)
        {
            const Real theta = std::atan2(v[1], v[0]);
            c += std::cos(theta);
            s += std::sin(theta);
        }
        return std::atan2(s, c);
    }

    static Real positiveAngleDiff(Real a, Real b)
    {
        Real d = a - b;
        while (d < Real(0.0))
            d += Real(2.0) * Pi;
        while (d > Real(2.0) * Pi)
            d -= Real(2.0) * Pi;
        return d;
    }

    static bool pointInsideClosedMeshByRay(const Rotor67DomainSurfaceMesh &mesh, const Vecd &pnt)
    {
        const Vecd dir(1.0, 0.0, 0.0);
        constexpr Real eps = 1.0e-12;
        std::vector<Real> hits;
        hits.reserve(16);
        for (const auto &face : mesh.faces_)
        {
            const auto &a0 = mesh.vertices_[static_cast<size_t>(face[0])];
            const auto &a1 = mesh.vertices_[static_cast<size_t>(face[1])];
            const auto &a2 = mesh.vertices_[static_cast<size_t>(face[2])];
            const Vecd v0(a0[0], a0[1], a0[2]);
            const Vecd v1(a1[0], a1[1], a1[2]);
            const Vecd v2(a2[0], a2[1], a2[2]);
            if (pnt[1] < std::min({v0[1], v1[1], v2[1]}) - eps ||
                pnt[1] > std::max({v0[1], v1[1], v2[1]}) + eps ||
                pnt[2] < std::min({v0[2], v1[2], v2[2]}) - eps ||
                pnt[2] > std::max({v0[2], v1[2], v2[2]}) + eps)
            {
                continue;
            }

            const Vecd edge1 = v1 - v0;
            const Vecd edge2 = v2 - v0;
            const Vecd h = dir.cross(edge2);
            const Real det = edge1.dot(h);
            if (std::abs(det) < eps)
            {
                continue;
            }
            const Real inv_det = Real(1.0) / det;
            const Vecd s = pnt - v0;
            const Real u = inv_det * s.dot(h);
            if (u < -eps || u > Real(1.0) + eps)
            {
                continue;
            }
            const Vecd q = s.cross(edge1);
            const Real v = inv_det * dir.dot(q);
            if (v < -eps || u + v > Real(1.0) + eps)
            {
                continue;
            }
            const Real t = inv_det * edge2.dot(q);
            if (t > eps)
            {
                hits.push_back(t);
            }
        }

        std::sort(hits.begin(), hits.end());
        size_t unique_hits = 0;
        Real last = -MaxReal;
        for (Real hit : hits)
        {
            if (unique_hits == 0 || std::abs(hit - last) > Real(1.0e-9))
            {
                unique_hits++;
                last = hit;
            }
        }
        return (unique_hits % 2) == 1;
    }

    std::string path_;
    std::array<Rotor67DomainSurfaceMesh, rotor67DomainSurfaceCount()> surfaces_;
    std::array<int, rotor67DomainSurfaceCount()> solid_block_counts_{};
    Rotor67DomainSurfaceMesh exterior_domain_mesh_;
    BoundingBoxd exterior_bounds_;
    Real exterior_inside_sign_ = Real(1.0);
    BoundingBoxd global_bounds_;
    Real periodic_theta_min_ = Real(0.0);
    Real periodic_theta_max_ = Real(0.0);
    Real periodic_delta_theta_ = Real(0.0);
};

class Rotor67DomainScope
{
  public:
    Rotor67DomainScope(const Rotor67Config &cfg, const Rotor67MeridionalCurve &curve,
                       const Rotor67DomainSTL &domain_stl)
        : cfg_(cfg), curve_(curve), domain_stl_(domain_stl),
          wall_offset_(cfg.particle_boundary_offset * cfg.global_resolution),
          theta0_(domain_stl.periodicThetaMin()),
          delta_theta_(domain_stl.periodicDeltaTheta()),
          z_in_(domain_stl.inletZ()),
          z_out_(domain_stl.outletZ()),
          reservoir_width_(cfg.sponge_width_factor * cfg.global_resolution),
          z_comp_min_(z_in_ - reservoir_width_),
          z_comp_max_(z_out_ + reservoir_width_)
    {
    }

    bool inSector(const Vecd &pnt, bool boundary_included = true) const
    {
        const Real r = std::sqrt(pnt[0] * pnt[0] + pnt[1] * pnt[1]);
        if (r < TinyReal)
        {
            return false;
        }
        const Real dtheta = rotor67NormalisedAngleDiff(std::atan2(pnt[1], pnt[0]), theta0_);
        return boundary_included ? (dtheta >= -Eps && dtheta <= delta_theta_ + Eps)
                                 : (dtheta > Eps && dtheta < delta_theta_ - Eps);
    }

    bool inAxialRange(const Vecd &pnt, bool boundary_included = true) const
    {
        return boundary_included ? (pnt[2] >= z_comp_min_ - Eps && pnt[2] <= z_comp_max_ + Eps)
                                 : (pnt[2] > z_comp_min_ + Eps && pnt[2] < z_comp_max_ - Eps);
    }

    bool inFluidRadialRange(const Vecd &pnt, bool boundary_included = true) const
    {
        const Real r = std::sqrt(pnt[0] * pnt[0] + pnt[1] * pnt[1]);
        const Real r_lo = curve_.hubR(pnt[2]) + wall_offset_;
        const Real r_hi = curve_.shroudR(pnt[2]) - wall_offset_;
        return boundary_included ? (r >= r_lo - Eps && r <= r_hi + Eps)
                                 : (r > r_lo + Eps && r < r_hi - Eps);
    }

    bool inBladeWallEnvelope(const Vecd &pnt) const
    {
        const Real r = std::sqrt(pnt[0] * pnt[0] + pnt[1] * pnt[1]);
        const Real half = Real(0.5) * static_cast<Real>(cfg_.wall_layers) * cfg_.global_resolution;
        return r >= curve_.hubR(pnt[2]) - half - Eps &&
               r <= curve_.shroudR(pnt[2]) + half + Eps;
    }

    bool inFluidDuct(const Vecd &pnt, bool boundary_included = true) const
    {
        return inSector(pnt, boundary_included) &&
               inAxialRange(pnt, boundary_included) &&
               inFluidRadialRange(pnt, boundary_included);
    }

    BoundingBoxd bounds(Real padding = Real(0.0)) const
    {
        const Real r_max =
            std::max(domain_stl_.bounds().upper_[0],
                     curve_.maxShroudR(z_comp_min_, z_comp_max_) + wall_offset_);
        BoundingBoxd box(Vecd(-r_max - padding, -r_max - padding, z_comp_min_ - padding),
                         Vecd(r_max + padding, r_max + padding, z_comp_max_ + padding));
        return box;
    }

    Real physicalInletZ() const { return z_in_; }
    Real physicalOutletZ() const { return z_out_; }
    Real computationalZMin() const { return z_comp_min_; }
    Real computationalZMax() const { return z_comp_max_; }
    Real reservoirWidth() const { return reservoir_width_; }

    const Rotor67Config &cfg_;
    const Rotor67MeridionalCurve &curve_;
    const Rotor67DomainSTL &domain_stl_;
    Real wall_offset_;
    Real theta0_;
    Real delta_theta_;
    Real z_in_;
    Real z_out_;
    Real reservoir_width_;
    Real z_comp_min_;
    Real z_comp_max_;
};

class Rotor67DomainFluidVolume : public Shape
{
  public:
    Rotor67DomainFluidVolume(const std::string &shape_name, const Rotor67Config &cfg,
                             const Rotor67MeridionalCurve &curve,
                             const Rotor67DomainSTL &domain_stl)
        : Shape(shape_name), domain_stl_(domain_stl),
          domain_scope_(cfg, curve, domain_stl)
    {
        bounding_box_ = domain_scope_.bounds();
    }

    bool checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED = true) override
    {
        if (!domain_scope_.inFluidDuct(pnt, BOUNDARY_INCLUDED))
        {
            return false;
        }
        return !domain_stl_.bladeContainsPoint(pnt);
    }

    Vecd findClosestPoint(const Vecd &probe_point) override
    {
        return domain_stl_.closestDomainPoint(probe_point);
    }

    BoundingBoxd findBounds() override { return bounding_box_; }

  private:
    const Rotor67DomainSTL &domain_stl_;
    Rotor67DomainScope domain_scope_;
};

struct Rotor67GeometryAudit
{
    std::array<size_t, rotor67DomainSurfaceCount()> named_surface_counts{};
    Real periodic_delta_theta = 0.0;
    Real inlet_hub_tip_ratio = 0.0;
    Real exit_hub_tip_ratio = 0.0;
    Real min_tip_clearance_m = 0.0;
    Real tip_clearance_over_dp = 0.0;
    Real reservoir_width_m = 0.0;
    Real computational_z_min = 0.0;
    Real computational_z_max = 0.0;
    std::string named_surface_status;
    std::string pitch_status;
    std::string tip_clearance_status;
    std::string literature_geometry_status;
    std::string method;
};

inline Rotor67GeometryAudit buildRotor67GeometryAudit(
    const Rotor67Config &cfg, const Rotor67MeridionalCurve &curve,
    const Rotor67DomainSTL &domain_stl)
{
    Rotor67GeometryAudit audit;
    bool named_ok = true;
    for (size_t i = 0; i != rotor67DomainSurfaceCount(); ++i)
    {
        const auto surface = static_cast<Rotor67DomainSurface>(i);
        audit.named_surface_counts[i] = domain_stl.faceCount(surface);
        named_ok = named_ok && audit.named_surface_counts[i] > 0;
    }

    audit.periodic_delta_theta = domain_stl.periodicDeltaTheta();
    const Real expected_pitch = Real(2.0) * Pi / static_cast<Real>(cfg.n_blades);
    audit.pitch_status =
        std::fabs(audit.periodic_delta_theta - expected_pitch) <= Real(1.0e-8)
            ? "pass"
            : "fail";
    audit.named_surface_status = named_ok ? "pass" : "fail";

    const Real z_in = domain_stl.inletZ();
    const Real z_out = domain_stl.outletZ();
    audit.inlet_hub_tip_ratio = curve.hubR(z_in) / std::max(curve.shroudR(z_in), TinyReal);
    audit.exit_hub_tip_ratio = curve.hubR(z_out) / std::max(curve.shroudR(z_out), TinyReal);
    audit.reservoir_width_m = cfg.sponge_width_factor * cfg.global_resolution;
    audit.computational_z_min = z_in - audit.reservoir_width_m;
    audit.computational_z_max = z_out + audit.reservoir_width_m;

    const auto &blade = domain_stl.surface(Rotor67DomainSurface::Blade);
    Real min_clearance = MaxReal;
    size_t candidate_count = 0;
    for (const auto &v : blade.vertices_)
    {
        const Vecd p(v[0], v[1], v[2]);
        const Real r = std::sqrt(p[0] * p[0] + p[1] * p[1]);
        const Real r_hub = curve.hubR(p[2]);
        const Real r_shroud = curve.shroudR(p[2]);
        const Real span = (r - r_hub) / std::max(r_shroud - r_hub, TinyReal);
        if (span >= Real(0.80) && span <= Real(1.05))
        {
            min_clearance = std::min(min_clearance, std::fabs(r_shroud - r));
            ++candidate_count;
        }
    }
    audit.min_tip_clearance_m = candidate_count > 0 ? min_clearance : Real(-1.0);
    audit.tip_clearance_over_dp =
        audit.min_tip_clearance_m > Real(0.0)
            ? audit.min_tip_clearance_m / cfg.global_resolution
            : Real(-1.0);
    audit.tip_clearance_status =
        audit.tip_clearance_over_dp >= Real(5.0)
            ? "literature-resolution-pass"
            : (audit.tip_clearance_over_dp >= Real(3.0) ? "trend-resolution-pass" : "smoke-only");
    audit.literature_geometry_status =
        (audit.named_surface_status == "pass" &&
         audit.pitch_status == "pass" &&
         audit.tip_clearance_over_dp >= Real(5.0))
            ? "literature-geometry-ready"
            : (audit.tip_clearance_over_dp >= Real(3.0) ? "trend-geometry-only" : "blocked");
    audit.method =
        "blade vertices with span>=0.80 against meridional shroud radius; coarse diagnostic, not CAD metrology";
    return audit;
}

inline void writeRotor67GeometryAuditJson(
    const std::string &path, const Rotor67GeometryAudit &audit)
{
    std::ofstream out(path);
    if (!out.is_open())
    {
        std::cerr << "[RotorDiag][GeometryAudit] failed to open " << path << "\n";
        return;
    }
    out << std::setprecision(12);
    out << "{\n";
    out << "  \"named_surface_counts\": {\n";
    for (size_t i = 0; i != rotor67DomainSurfaceCount(); ++i)
    {
        const auto surface = static_cast<Rotor67DomainSurface>(i);
        out << "    \"" << rotor67DomainSurfaceName(surface) << "\": "
            << audit.named_surface_counts[i]
            << (i + 1 == rotor67DomainSurfaceCount() ? "\n" : ",\n");
    }
    out << "  },\n";
    out << "  \"periodic_delta_theta\": " << audit.periodic_delta_theta << ",\n";
    out << "  \"inlet_hub_tip_ratio\": " << audit.inlet_hub_tip_ratio << ",\n";
    out << "  \"exit_hub_tip_ratio\": " << audit.exit_hub_tip_ratio << ",\n";
    out << "  \"min_tip_clearance_m\": " << audit.min_tip_clearance_m << ",\n";
    out << "  \"literature_tip_clearance_m\": " << rotor67_literature_tip_clearance << ",\n";
    out << "  \"tip_clearance_over_dp\": " << audit.tip_clearance_over_dp << ",\n";
    out << "  \"reservoir_width_m\": " << audit.reservoir_width_m << ",\n";
    out << "  \"computational_z_min\": " << audit.computational_z_min << ",\n";
    out << "  \"computational_z_max\": " << audit.computational_z_max << ",\n";
    out << "  \"named_surface_status\": \"" << audit.named_surface_status << "\",\n";
    out << "  \"pitch_status\": \"" << audit.pitch_status << "\",\n";
    out << "  \"tip_clearance_status\": \"" << audit.tip_clearance_status << "\",\n";
    out << "  \"literature_geometry_status\": \"" << audit.literature_geometry_status << "\",\n";
    out << "  \"method\": \"" << audit.method << "\"\n";
    out << "}\n";
}

class Rotor67BladeWallVolume : public Shape
{
  public:
    Rotor67BladeWallVolume(const std::string &shape_name, const Rotor67Config &cfg,
                           const Rotor67MeridionalCurve &,
                           const Rotor67DomainSTL &domain_stl)
        : Shape(shape_name), domain_stl_(domain_stl),
          boundary_tolerance_(cfg.particle_boundary_offset * cfg.global_resolution)
    {
        const BoundingBoxd blade_bounds = domain_stl_.surface(Rotor67DomainSurface::Blade).bounds_;
        bounding_box_ = BoundingBoxd(
            blade_bounds.lower_ - Vecd(boundary_tolerance_, boundary_tolerance_, boundary_tolerance_),
            blade_bounds.upper_ + Vecd(boundary_tolerance_, boundary_tolerance_, boundary_tolerance_));
    }

    bool checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED = true) override
    {
        return domain_stl_.bladeContainsPoint(pnt);
    }

    Vecd findClosestPoint(const Vecd &probe_point) override
    {
        return domain_stl_.closestPoint(Rotor67DomainSurface::Blade, probe_point);
    }

    BoundingBoxd findBounds() override { return bounding_box_; }

  private:
    const Rotor67DomainSTL &domain_stl_;
    Real boundary_tolerance_;
};

class Rotor67BladeContactWallVolume : public Shape
{
  public:
    Rotor67BladeContactWallVolume(const std::string &shape_name, const Rotor67Config &cfg,
                                  const Rotor67MeridionalCurve &,
                                  const Rotor67DomainSTL &domain_stl)
        : Shape(shape_name), domain_stl_(domain_stl),
          first_layer_offset_(cfg.particle_boundary_offset * cfg.global_resolution),
          contact_depth_(Real(0.5) * static_cast<Real>(cfg.wall_layers) * cfg.global_resolution)
    {
        const BoundingBoxd blade_bounds = domain_stl_.surface(Rotor67DomainSurface::Blade).bounds_;
        bounding_box_ = BoundingBoxd(
            blade_bounds.lower_ - Vecd(contact_depth_, contact_depth_, contact_depth_),
            blade_bounds.upper_ + Vecd(contact_depth_, contact_depth_, contact_depth_));
    }

    bool checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED = true) override
    {
        if (!domain_stl_.bladeContainsPoint(pnt))
        {
            return false;
        }
        const Real depth = domain_stl_.bladeUnsignedDistance(pnt);
        return BOUNDARY_INCLUDED ? (depth >= first_layer_offset_ - Eps &&
                                    depth <= contact_depth_ + Eps)
                                 : (depth > first_layer_offset_ + Eps &&
                                    depth < contact_depth_ - Eps);
    }

    Vecd findClosestPoint(const Vecd &probe_point) override
    {
        return domain_stl_.closestPoint(Rotor67DomainSurface::Blade, probe_point);
    }

    BoundingBoxd findBounds() override { return bounding_box_; }

  private:
    const Rotor67DomainSTL &domain_stl_;
    Real first_layer_offset_;
    Real contact_depth_;
};

//----------------------------------------------------------------------
//  Hub wall: thin shell around the tabulated r_hub(z) meridional curve,
//  spanning the sector. The shell is r in
//  [r_hub(z) - thickness/2, r_hub(z) + thickness/2]
//  for volumetric wall particles. r_hub(z) is loaded from the CFX
//  Rot_Hub.curve file via Rotor67MeridionalCurve (Phase B+ true curved
//  geometry — replaces the Phase A constant-radius simplification).
//----------------------------------------------------------------------
class Rotor67HubShell : public Shape
{
  public:
    Rotor67HubShell(const std::string &shape_name, const Rotor67Config &cfg,
                    Real thickness, const Rotor67MeridionalCurve &curve,
                    const Rotor67DomainSTL &domain_stl)
        : Shape(shape_name),
          theta0_(domain_stl.periodicThetaMin()), delta_theta_(domain_stl.periodicDeltaTheta()),
          thickness_(thickness),
          curve_(curve),
          z_in_(domain_stl.inletZ()), z_out_(domain_stl.outletZ())
    {
        // Bounding box uses the maximum (r_hub + thickness/2) over the SPH z range.
        // We walk the curve z range, query r_hub(z) at the bracket endpoints,
        // and combine with the clamped projections. Since r_hub(z) is bounded
        // over the curve z range, sampling at the endpoints is sufficient.
        Real r_outer_max = Real(0.0);
        const Real z_lo = std::max(z_in_, curve_.hubZMin());
        const Real z_hi = std::min(z_out_, curve_.hubZMax());
        if (z_lo < z_hi)
        {
            // Sample several points across [z_lo, z_hi] for a safe bound.
            for (Real t = Real(0.0); t <= Real(1.0); t += Real(0.05))
            {
                const Real z = z_lo + t * (z_hi - z_lo);
                r_outer_max = std::max(r_outer_max, curve_.hubR(z) + Real(0.5) * thickness_);
            }
        }
        // Always consider the clamp-extended endpoint projections.
        r_outer_max = std::max(r_outer_max, curve_.hubR(z_in_) + Real(0.5) * thickness_);
        r_outer_max = std::max(r_outer_max, curve_.hubR(z_out_) + Real(0.5) * thickness_);
        bounding_box_ = BoundingBoxd(Vecd(-r_outer_max, -r_outer_max, z_in_),
                                     Vecd(r_outer_max, r_outer_max, z_out_));
    };
    ~Rotor67HubShell() {};

    bool checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED = true) override
    {
        const Real x = pnt[0], y = pnt[1], z = pnt[2];
        const Real r = std::sqrt(x * x + y * y);
        if (r < TinyReal)
            return false;
        const Real theta = std::atan2(y, x);
        const Real dtheta = normalisedAngleDiff(theta, theta0_);
        const bool in_theta = (dtheta >= -Eps && dtheta <= delta_theta_ + Eps);
        // Curved shell: r_hub(z) is tabulated, shell is r_hub(z) ± thickness/2.
        const Real r_hub = curve_.hubR(z);
        const Real r_lo = std::max(Real(0.0), r_hub - Real(0.5) * thickness_);
        const Real r_hi = r_hub + Real(0.5) * thickness_;
        const bool in_r = (r >= r_lo - Eps && r <= r_hi + Eps);
        const bool in_z = (z >= z_in_ - Eps && z <= z_out_ + Eps);
        return in_theta && in_r && in_z;
    };
    Vecd findClosestPoint(const Vecd &probe_point) override
    {
        // The hub wall is a curved thin shell around the r_hub(z) curve.
        // We project the probe onto the mid-surface r = r_hub(z) (not the
        // shell boundary) so that Shape::findNormalDirection gets a non-zero
        // surface displacement for particles sitting on either the inner or
        // outer shell face. (Returning the probe for shell-boundary hits
        // would deadlock the jitter while-loop in findNormalDirection because
        // the jitter is ~1e-14, far too small to escape the shell.)
        Real x = probe_point[0], y = probe_point[1], z = probe_point[2];
        Real r = std::sqrt(x * x + y * y);
        if (r < TinyReal)
        {
            const Real r0 = curve_.hubR(Real(0.5) * (z_in_ + z_out_));
            x = r0 * std::cos(theta0_);
            y = r0 * std::sin(theta0_);
            r = r0;
        }
        // Clamp z first so the mid-surface r is evaluated at a valid z.
        z = std::max(z_in_, std::min(z_out_, z));
        const Real r_hub = curve_.hubR(z);
        // Project to the mid-surface r = r_hub(z) (not the shell boundary).
        const Real scale = r_hub / r;
        x *= scale;
        y *= scale;
        // theta clamping to the wedge sector.
        Real theta = std::atan2(y, x);
        Real dtheta = normalisedAngleDiff(theta, theta0_);
        if (dtheta < 0.0)
        {
            x = r_hub * std::cos(theta0_);
            y = r_hub * std::sin(theta0_);
        }
        else if (dtheta > delta_theta_)
        {
            const Real tu = theta0_ + delta_theta_;
            x = r_hub * std::cos(tu);
            y = r_hub * std::sin(tu);
        }
        return Vecd(x, y, z);
    };
    BoundingBoxd findBounds() override { return bounding_box_; };

  private:
    Real normalisedAngleDiff(Real a, Real b) const
    {
        Real d = a - b;
        while (d > Pi)
            d -= Real(2.0) * Pi;
        while (d < -Pi)
            d += Real(2.0) * Pi;
        return d;
    }
    Real theta0_, delta_theta_, thickness_;
    const Rotor67MeridionalCurve &curve_;
    Real z_in_, z_out_;
};

//----------------------------------------------------------------------
//  Casing wall: thin shell around the tabulated r_shd(z) meridional curve.
//  Mirrors Rotor67HubShell with curve_.shroudR(z) in place of curve_.hubR(z).
//----------------------------------------------------------------------
class Rotor67CasingShell : public Shape
{
  public:
    Rotor67CasingShell(const std::string &shape_name, const Rotor67Config &cfg,
                       Real thickness, const Rotor67MeridionalCurve &curve,
                       const Rotor67DomainSTL &domain_stl)
        : Shape(shape_name),
          theta0_(domain_stl.periodicThetaMin()), delta_theta_(domain_stl.periodicDeltaTheta()),
          thickness_(thickness),
          curve_(curve),
          z_in_(domain_stl.inletZ()), z_out_(domain_stl.outletZ())
    {
        Real r_outer_max = Real(0.0);
        const Real z_lo = std::max(z_in_, curve_.shroudZMin());
        const Real z_hi = std::min(z_out_, curve_.shroudZMax());
        if (z_lo < z_hi)
        {
            for (Real t = Real(0.0); t <= Real(1.0); t += Real(0.05))
            {
                const Real z = z_lo + t * (z_hi - z_lo);
                r_outer_max = std::max(r_outer_max, curve_.shroudR(z) + Real(0.5) * thickness_);
            }
        }
        r_outer_max = std::max(r_outer_max, curve_.shroudR(z_in_) + Real(0.5) * thickness_);
        r_outer_max = std::max(r_outer_max, curve_.shroudR(z_out_) + Real(0.5) * thickness_);
        bounding_box_ = BoundingBoxd(Vecd(-r_outer_max, -r_outer_max, z_in_),
                                     Vecd(r_outer_max, r_outer_max, z_out_));
    };
    ~Rotor67CasingShell() {};

    bool checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED = true) override
    {
        const Real x = pnt[0], y = pnt[1], z = pnt[2];
        const Real r = std::sqrt(x * x + y * y);
        if (r < TinyReal)
            return false;
        const Real theta = std::atan2(y, x);
        const Real dtheta = normalisedAngleDiff(theta, theta0_);
        const bool in_theta = (dtheta >= -Eps && dtheta <= delta_theta_ + Eps);
        const Real r_shd = curve_.shroudR(z);
        const Real r_lo = std::max(Real(0.0), r_shd - Real(0.5) * thickness_);
        const Real r_hi = r_shd + Real(0.5) * thickness_;
        const bool in_r = (r >= r_lo - Eps && r <= r_hi + Eps);
        const bool in_z = (z >= z_in_ - Eps && z <= z_out_ + Eps);
        return in_theta && in_r && in_z;
    };
    Vecd findClosestPoint(const Vecd &probe_point) override
    {
        // See Rotor67HubShell::findClosestPoint for the mid-surface projection
        // rationale (avoids the findNormalDirection jitter deadlock).
        Real x = probe_point[0], y = probe_point[1], z = probe_point[2];
        Real r = std::sqrt(x * x + y * y);
        if (r < TinyReal)
        {
            const Real r0 = curve_.shroudR(Real(0.5) * (z_in_ + z_out_));
            x = r0 * std::cos(theta0_);
            y = r0 * std::sin(theta0_);
            r = r0;
        }
        z = std::max(z_in_, std::min(z_out_, z));
        const Real r_shd = curve_.shroudR(z);
        const Real scale = r_shd / r;
        x *= scale;
        y *= scale;
        Real theta = std::atan2(y, x);
        Real dtheta = normalisedAngleDiff(theta, theta0_);
        if (dtheta < 0.0)
        {
            x = r_shd * std::cos(theta0_);
            y = r_shd * std::sin(theta0_);
        }
        else if (dtheta > delta_theta_)
        {
            const Real tu = theta0_ + delta_theta_;
            x = r_shd * std::cos(tu);
            y = r_shd * std::sin(tu);
        }
        return Vecd(x, y, z);
    };
    BoundingBoxd findBounds() override { return bounding_box_; };

  private:
    Real normalisedAngleDiff(Real a, Real b) const
    {
        Real d = a - b;
        while (d > Pi)
            d -= Real(2.0) * Pi;
        while (d < -Pi)
            d += Real(2.0) * Pi;
        return d;
    }
    Real theta0_, delta_theta_, thickness_;
    const Rotor67MeridionalCurve &curve_;
    Real z_in_, z_out_;
};

//----------------------------------------------------------------------
//  Initial condition: uniform inlet static state + annular-swirl relative
//  velocity w(pos) = (Omega*y, -Omega*x, v_x). Density/pressure from the
//  isentropic inlet static conditions. E = Vol*(rho_e + 0.5*rho*|w|^2).
//----------------------------------------------------------------------
class Rotor67InitialCondition : public fluid_dynamics::CompressibleFluidInitialCondition
{
  public:
    explicit Rotor67InitialCondition(SPHBody &sph_body, const Rotor67Config &cfg,
                                     const Rotor67MeridionalCurve &curve)
        : fluid_dynamics::CompressibleFluidInitialCondition(sph_body),
          cfg_(cfg),
          curve_(curve),
          gamma_(cfg.gamma),
          // registerStateVariableData is idempotent.
          mass_(particles_->registerStateVariableData<Real>("Mass")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")) {};

    void update(size_t index_i, Real dt)
    {
        const Rotor67InletPrimitiveState state =
            rotor67InletStaticState(pos_[index_i], cfg_, curve_);
        rho_[index_i] = state.rho_static;
        p_[index_i] = state.p_static;
        vel_[index_i] = state.velocity_relative;
        mass_[index_i] = rho_[index_i] * Vol_[index_i];
        mom_[index_i] = mass_[index_i] * vel_[index_i];
        const Real rho_e = p_[index_i] / (gamma_ - Real(1.0));
        E_[index_i] = rho_e * Vol_[index_i] +
                      Real(0.5) * mass_[index_i] * vel_[index_i].squaredNorm();
    }

  protected:
    const Rotor67Config &cfg_;
    const Rotor67MeridionalCurve &curve_;
    Real gamma_;
    Real *mass_, *Vol_;
};

//----------------------------------------------------------------------
//  Finiteness diagnostics (mirrors eulerian_channel_geometry.hpp).
//----------------------------------------------------------------------
inline bool rotorIsFiniteReal(Real value)
{
    return std::isfinite(static_cast<double>(value));
}
inline bool rotorIsFiniteVector(const Vecd &value)
{
    for (int axis = 0; axis != Dimensions; ++axis)
    {
        if (!rotorIsFiniteReal(value[axis]))
        {
            return false;
        }
    }
    return true;
}

/**
 * @brief Scan the fluid heptad (rho/p/mass/vel/mom/E/Vol) for finiteness.
 *        Prints the first bad particle record. Returns true if all finite.
 */
inline bool reportRotorFiniteState(BaseParticles &particles, const std::string &stage,
                                   size_t *first_bad_particle = nullptr)
{
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *rho = particles.getVariableDataByName<Real>("Density");
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    Real *pressure = particles.getVariableDataByName<Real>("Pressure");
    Real *E = particles.getVariableDataByName<Real>("TotalEnergy");
    Vecd *velocity = particles.getVariableDataByName<Vecd>("Velocity");
    Vecd *momentum = particles.getVariableDataByName<Vecd>("Momentum");
    const size_t total_particles = particles.TotalRealParticles();

    Real rho_min = std::numeric_limits<Real>::max();
    Real rho_max = -std::numeric_limits<Real>::max();
    Real p_min = std::numeric_limits<Real>::max();
    Real p_max = -std::numeric_limits<Real>::max();
    Real speed_max = Real(0.0);
    size_t first_bad = total_particles;

    for (size_t i = 0; i != total_particles; ++i)
    {
        const bool finite =
            rotorIsFiniteVector(pos[i]) && rotorIsFiniteReal(vol[i]) &&
            rotorIsFiniteReal(rho[i]) && rotorIsFiniteReal(mass[i]) &&
            rotorIsFiniteReal(pressure[i]) && rotorIsFiniteReal(E[i]) &&
            rotorIsFiniteVector(velocity[i]) && rotorIsFiniteVector(momentum[i]);
        if (!finite && first_bad == total_particles)
        {
            first_bad = i;
        }
        if (rotorIsFiniteReal(rho[i]))
        {
            rho_min = SMIN(rho_min, rho[i]);
            rho_max = SMAX(rho_max, rho[i]);
        }
        if (rotorIsFiniteReal(pressure[i]))
        {
            p_min = SMIN(p_min, pressure[i]);
            p_max = SMAX(p_max, pressure[i]);
        }
        if (rotorIsFiniteVector(velocity[i]))
        {
            speed_max = SMAX(speed_max, velocity[i].norm());
        }
    }

    if (first_bad != total_particles)
    {
        std::cout << "[RotorDiag][State] " << stage
                  << ": rho=[" << rho_min << ", " << rho_max << "]"
                  << ", p=[" << p_min << ", " << p_max << "]"
                  << ", max|u|=" << speed_max
                  << ", first_bad_particle=" << first_bad << std::endl;
        std::cout << "[RotorDiag][BadParticle] " << stage
                  << ": i=" << first_bad
                  << ", pos=(" << pos[first_bad][0] << ", " << pos[first_bad][1] << ", "
                  << pos[first_bad][2] << ")"
                  << ", Vol=" << vol[first_bad]
                  << ", rho=" << rho[first_bad]
                  << ", mass=" << mass[first_bad]
                  << ", p=" << pressure[first_bad]
                  << ", E=" << E[first_bad]
                  << ", vel=(" << velocity[first_bad][0] << ", " << velocity[first_bad][1]
                  << ", " << velocity[first_bad][2] << ")"
                  << ", mom=(" << momentum[first_bad][0] << ", " << momentum[first_bad][1]
                  << ", " << momentum[first_bad][2] << ")" << std::endl;
    }
    if (first_bad_particle != nullptr)
    {
        *first_bad_particle = first_bad;
    }
    return first_bad == total_particles;
}

} // namespace SPH

#endif // ROTOR67_GEOMETRY_H
