/**
 * @file 	shock_tube.h
 * @brief 	2D Eulerian SPH shock tube (Lax problem), solid-wall paradigm.
 *          Single-phase compressible ideal gas, membrane at x=0.5.
 *          x ends are real SolidBody reflective walls via ContactRelation +
 *          compressible MUSCL-WithWall integrators (the project standard;
 *          the ESPH ghost-mirror + boundary_type paradigm is retired).
 *          y direction is periodic, which on a u_y == 0 slab gives the 1D solution.
 * @author 	KIYOYOZU
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real particle_spacing_ref = 1.0 / 400.0; /**< Initial reference particle spacing. */
Real BW = 4.0 * particle_spacing_ref;    /**< System-domain padding around the slab. */
Real L = 1.0;                            /**< Tube length in x. */
Real H = 10.0 * particle_spacing_ref;    /**< Tube height in y (periodic slab). */
Real wall_thickness = 3.0 * particle_spacing_ref; /**< Solid wall slab thickness. */
Real x_membrane = 0.5;                   /**< Membrane position. */
BoundingBoxd system_domain_bounds(Vec2d(-wall_thickness - BW, -BW),
                                  Vec2d(L + wall_thickness + BW, H + BW));
//----------------------------------------------------------------------
//	Material properties of the fluid (Lax problem, ideal gas gamma=1.4).
//----------------------------------------------------------------------
Real heat_capacity_ratio = 1.4;
Real rho_L = 1.0;   /**< Left density. */
Real p_L = 0.425;   /**< Left pressure. */
Real u_L = 0.0;     /**< Left velocity (x component). */
Real rho_R = 0.125; /**< Right density. */
Real p_R = 0.1;     /**< Right pressure. */
Real u_R = 0.0;     /**< Right velocity (x component). */
Real end_time = 0.2;  /**< Compare with exact solution at t=0.2. */
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
//	Case-dependent initial condition: left/right states separated by membrane.
//----------------------------------------------------------------------
class ShockTubeInitialCondition : public fluid_dynamics::CompressibleFluidInitialCondition
{
  public:
    explicit ShockTubeInitialCondition(SPHBody &sph_body)
        : fluid_dynamics::CompressibleFluidInitialCondition(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")) {};

    void update(size_t index_i, Real dt)
    {
        const Real x = pos_[index_i][0];
        if (x < x_membrane)
        {
            rho_[index_i] = rho_L;
            p_[index_i] = p_L;
            vel_[index_i][0] = u_L;
        }
        else
        {
            rho_[index_i] = rho_R;
            p_[index_i] = p_R;
            vel_[index_i][0] = u_R;
        }
        vel_[index_i][1] = 0.0;
        mass_[index_i] = rho_[index_i] * Vol_[index_i];
        mom_[index_i] = mass_[index_i] * vel_[index_i];
        Real rho_e = p_[index_i] / (gamma_ - 1.0);
        E_[index_i] = rho_e * Vol_[index_i] + 0.5 * mass_[index_i] * vel_[index_i].squaredNorm();
    }

  protected:
    Real gamma_ = heat_capacity_ratio;
    Vecd *pos_;
    Vecd *vel_;
};
