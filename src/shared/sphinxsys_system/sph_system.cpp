#include "sph_system.hpp"

#include "all_body_relations.h"
#include "io_log.h"
#include "predefined_bodies.h"

namespace SPH
{
//=================================================================================================//
SPHSystem::SPHSystem(BoundingBoxd system_domain_bounds, Real global_resolution, size_t number_of_threads)
    : SPHSystem(true, system_domain_bounds, global_resolution, number_of_threads) {}
//=================================================================================================//
SPHSystem::SPHSystem(bool is_physical, BoundingBoxd system_domain_bounds,
                     Real global_resolution, size_t number_of_threads)
    : system_domain_bounds_(system_domain_bounds),
      global_resolution_(global_resolution),
      tbb_global_control_(tbb::global_control::max_allowed_parallelism, number_of_threads),
      is_physical_(is_physical),
      io_environment_(io_keeper_.createPtr<IOEnvironment>(*this)),
      run_particle_relaxation_(false), reload_particles_(false),
      restart_step_(0), generate_regression_data_(false), state_recording_(true)
{
    Log::init();
    sv_physical_time_ = registerSystemVariable<Real>("PhysicalTime", 0.0);
    Log::get()->info("The reference resolution of the SPHSystem is " + std::to_string(global_resolution_) + ".");
}
//=================================================================================================//
void SPHSystem::setLogLevel(size_t log_level)
{
    if (log_level > 6)
    {
        std::cerr << "Log level must be between 0 and 6.\n";
        exit(1);
    }
    log_level_ = log_level;
    Log::get()->setLevel(static_cast<int>(log_level_));
}
//=================================================================================================//
IOEnvironment &SPHSystem::getIOEnvironment()
{
    checkPointer(io_environment_, "io_environment_", "SPHSystem");
    return *io_environment_;
}
//=================================================================================================//
void SPHSystem::addRealBody(RealBody *real_body)
{
    real_bodies_.push_back(real_body);
}
//=================================================================================================//

void SPHSystem::initializeSystemCellLinkedLists()
{
    for (auto &body : real_bodies_)
    {
        DynamicCast<RealBody>(this, body)->updateCellLinkedList();
    }
}
//=================================================================================================//
void SPHSystem::initializeSystemConfigurations()
{
    for (auto &body : sph_bodies_)
    {
        StdVec<SPHRelation *> &body_relations = body->getBodyRelations();
        for (size_t i = 0; i < body_relations.size(); i++)
        {
            body_relations[i]->updateConfiguration();
        }
    }
}
//=================================================================================================//
SPHSystem *SPHSystem::handleCommandlineOptions(int ac, char *av[])
{
    for (int i = 1; i < ac; ++i)
    {
        std::string arg(av[i]);
        if (arg == "--relax=true" || arg == "--relax=1")
        {
            run_particle_relaxation_ = true;
            std::cout << "Particle relaxation set to true.\n";
        }
        else if (arg == "--relax=false" || arg == "--relax=0")
            run_particle_relaxation_ = false;
        else if (arg == "--reload=true" || arg == "--reload=1")
        {
            reload_particles_ = true;
            std::cout << "Particle reload set to true.\n";
        }
        else if (arg == "--reload=false" || arg == "--reload=0")
            reload_particles_ = false;
        else if (arg == "--regression=true" || arg == "--regression=1")
            generate_regression_data_ = true;
        else if (arg.substr(0, 15) == "--restart_step=")
            restart_step_ = static_cast<size_t>(std::stoi(arg.substr(15)));
        else if (arg == "--help")
        {
            std::cout << "Options: --relax=<bool> --reload=<bool> --regression=<bool> --restart_step=<int>\n";
            exit(0);
        }
    }
    return this;
}
//=================================================================================================//
} // namespace SPH
