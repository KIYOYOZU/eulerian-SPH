#include "csv_writer.h"

#include <fstream>
#include <stdexcept>

namespace eulerian_sph
{
namespace io
{
void writeCsvFrame(const SimulationContext &context, const std::filesystem::path &file_path)
{
    std::ofstream output(file_path, std::ios::out | std::ios::trunc);
    if (!output.is_open())
    {
        throw std::runtime_error("Failed to open CSV output file: " + file_path.string());
    }

    output << "x,y,rho,p,u,v,mass,volume,energy,indicator\n";
    for (const ParticleState &particle : context.particles.data())
    {
        output << particle.pos.x() << ','
               << particle.pos.y() << ','
               << particle.rho << ','
               << particle.p << ','
               << particle.vel.x() << ','
               << particle.vel.y() << ','
               << particle.mass << ','
               << particle.volume << ','
               << particle.energy << ','
               << particle.indicator << '\n';
    }
}
} // namespace io
} // namespace eulerian_sph
