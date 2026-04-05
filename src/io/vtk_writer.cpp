#include "vtk_writer.h"

#include <fstream>
#include <stdexcept>

namespace eulerian_sph
{
namespace io
{
void writeVtkFrame(const SimulationContext &context, const std::filesystem::path &file_path)
{
    std::ofstream output(file_path, std::ios::out | std::ios::trunc);
    if (!output.is_open())
    {
        throw std::runtime_error("Failed to open VTK output file: " + file_path.string());
    }

    const std::size_t particle_count = context.particles.size();

    output << "<?xml version=\"1.0\"?>\n";
    output << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    output << "  <PolyData>\n";
    output << "    <Piece NumberOfPoints=\"" << particle_count
           << "\" NumberOfVerts=\"" << particle_count
           << "\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
    output << "      <PointData Scalars=\"Density\">\n";
    output << "        <DataArray type=\"Float64\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (const ParticleState &particle : context.particles.data())
    {
        output << "          " << particle.rho << "\n";
    }
    output << "        </DataArray>\n";
    output << "        <DataArray type=\"Float64\" Name=\"Pressure\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (const ParticleState &particle : context.particles.data())
    {
        output << "          " << particle.p << "\n";
    }
    output << "        </DataArray>\n";
    output << "        <DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const ParticleState &particle : context.particles.data())
    {
        output << "          " << particle.vel.x() << ' ' << particle.vel.y() << " 0.0\n";
    }
    output << "        </DataArray>\n";
    output << "        <DataArray type=\"Int32\" Name=\"Indicator\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (const ParticleState &particle : context.particles.data())
    {
        output << "          " << particle.indicator << "\n";
    }
    output << "        </DataArray>\n";
    output << "        <DataArray type=\"Int32\" Name=\"SmearedSurface\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (const ParticleState &particle : context.particles.data())
    {
        output << "          " << particle.smeared_surface << "\n";
    }
    output << "        </DataArray>\n";
    output << "        <DataArray type=\"Float64\" Name=\"NormalDirection\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const ParticleState &particle : context.particles.data())
    {
        output << "          " << particle.normal.x() << ' ' << particle.normal.y() << " 0.0\n";
    }
    output << "        </DataArray>\n";
    output << "      </PointData>\n";
    output << "      <Points>\n";
    output << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const ParticleState &particle : context.particles.data())
    {
        output << "          " << particle.pos.x() << ' ' << particle.pos.y() << " 0.0\n";
    }
    output << "        </DataArray>\n";
    output << "      </Points>\n";
    output << "      <Verts>\n";
    output << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    output << "          ";
    for (std::size_t i = 0; i != particle_count; ++i)
    {
        output << i << ' ';
    }
    output << "\n        </DataArray>\n";
    output << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    output << "          ";
    for (std::size_t i = 0; i != particle_count; ++i)
    {
        output << (i + 1) << ' ';
    }
    output << "\n        </DataArray>\n";
    output << "      </Verts>\n";
    output << "    </Piece>\n";
    output << "  </PolyData>\n";
    output << "</VTKFile>\n";
}
} // namespace io
} // namespace eulerian_sph
