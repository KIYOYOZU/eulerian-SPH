#pragma once

#include "../core/simulation_context.h"

namespace eulerian_sph
{
namespace io
{
void writeCsvFrame(const SimulationContext &context, const std::filesystem::path &file_path);
} // namespace io
} // namespace eulerian_sph
