#pragma once

#include "../core/simulation_context.h"

#include <functional>

namespace eulerian_sph
{
namespace io
{
class DatRecorder
{
  public:
    using SampleFunction = std::function<std::vector<Scalar>(const SimulationContext &)>;

    DatRecorder(std::filesystem::path file_path, std::string quantity_name, SampleFunction sampler,
                std::vector<std::string> column_names = {});

    void write(const SimulationContext &context);

  private:
    std::vector<std::string> buildColumnNames(std::size_t value_count) const;

    std::filesystem::path file_path_;
    std::string quantity_name_;
    SampleFunction sampler_;
    std::vector<std::string> column_names_;
    bool header_written_ = false;
};
} // namespace io
} // namespace eulerian_sph
