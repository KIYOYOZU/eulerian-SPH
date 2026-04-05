#include "dat_recorder.h"

#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace eulerian_sph
{
namespace io
{
DatRecorder::DatRecorder(std::filesystem::path file_path, std::string quantity_name, SampleFunction sampler,
                         std::vector<std::string> column_names)
    : file_path_(std::move(file_path)),
      quantity_name_(std::move(quantity_name)),
      sampler_(std::move(sampler)),
      column_names_(std::move(column_names)) {}

std::vector<std::string> DatRecorder::buildColumnNames(const std::size_t value_count) const
{
    if (!column_names_.empty())
    {
        return column_names_;
    }

    if (value_count == 1)
    {
        return {quantity_name_};
    }

    std::vector<std::string> names;
    names.reserve(value_count);
    for (std::size_t i = 0; i != value_count; ++i)
    {
        names.push_back(quantity_name_ + "[" + std::to_string(i) + "]");
    }
    return names;
}

void DatRecorder::write(const SimulationContext &context)
{
    const std::vector<Scalar> values = sampler_(context);
    const std::vector<std::string> names = buildColumnNames(values.size());

    if (!header_written_)
    {
        std::ofstream header(file_path_, std::ios::out | std::ios::trunc);
        if (!header.is_open())
        {
            throw std::runtime_error("Failed to open DAT header file: " + file_path_.string());
        }

        header << "\"run_time\"   ";
        for (const std::string &name : names)
        {
            header << "\"" << name << "\"   ";
        }
        header << "\n";
        header_written_ = true;
    }

    std::ofstream output(file_path_, std::ios::out | std::ios::app);
    if (!output.is_open())
    {
        throw std::runtime_error("Failed to open DAT file: " + file_path_.string());
    }

    output << std::fixed << std::setprecision(9) << context.current_time << "   ";
    for (const Scalar value : values)
    {
        output << std::fixed << std::setprecision(9) << value << "   ";
    }
    output << "\n";
}
} // namespace io
} // namespace eulerian_sph
