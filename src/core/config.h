#pragma once

#include "types.h"

namespace eulerian_sph
{
struct SpatialConfig
{
    Vec2 domain_lower = Vec2::Zero();
    Vec2 domain_upper = Vec2::Ones();
    Scalar particle_spacing = 0.02;
    Scalar smoothing_length_ratio = 1.3;
    Scalar ghost_band_width = 0.0;
};

struct MaterialConfig
{
    Scalar reference_density = 1.0;
    Scalar sound_speed = 1.0;
    Scalar dynamic_viscosity = 0.0;
    Scalar reference_pressure = 0.0;
    Scalar heat_capacity_ratio = 1.4;
};

struct TimeControl
{
    Scalar end_time = 1.0;
    Scalar output_interval = 0.0;
    Scalar acoustic_cfl = 0.25;
    std::size_t max_steps = 1000;
    std::size_t screen_output_interval = 0;
};

struct OutputConfig
{
    std::filesystem::path directory = "output";
    bool write_csv = false;
    bool write_vtp = true;
    std::size_t vtp_output_count = 10;
};

struct SimulationConfig
{
    std::string case_name = "unnamed_case";
    SpatialConfig spatial;
    MaterialConfig material;
    TimeControl time;
    OutputConfig output;
};
} // namespace eulerian_sph
