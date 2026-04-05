/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,           *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file    io_log.h
 * @brief   Lightweight logger (no spdlog/fmt dependency).
 * @author  Chi Zhang, Shuoguo Zhang, Zhenxi Zhao and Xiangyu Hu
 */

#ifndef IO_LOG_H
#define IO_LOG_H

#include "io_base.h"

#include <memory>
#include <string>

namespace SPH
{
namespace Log
{

/**
 * @brief Minimal logger that writes to stdout and an optional rotating log file.
 *        Replaces spdlog to remove the fmt/spdlog dependency.
 */
class SimpleLogger
{
  public:
    explicit SimpleLogger(const std::string &filename = "sphinxsys.log");
    ~SimpleLogger();

    template <typename... Args>
    void info(const std::string &msg, Args &&...) { write("INFO", msg); }

    template <typename... Args>
    void warn(const std::string &msg, Args &&...) { write("WARN", msg); }

    template <typename... Args>
    void error(const std::string &msg, Args &&...) { write("ERROR", msg); }

    void setLevel(int /*level*/) {} // no-op: kept for API compatibility

  private:
    void write(const char *level, const std::string &msg);
    std::string filename_;
};

/** Initialize the global logger (call once at startup). */
std::shared_ptr<SimpleLogger> init();

/** Access the global logger. */
std::shared_ptr<SimpleLogger> get();

} // namespace Log
} // namespace SPH

#endif // IO_LOG_H
