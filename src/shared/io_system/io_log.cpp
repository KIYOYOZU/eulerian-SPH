#include "io_log.h"

#include <chrono>
#include <ctime>
#include <fstream>
#include <iostream>

namespace SPH
{
namespace Log
{
//=================================================================================================//
SimpleLogger::SimpleLogger(const std::string &filename) : filename_(filename) {}
//=================================================================================================//
SimpleLogger::~SimpleLogger() {}
//=================================================================================================//
void SimpleLogger::write(const char *level, const std::string &msg)
{
    // Format: [HH:MM:SS] [LEVEL] message
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    char buf[20];
    std::strftime(buf, sizeof(buf), "%H:%M:%S", std::localtime(&t));

    std::string line = std::string("[") + buf + "] [" + level + "] " + msg + "\n";
    std::cout << line;

    // Append to log file
    std::ofstream ofs(filename_, std::ios::app);
    if (ofs)
        ofs << line;
}
//=================================================================================================//
static std::shared_ptr<SimpleLogger> g_logger;
//=================================================================================================//
std::shared_ptr<SimpleLogger> init()
{
    if (!g_logger)
        g_logger = std::make_shared<SimpleLogger>("sphinxsys.log");
    return g_logger;
}
//=================================================================================================//
std::shared_ptr<SimpleLogger> get()
{
    if (!g_logger)
        throw std::runtime_error("Logger not initialized. Call Log::init() first.");
    return g_logger;
}
//=================================================================================================//
} // namespace Log
} // namespace SPH
