cmake_minimum_required(VERSION 3.20)

project(eulerian_sph LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(EULERIAN_SPH_BUILD_CASES "Build standalone Eulerian-SPH example cases." ON)
find_package(Eigen3 REQUIRED)
find_package(TBB REQUIRED)
find_package(Boost QUIET)  # boost::geometry (header-only) required by MultiPolygonShape

add_subdirectory(src)

if(EULERIAN_SPH_BUILD_CASES)
    add_subdirectory(cases)
endif()

