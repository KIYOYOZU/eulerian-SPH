cmake_minimum_required(VERSION 3.20)

project(eulerian_sph LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(EULERIAN_SPH_BUILD_CASES "Build standalone Eulerian-SPH example cases." ON)
option(EULERIAN_SPH_ENABLE_OPENMP "Enable OpenMP parallel loops for runtime hotspots." ON)

find_package(Eigen3 REQUIRED)

if(EULERIAN_SPH_ENABLE_OPENMP)
    find_package(OpenMP QUIET)
endif()

add_subdirectory(src)

if(EULERIAN_SPH_BUILD_CASES)
    add_subdirectory(cases)
endif()
