cmake_minimum_required(VERSION 3.20)

project(eulerian_sph LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(EULERIAN_SPH_BUILD_TESTS "Build eulerian test cases." ON)
option(TEST_STATE_RECORDING "State recording when run Ctest" ON)

set(EULERIAN_SPH_PROJECT_DIR ${PROJECT_SOURCE_DIR})

find_package(Eigen3 REQUIRED)
find_package(TBB REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)  # program_options + boost::geometry

# Define BOOST_AVAILABLE since we always require program_options
add_compile_definitions(BOOST_AVAILABLE)

add_subdirectory(src)

if(EULERIAN_SPH_BUILD_TESTS)
    add_subdirectory(tests)
endif()
