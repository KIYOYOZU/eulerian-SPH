set(EULERIAN_SPH_SOURCES
    core/simulation.cpp
    io/csv_writer.cpp
    io/dat_recorder.cpp
    io/vtk_writer.cpp
)

add_library(eulerian_sph_core ${EULERIAN_SPH_SOURCES})
add_library(eulerian_sph::core ALIAS eulerian_sph_core)

target_include_directories(eulerian_sph_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(eulerian_sph_core
    PUBLIC
        Eigen3::Eigen
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(eulerian_sph_core PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(eulerian_sph_core PUBLIC EULERIAN_SPH_USE_OPENMP=1)
endif()
