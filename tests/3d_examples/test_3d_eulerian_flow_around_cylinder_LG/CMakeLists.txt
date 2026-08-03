cmake_minimum_required(VERSION 3.13)
project(test_3d_eulerian_flow_around_cylinder_LG)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin/")

add_executable(${PROJECT_NAME})
target_sources(${PROJECT_NAME} PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/test_3d_eulerian_flow_around_cylinder_LG.cpp")
target_link_libraries(${PROJECT_NAME} sphinxsys_3d)
set_target_properties(${PROJECT_NAME} PROPERTIES
    VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/config.ini"
            "${EXECUTABLE_OUTPUT_PATH}/config.ini")
