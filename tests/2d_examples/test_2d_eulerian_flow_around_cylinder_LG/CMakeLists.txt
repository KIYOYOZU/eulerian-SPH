STRING(REGEX REPLACE ".*/(.*)" "\\1" CURRENT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR})
PROJECT("${CURRENT_FOLDER}")

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin/")
SET(BUILD_INPUT_PATH "${EXECUTABLE_OUTPUT_PATH}/input")
SET(BUILD_RELOAD_PATH "${EXECUTABLE_OUTPUT_PATH}/reload")

file(MAKE_DIRECTORY ${BUILD_INPUT_PATH})
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/regression_test_tool")
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/regression_test_tool/
         DESTINATION ${BUILD_INPUT_PATH})
endif()

add_executable(${PROJECT_NAME})
target_sources(${PROJECT_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/2d_eulerian_flow_around_cylinder_LG.cpp")
target_link_libraries(${PROJECT_NAME} sphinxsys_2d)
set_target_properties(${PROJECT_NAME} PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")

add_executable(${PROJECT_NAME}_remap "${CMAKE_CURRENT_SOURCE_DIR}/local_refinement_remap.cpp")
target_link_libraries(${PROJECT_NAME}_remap sphinxsys_2d)
set_target_properties(${PROJECT_NAME}_remap PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")

add_test(NAME ${PROJECT_NAME}_particle_relaxation
     COMMAND ${PROJECT_NAME} --relax=true --state_recording=${TEST_STATE_RECORDING}
     WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
add_test(NAME ${PROJECT_NAME}
     COMMAND ${PROJECT_NAME} --reload=true --state_recording=${TEST_STATE_RECORDING}
     WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})

set_tests_properties(${PROJECT_NAME}_particle_relaxation PROPERTIES LABELS "Eulerian")
set_tests_properties(${PROJECT_NAME} PROPERTIES DEPENDS "${PROJECT_NAME}_particle_relaxation")
set_tests_properties(${PROJECT_NAME} PROPERTIES LABELS "Eulerian")

