cmake_minimum_required(VERSION 3.13)
project(test_3d_eulerian_channel)

# 尝试查找系统安装的GTest（超算环境）
find_package(GTest QUIET)

# 如果没找到，使用FetchContent下载（本地环境）
if(NOT GTest_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://gitee.com/mirrors/googletest/repository/archive/v1.14.0.zip
  )
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
endif()

include(GoogleTest)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin/")

# 主测试程序（包含GTest）
add_executable(${PROJECT_NAME} test_3d_eulerian_channel.cpp)
if(GTest_FOUND)
  target_link_libraries(${PROJECT_NAME} sphinxsys_3d GTest::gtest GTest::gtest_main)
else()
  target_link_libraries(${PROJECT_NAME} sphinxsys_3d gtest gtest_main)
endif()
set_target_properties(${PROJECT_NAME} PROPERTIES
  OUTPUT_NAME "eulerian_channel"
  RUNTIME_OUTPUT_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}"
  VS_DEBUGGER_WORKING_DIRECTORY "${EXECUTABLE_OUTPUT_PATH}")

# 将 config.ini 拷贝到运行目录，便于直接执行（resolve_default_config_path 也会回退到源目录）
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          "${CMAKE_CURRENT_SOURCE_DIR}/config.ini"
          "${EXECUTABLE_OUTPUT_PATH}/config.ini")

# GTest测试配置
add_test(NAME ${PROJECT_NAME}
    COMMAND ${PROJECT_NAME}
    WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})

gtest_discover_tests(${PROJECT_NAME}
    WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
    DISCOVERY_MODE PRE_TEST)
