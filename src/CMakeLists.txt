# -----------------------------------------------------------------------
# eulerian-SPH 核心库：shared + for_2D_build，嵌入 tinyxml2 和 kernel_tabulated_ck
# 依赖：Eigen3 + TBB + Boost(header-only for geometry)
# -----------------------------------------------------------------------

# 3rd-party: tinyxml2 (内嵌，无需外部安装)
add_library(tinyxml2_embedded STATIC
    3rd_party/tinyxml2/tinyxml2.cpp)
target_include_directories(tinyxml2_embedded PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/tinyxml2>)
target_compile_features(tinyxml2_embedded PUBLIC cxx_std_17)
if(MSVC)
    target_compile_options(tinyxml2_embedded PRIVATE /utf-8)
endif()

# 收集 shared 源文件，排除 simbody 和 shared_ck
file(GLOB_RECURSE SPHINXSYS_SHARED_HEADERS CONFIGURE_DEPENDS
    shared/*.h shared/*.hpp shared/*.hxx)
file(GLOB_RECURSE SPHINXSYS_SHARED_SOURCES CONFIGURE_DEPENDS
    shared/*.cpp)

# 收集 for_2D_build 源文件，排除 simbody
file(GLOB_RECURSE SPHINXSYS_2D_HEADERS CONFIGURE_DEPENDS
    for_2D_build/*.h for_2D_build/*.hpp)
file(GLOB_RECURSE SPHINXSYS_2D_SOURCES CONFIGURE_DEPENDS
    for_2D_build/*.cpp)

# 排除 simbody 相关（防止有残留）
list(FILTER SPHINXSYS_SHARED_HEADERS EXCLUDE REGEX ".*simbody.*")
list(FILTER SPHINXSYS_SHARED_SOURCES EXCLUDE REGEX ".*simbody.*")
list(FILTER SPHINXSYS_2D_HEADERS EXCLUDE REGEX ".*simbody.*")
list(FILTER SPHINXSYS_2D_SOURCES EXCLUDE REGEX ".*simbody.*")

add_library(sphinxsys_2d)
add_library(SPHinXsys::sphinxsys_2d ALIAS sphinxsys_2d)

target_sources(sphinxsys_2d
    PUBLIC  ${SPHINXSYS_SHARED_HEADERS} ${SPHINXSYS_2D_HEADERS}
    PRIVATE ${SPHINXSYS_SHARED_SOURCES} ${SPHINXSYS_2D_SOURCES}
)

# 收集所有头文件目录作为 include 路径
foreach(FILEPATH ${SPHINXSYS_SHARED_HEADERS})
    get_filename_component(DIR ${FILEPATH} DIRECTORY)
    list(APPEND SPHINXSYS_2D_INCLUDE_DIRS ${DIR})
endforeach()
foreach(FILEPATH ${SPHINXSYS_2D_HEADERS})
    get_filename_component(DIR ${FILEPATH} DIRECTORY)
    list(APPEND SPHINXSYS_2D_INCLUDE_DIRS ${DIR})
endforeach()
list(REMOVE_DUPLICATES SPHINXSYS_2D_INCLUDE_DIRS)

foreach(DIR ${SPHINXSYS_2D_INCLUDE_DIRS})
    target_include_directories(sphinxsys_2d PUBLIC $<BUILD_INTERFACE:${DIR}>)
endforeach()

# 链接依赖
target_link_libraries(sphinxsys_2d PUBLIC
    Eigen3::Eigen
    TBB::tbb
    TBB::tbbmalloc
    Boost::boost
    Boost::program_options
    tinyxml2_embedded
)

# 编译选项
target_compile_features(sphinxsys_2d PUBLIC cxx_std_20)
if(MSVC)
    target_compile_definitions(sphinxsys_2d PUBLIC
        _USE_MATH_DEFINES
        NOMINMAX
        _SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING
    )
    target_compile_options(sphinxsys_2d PUBLIC /MP /permissive- /utf-8)
endif()
target_compile_definitions(sphinxsys_2d PUBLIC
    SPHINXSYS_USE_SYCL=0
    SPHINXSYS_USE_FLOAT=0
    SPHINXSYS_USE_ONEDPL_SORTING=0
)
