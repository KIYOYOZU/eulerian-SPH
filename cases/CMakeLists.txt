add_executable(taylor_green_case
    taylor_green/main.cpp
)
target_link_libraries(taylor_green_case PRIVATE eulerian_sph::core)

add_executable(flow_around_cylinder_lg_case
    flow_around_cylinder_lg/main.cpp
)
target_link_libraries(flow_around_cylinder_lg_case PRIVATE eulerian_sph::core)

add_executable(supersonic_flow_new_bc_case
    supersonic_flow_new_bc/main.cpp
)
target_link_libraries(supersonic_flow_new_bc_case PRIVATE eulerian_sph::core)
