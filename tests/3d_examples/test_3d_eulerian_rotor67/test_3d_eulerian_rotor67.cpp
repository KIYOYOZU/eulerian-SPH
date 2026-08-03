/**
 * @file 	test_3d_eulerian_rotor67.cpp
 * @brief 	Placeholder. The Rotor67 case is being rewritten.
 * @details Geometry construction (Rotor67DomainSTL, meridional curve, fluid
 *          volume and wall shells) is preserved in rotor67_geometry.hpp,
 *          and its configuration/INI parsing in rotor67_data.hpp. The
 *          previous force/boundary layers and the main integration loop
 *          were removed together with the shared compressible/rotating
 *          src helpers they depended on; they will be rewritten from
 *          scratch on top of the preserved geometry layer.
 *          This file is intentionally left as a GTest placeholder so the
 *          target still builds and registers with CTest during the rewrite.
 * @author 	KIYOYOZU
 */
#include <gtest/gtest.h>

// Placeholder test. The real Rotor67 main loop will be rewritten on top of
// rotor67_geometry.hpp / rotor67_data.hpp.
TEST(Rotor67Placeholder, UnderConstruction)
{
    SUCCEED();
}
