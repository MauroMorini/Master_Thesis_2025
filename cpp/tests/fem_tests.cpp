#include <catch2/catch_test_macros.hpp>
#include <Eigen/Sparse>
#include "stiffness_matrix_1d.hpp"

using namespace Eigen;

// test local stiffness matrix for linear elements
TEST_CASE("Element stiffness matrix for element [0,1]", "[el_stiff]") {
    const int dof = 2;
    Vector2d element(0.0, 1.0);
    MatrixXd K_loc_test = element_stiffness_matrix_1d(element, dof);
    MatrixXd K_loc_ref(2,2);
    K_loc_ref << 1, -1, 
                -1, 1;
    REQUIRE(K_loc_test.isApprox(K_loc_ref));
}