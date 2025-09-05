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

// test assembly of trivial stiffness matrix linear elements
TEST_CASE("global stiffness matrix for element [0, 1], N=6", "[glob_stiff_triv]") {
    VectorXd nodes(6);
    nodes << 0.0, 0.2, 0.4, 0.6, 0.8, 1.0;
    Matrix<int, Dynamic, Dynamic> elements(5,2); 
    elements << 1, 2, 2, 3, 3, 4 ,4, 5, 5, 6;
    MatrixXd K_ref(6,6);
    K_ref.setZero();
    K_ref(0,0) = 5.0;
    K_ref(0,1) = -5.0;
    K_ref(1,0) = -5.0;
    K_ref(1,1) = 10.0;
    K_ref(1,2) = -5.0;
    K_ref(2,1) = -5.0;
    K_ref(2,2) = 5.0;
    K_ref(3,3) = 5.0;
    K_ref(3,4) = -5.0;
    K_ref(4,3) = -5.0;
    K_ref(4,4) = 10.0;
    K_ref(4,5) = -5.0;
    K_ref(5,4) = -5.0;
    K_ref(5,5) = 5.0;
    // TODO: throws exception!
    SparseMatrix<double> K_test = stiffness_matrix_1d(nodes, elements);
    REQUIRE(K_ref.isApprox(K_test.toDense()));
}