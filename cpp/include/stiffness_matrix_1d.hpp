#pragma once

#include <Eigen/Sparse>
#include <vector>
#include <cmath>
using namespace Eigen;

Matrix<double, Dynamic, Dynamic> der_shape_func_deg1(VectorXd &element_nodes, int shape_func_idx);

Matrix<double, Dynamic, Dynamic> element_stiffness_matrix_1d(const VectorXd &element_nodes, const int &dof);

SparseMatrix<double> stiffness_matrix_1d(VectorXd &nodes, Matrix<int, Dynamic, Dynamic> &elements);




