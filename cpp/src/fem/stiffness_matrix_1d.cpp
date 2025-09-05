#include <Eigen/Sparse>
#include <vector>
#include <cmath>
using namespace Eigen;

Matrix<double, Dynamic, Dynamic> der_shape_func_deg1(VectorXd &element_nodes, int shape_func_idx) {

    switch (shape_func_idx)
    {
    case 0:
        
        break;
    
    default:
        break;
    }
} 

/**
 * @brief Computes the local stiffness matrix for a 1D finite element.
 * 
 * This function calculates the local stiffness matrix for a single element
 * using numerical integration (quadrature) based on the degree of freedom (dof).
 * 
 * @param element_nodes A vector containing the coordinates of the nodes in the element.
 * @param dof The degree of freedom, i.e., the number of nodes per element.
 * @return Matrix<double, Dynamic, Dynamic> The local stiffness matrix for the element.
 */
Matrix<double, Dynamic, Dynamic> element_stiffness_matrix_1d(const VectorXd &element_nodes, const int &dof) {

    // Initializations
    MatrixXd K_loc(dof, dof);
    MatrixXd shape_funct_val(dof, dof);
    RowVectorXd quad_nodes(dof);
    RowVectorXd quad_weights(dof); 

    switch (dof)
    {
    case 2:
        quad_nodes << -1, 1;
        quad_weights << 1, 1;
        shape_funct_val.col(0) = (-0.5)*Vector2d::Ones();
        shape_funct_val.col(1) = 0.5*Vector2d::Ones();
        break;
    
    default:
        break;
    }
    
    double element_size = std::abs(element_nodes(0) - element_nodes(element_nodes.size() - 1));
    for (int i = 0; i < dof; i++) {
        for (int j = 0; j <= i; j++) {
            K_loc(i,j) = quad_weights*(shape_funct_val.col(i).array() * shape_funct_val.col(j).array()).matrix();
            K_loc(j,i) = K_loc(i,j);
        }
    }
    K_loc *= 2/element_size;
    return K_loc;
}


/**
 * @brief Computes the local stiffness matrix for a 1D finite element.
 * 
 * This function calculates the local stiffness matrix for a single element
 * using numerical integration (quadrature) based on the degree of freedom (dof).
 * 
 * @param element_nodes A vector containing the coordinates of the nodes in the element.
 * @param dof The degree of freedom, i.e., the number of nodes per element.
 * @return Matrix<double, Dynamic, Dynamic> The local stiffness matrix for the element.
 */
SparseMatrix<double> stiffness_matrix_1d(VectorXd &nodes, Matrix<int, Dynamic, Dynamic> &elements) {

    // Initializations
    int num_nodes = nodes.size();
    int num_of_elements = elements.rows();
    const int dof = elements.cols();
    SparseMatrix<double> K(num_nodes, num_nodes);
    typedef Triplet<double> T;
    std::vector<T> triplet_list;
    triplet_list.reserve(dof * dof * num_of_elements);

    for (int k = 0; k < num_of_elements; k++) {

        // collect element nodes
        VectorXd element_nodes; 
        for (int i = 0; i < dof; i++) {
            element_nodes(i) = nodes(elements(k,i));
        }

        // get local element matrix and write into K
        MatrixXd K_loc = element_stiffness_matrix_1d(element_nodes, dof); 
        for (int i = 0; i < K_loc.rows(); i++) {
            for (int j = 0; j < K_loc.cols(); j++) {
                triplet_list.push_back(T(elements(k,i), elements(k,j), K_loc(i,j)));
            }   
        }
    }
    K.setFromTriplets(triplet_list.begin(), triplet_list.end());
    return K;
}




