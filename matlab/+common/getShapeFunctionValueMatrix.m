function [phi_val, dphi_val, quad_weights] = getShapeFunctionValueMatrix(dof)
    % calculates matrix which contains the values of the shape functions
    % and their derivatives at the quadrature nodes
    %
    % Inputs:      
    %       dof:            scalar degree of freedom
    % Outputs:
    %       phi_val:        (dof, dof) matrix with entry (i,j) phi_i(x_j)
    %       dphi_val:       (dof, dof) matrix with entry (i,j) dphi_i(x_j)
    %       quad_weights:   (1, dof) weight vector
    phi_val = zeros(dof,dof);
    dphi_val = zeros(dof, dof);
    [quad_nodes, quad_weights] = common.QuadratureFEM.getLobattoQuadrature(dof);
    [phi_cell, dphi_cell] = common.QuadratureFEM.getLagrangeBasisFun(dof);
    for i = 1:dof
        phi_val(i,:) = phi_cell{i}(quad_nodes);
        dphi_val(i,:) = dphi_cell{i}(quad_nodes);
    end
end