function [quad_nodes, quad_weights] = getQuadrature(dof, type_string)
    % depending on dof parameter calculates quadrature points and weights
    arguments (Input)
        dof                                     % scalar degrees of freedom
        type_string string = "Gauss-Lobatto"    % type of chosen quadrature
    end
    arguments (Output)
        quad_nodes      % (1, dof) node vector
        quad_weights    % (1, dof) weight vector
    end

    switch type_string
        case "Gauss-Lobatto"
            [quad_nodes, quad_weights] = getLobattoQuadrature(dof);
        case "Gauss-Legendre"
            [quad_nodes, quad_weights] = getGaussLegendreQuadrature(dof);
        otherwise
            error("the quadrature:  --" + type_string  + "--  has not been implemented")
    end
end