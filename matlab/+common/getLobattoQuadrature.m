function [quad_nodes, quad_weights] = getLobattoQuadrature(dof)
    % depending on dof parameter yields quadrature points and weights, currently only yields 
    % uniformly distributed nodal basis
    arguments (Input)
        dof             % scalar degrees of freedom
    end
    arguments (Output)
        quad_nodes      % (1, dof) node vector
        quad_weights    % (1, dof) weight vector
    end
    switch dof
        case 2
            quad_nodes = [-1, 1];
            quad_weights = [1, 1];
        case 3
            quad_nodes = [-1, 0, 1];
            quad_weights = [1, 4, 1]/3;    
        case 4
            quad_nodes = [-1, -1/sqrt(5), 1/sqrt(5), 1];
            quad_weights = [1, 5, 5, 1]/6;            
        otherwise
            error("method has not been implemented for dof = " + dof)
    end
end