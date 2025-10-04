function [phi_val, dphi_val, quad_weights] = getShapeFunctionValueMatrix(dof, num_quad_nodes)
    % calculates matrix which contains the values of the shape functions
    % and their derivatives at the quadrature nodes
    arguments (Input)
        dof (1,1) double                    % scalar degrees of freedom
        num_quad_nodes (1,1) double = 0     % number of quadrature nodes (default coincides with dof)
    end
    arguments (Output)
        phi_val         % (dof, dof) value matrix for phi
        dphi_val        % (dof, dof) value matrix for dphi
        quad_weights    % (1, dof) weight vector
    end

    if num_quad_nodes == 0
        num_quad_nodes = dof;
    end

    % cache values for speed
    maxDof = 20;
    max_num_quad_nodes = maxDof + 1;
    persistent cache

    if isempty(cache)
        cache = cell(maxDof, max_num_quad_nodes);
    end

    if ~isempty(cache{dof, num_quad_nodes})
        phi_val = cache{dof, num_quad_nodes}.phi_val;
        dphi_val = cache{dof, num_quad_nodes}.dphi_val;
        quad_weights = cache{dof, num_quad_nodes}.quad_weights;
        return
    end
    if dof == 2
        phi_val = [1, 0;
                    0, 1];
        dphi_val = (1/2)*[-1, -1;
                            1, 1];
        quad_weights = [1, 1];
    else
        phi_val = zeros(dof,num_quad_nodes);
        dphi_val = zeros(dof, num_quad_nodes);
        [quad_nodes, quad_weights] = common.getLobattoQuadrature(num_quad_nodes);
        [basis_nodes, ~] = common.getLobattoQuadrature(num_quad_nodes);
        [phi_cell, dphi_cell] = common.getLagrangeBasisFun(basis_nodes);
        for i = 1:dof
            phi_val(i,:) = phi_cell{i}(quad_nodes);
            dphi_val(i,:) = dphi_cell{i}(quad_nodes);
        end
    end
    cache{dof, num_quad_nodes} = struct("phi_val", phi_val, "dphi_val", dphi_val,"quad_weights",quad_weights);
end