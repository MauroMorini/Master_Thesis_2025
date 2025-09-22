function [phi_val, dphi_val, quad_weights] = getShapeFunctionValueMatrix(dof)
    % calculates matrix which contains the values of the shape functions
    % and their derivatives at the quadrature nodes
    arguments (Input)
        dof             % scalar degrees of freedom
    end
    arguments (Output)
        phi_val         % (dof, dof) value matrix for phi
        dphi_val        % (dof, dof) value matrix for dphi
        quad_weights    % (1, dof) weight vector
    end
    % cache values for speed
    maxDof = 20;
    persistent cache

    if isempty(cache)
        cache = cell(1,maxDof);
    end

    if ~isempty(cache{dof})
        phi_val = cache{dof}.phi_val;
        dphi_val = cache{dof}.dphi_val;
        quad_weights = cache{dof}.quad_weights;
        return
    end
    if dof == 2
        phi_val = [1, 0;
                    0, 1];
        dphi_val = (1/2)*[-1, -1;
                            1, 1];
        quad_weights = [1, 1];
    else
        phi_val = zeros(dof,dof);
        dphi_val = zeros(dof, dof);
        [quad_nodes, quad_weights] = common.getLobattoQuadrature(dof);
        [phi_cell, dphi_cell] = common.getLagrangeBasisFun(quad_nodes);
        for i = 1:dof
            phi_val(i,:) = phi_cell{i}(quad_nodes);
            dphi_val(i,:) = dphi_cell{i}(quad_nodes);
        end
    end
    cache{dof} = struct("phi_val", phi_val, "dphi_val", dphi_val,"quad_weights",quad_weights);
end