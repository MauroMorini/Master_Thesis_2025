function [quad_nodes, quad_weights] = getGaussLegendreQuadrature(dof)
    % depending on dof parameter calculates quadrature points and weights
    arguments (Input)
        dof             % scalar degrees of freedom
    end
    arguments (Output)
        quad_nodes      % (1, dof) node vector
        quad_weights    % (1, dof) weight vector
    end
    % check for unsuitable dofs
    maxDof = 20;
    if dof < 2 || dof > maxDof
        error("dof=" + dof+" is unsuitable as degrees of freedom")
    end
    % cache for multiple use
    persistent cache
    if isempty(cache)
        cache = cell(1, maxDof);
    end

    if ~isempty(cache{dof})
        quad_nodes = cache{dof}.quad_nodes;
        quad_weights = cache{dof}.quad_weights;
        return
    end

    % get gauss nodes
    quad_nodes = common.getGaussLegendreQuadratureNodes(dof);

    % quad functions
    [phi_cell, ~] = common.getLagrangeBasisFun(quad_nodes);

    % get quad weights
    quad_weights = zeros(1, dof);
    for i = 1:length(phi_cell)
        quad_weights(i) = integral(phi_cell{i}, -1, 1);
    end

    cache{dof} = struct("quad_nodes", quad_nodes, "quad_weights", quad_weights);
end