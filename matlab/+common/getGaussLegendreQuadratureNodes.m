function quad_nodes = getGaussLegendreQuadratureNodes(dof)
    % depending on dof parameter yields quadrature points and weights, currently only yields 
    % uniformly distributed nodal basis
    arguments (Input)
        dof             % scalar degrees of freedom
    end
    arguments (Output)
        quad_nodes      % (1, dof) node vector
    end

    % initialize
    quad_nodes = zeros(1, dof);
    quad_nodes(1) = -1;
    quad_nodes(end) = 1;
    if dof == 2
        return
    end

    % get roots of legendre polynomial (inner nodes)
    syms x 
    legendre_poly = legendreP(dof-1, x);
    dlegendre_poly = diff(legendre_poly, 1);
    roots = double(vpasolve(dlegendre_poly == 0));

    quad_nodes(2:end-1) = roots; 
end