function energy_error = energyNormError1D(nodes, elements, uh_vals, c_vals, sigma, u_exact_vals, du_exact_vals)
    % calculates error of numerical solution in energy norm 
    arguments (Input)
        nodes                   % (num_nodes, 1) node vector
        elements                % (num_el, dof) connectivity matrix
        uh_vals                 % (num_nodes, 1) values of numerical sol (coefficient vector of FEM sol)
        c_vals                  % (num_nodes, 1) values of coefficient matrix 
        sigma double            % penalization parameter
        u_exact_vals            % (num_nodes, 1) values of exact solution at nodes (interpolant)
        du_exact_vals           % (num_nodes, 1) values of exact solution at nodes 
    end
    arguments (Output)
        energy_error double     % scalar energy error 
    end

    % initializations
    der_error = 0;
    penalty_error = 0;
    [num_el, dof] = size(elements);

    % collect quadrature information
    [phi_val, dphi_val, quad_weights] = common.getShapeFunctionValueMatrix(dof);

    % calculate derivative error part 
    for k = 1:num_el
    
        % element
        K = nodes(elements(k, :));
    
        % element length
        h = abs(K(1) - K(end));
        
        % values of u_h at quadrature points
        duh_loc_val = (uh_vals(elements(k,:)).*sqrt(c_vals(elements(k,:)))).'*dphi_val;

        % apply quadrature
        der_error = der_error + (h/2)*(du_exact_vals(elements(k,:)).' - (2/h)*duh_loc_val).^2*quad_weights.';
    end
    
    % calculate penalty error part first for interior faces
    meshsizes = abs(nodes(elements(:, 1)) - nodes(elements(:, end)));
    for n = 2:num_el
        h_min = min(meshsizes(n-1), meshsizes(n));
        c_max = max([c_vals(elements(n-1,:)); c_vals(elements(n, :))]);
        uh_loc_lower = uh_vals(elements(n-1, :)).'*phi_val(:, end);
        uh_loc_upper = uh_vals(elements(n, :)).'*phi_val(:, 1);
        u_exact_loc_lower = u_exact_vals(elements(n-1, end));
        u_exact_loc_upper = u_exact_vals(elements(n, 1));
        penalty_error = penalty_error + sigma*c_max/h_min*( (u_exact_loc_lower-u_exact_loc_upper) - (uh_loc_lower - uh_loc_upper) )^2;
    end

    % calculate penalty error part for boundary faces
    boundary_el_idx = [1, num_el];

    % lower boundary face
    h_min = meshsizes(boundary_el_idx(1));
    c_max = max(c_vals(boundary_el_idx(1), :));
    uh_loc = -uh_vals(elements(boundary_el_idx(1), :)).'*phi_val(:, 1);
    u_exact_loc = -u_exact_vals(elements(boundary_el_idx(1), 1));
    penalty_error = penalty_error + sigma*c_max/h_min * ( u_exact_loc - uh_loc)^2;

    % upper boundary face
    h_min = meshsizes(boundary_el_idx(2));
    c_max = max(c_vals(boundary_el_idx(2), :));
    uh_loc = uh_vals(elements(boundary_el_idx(2), :)).'*phi_val(:, end);
    u_exact_loc = u_exact_vals(elements(boundary_el_idx(2), end));
    penalty_error = penalty_error + sigma*c_max/h_min * ( u_exact_loc - uh_loc)^2;    

    energy_error = sqrt(der_error + penalty_error);
end
