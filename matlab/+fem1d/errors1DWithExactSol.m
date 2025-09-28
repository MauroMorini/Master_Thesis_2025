function [l2_error, h1_error] = errors1DWithExactSol(nodes, elements, uh_vals, u_exact_vals, du_exact_vals)
    % calculates the error in L^2 and H^1 norm for FEM. For this to work u_exact_vals
    % has to contain the values of the exact solution at the exact quadrature points
    % (so far just at the nodes since quadrature and nodes coincide for Gauss-Lobatto)
    arguments (Input)
        nodes                   % (num_nodes, 1) node vector
        elements                % (num_nodes, num_el) connectivity matrix
        uh_vals                 % (num_nodes, 1) values (coefficients) of numerical sol 
        u_exact_vals            % (num_nodes, 1) values of exact sol AT quadrature points in elements
        du_exact_vals           % (num_nodes, 1) values of derivative of exact sol AT quadrature points in elements
    end
    arguments (Output)
        l2_error                % scalar value of L^2 error over domain
        h1_error                % scalar value of H^1 error over domain     
    end
    % initialization
    [num_el, dof] = size(elements);

    % preallocation
    l2_error = 0;
    h1_error = 0;

    % collect quadrature information
    [phi_val, dphi_val, quad_weights] = common.getShapeFunctionValueMatrix(dof);
    
    % iterate over elements
    for k = 1:num_el
    
        % element
        K = nodes(elements(k, :));
    
        % element length
        h = abs(K(1) - K(end));
        
        % values of u_h, du_h at quadrature points
        uh_loc_val = (uh_vals(elements(k,:)).'*phi_val);
        duh_loc_val = (uh_vals(elements(k,:)).'*dphi_val);

        l2_error = l2_error + (h/2)*(u_exact_vals(elements(k,:)).' - uh_loc_val).^2*quad_weights.';
        h1_error = h1_error + (h/2)*(du_exact_vals(elements(k,:)).' - (2/h)*duh_loc_val).^2*quad_weights.';
       
    end
    h1_error = sqrt(h1_error + l2_error);
    l2_error = sqrt(l2_error);
        
end