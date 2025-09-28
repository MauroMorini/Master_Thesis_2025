function [l2_error, h1_error] = errors1DBetweenFemSol(nodes, elements, uh1_vals, uh2_vals)
    % calculates the  relative error in L^2 and H^1 norm for FEM between two fem solutions on the same mesh
    % 
    arguments (Input)
        nodes                   % (num_nodes, 1) node vector
        elements                % (num_nodes, num_el) connectivity matrix
        uh1_vals                % (num_nodes, 1) values (coefficients) of numerical sol 1
        uh2_vals                % (num_nodes, 1) values (coefficients) of numerical sol 1
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
        uh1_loc_val = (uh1_vals(elements(k,:)).'*phi_val);
        duh1_loc_val = (uh1_vals(elements(k,:)).'*dphi_val);

        uh2_loc_val = (uh2_vals(elements(k,:)).'*phi_val);
        duh2_loc_val = (uh2_vals(elements(k,:)).'*dphi_val);

        l2_error = l2_error + (h/2)*(uh1_loc_val - uh2_loc_val).^2*quad_weights.';
        h1_error = h1_error + (2/h)*(duh1_loc_val - duh2_loc_val).^2*quad_weights.';
       
    end
    h1_error = sqrt(h1_error + l2_error);
    l2_error = sqrt(l2_error);
        
end