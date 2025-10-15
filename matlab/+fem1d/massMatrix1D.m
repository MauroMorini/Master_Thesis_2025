function M = massMatrix1D(nodes, elements, c_vals)
    % calculates global mass matrix for linear or quadratic 
    % FE in 1D M(i,j) = int_Omega phi_i*phi_j*c
    arguments (Input)
        nodes                               % (num_nodes,1) node vector
        elements                            % (num_el, dof) connectivity matrix
        c_vals double                       % (num_el, num_quad) values of coefficient function at quadrature nodes    
    end

    arguments (Output)
        M               % (num_nodes, num_nodes) sparse mass matrix
    end
    
    % initializations
    nEl = size(elements, 1); 
    num_nodes = length(nodes);
    dof = size(elements, 2);
    triplet_list_iterator = 1;

    num_quad = size(c_vals, 2);

    % preallocation
    M_max_entries = nEl*dof^2;
    triplet_list_rows = zeros(M_max_entries, 1);
    triplet_list_cols = zeros(M_max_entries, 1);
    triplet_list_entries = zeros(M_max_entries, 1);
    
    % collect quadrature information
    [phi_val, ~, quad_weights] = common.getShapeFunctionValueMatrix(dof, num_quad);

    % iterate over all elements
    for k = 1:nEl
        
        % element
        K = nodes(elements(k,:));
        
        % element length
        h = abs(K(end) - K(1));
    
        % calculate local mass matrix
        M_loc = zeros(dof);
        for p = 1:dof
            for q = 1:p
                M_loc(p,q) = (phi_val(p,:).*phi_val(q,:).*c_vals(k,:))*quad_weights.';
                M_loc(q,p) = M_loc(p,q);
            end
        end
        M_loc = M_loc*(h/2);
    
        % assembling of mass matrix    TODO: vectorize
        for i = 1:dof
            for j = 1:dof
                triplet_list_rows(triplet_list_iterator) = elements(k,i);
                triplet_list_cols(triplet_list_iterator) = elements(k,j);
                triplet_list_entries(triplet_list_iterator) = M_loc(i,j);
                triplet_list_iterator = triplet_list_iterator + 1;
            end
        end
    end

    M = sparse(triplet_list_rows, triplet_list_cols, triplet_list_entries, num_nodes, num_nodes);
end