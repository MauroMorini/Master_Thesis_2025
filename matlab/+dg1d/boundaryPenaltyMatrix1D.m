function B_penalty = boundaryPenaltyMatrix1D(nodes, elements, c_vals, sigma)
    % this method assembles the symmetric penalty part of the SIP-DG method in 1d
    % for boundary faces.
    arguments (Input)
        nodes               % (num_nodes, 1) node value matrix
        elements            % (num_el, dof) connectivity (element index) matrix 
        c_vals              % (num_el, num_quad) vector with values of c at nodes (c(nodes))  
        sigma               % scalar penalty constant 
    end
    arguments (Output)
        B_penalty              % (num_nodes, num_nodes) sparse matrix 
    end

% initializations 
    num_nodes = length(nodes);
    dof = size(elements, 2);
    B_flux_max_entries = 2*dof^2;
    triplet_list_iterator = 1;

    % preallocation
    triplet_list_rows = zeros(B_flux_max_entries, 1);
    triplet_list_cols = zeros(B_flux_max_entries, 1);
    triplet_list_entries = zeros(B_flux_max_entries, 1);

    % get basis function values
    [phi_val, ~, ~] = common.getShapeFunctionValueMatrix(dof);

    % lower boundary face contribution     
    lower_boundary_element_idx = 1;
    h = abs(nodes(elements(lower_boundary_element_idx,end))-nodes(elements(lower_boundary_element_idx,1)));
    B_loc_1 = c_vals(lower_boundary_element_idx,1)*sigma/h*(1)*phi_val(:,1)*phi_val(:,1).';
    B_loc_1_rows = repmat(elements(lower_boundary_element_idx,:).',1,dof);
    B_loc_1_cols = repmat(elements(lower_boundary_element_idx,:), dof, 1);

    % fill local matrices into triplet lists
    triplet_list_rows(triplet_list_iterator:(triplet_list_iterator + dof^2 - 1)) = B_loc_1_rows(:);
    triplet_list_cols(triplet_list_iterator:(triplet_list_iterator + dof^2 - 1)) = B_loc_1_cols(:);
    triplet_list_entries(triplet_list_iterator:(triplet_list_iterator + dof^2 - 1)) = B_loc_1(:);
    triplet_list_iterator = triplet_list_iterator + dof^2;

    % upper boundary face contribution     
    upper_boundary_element_idx = size(elements,1);
    h = abs(nodes(elements(upper_boundary_element_idx,end))-nodes(elements(upper_boundary_element_idx,1)));
    B_loc_2 = c_vals(upper_boundary_element_idx,end)*sigma/h*(1)*phi_val(:,end)*phi_val(:,end).';
    B_loc_2_rows = repmat(elements(upper_boundary_element_idx,:).',1,dof);
    B_loc_2_cols = repmat(elements(upper_boundary_element_idx,:), dof, 1);

    % fill local matrices into triplet lists
    triplet_list_rows(triplet_list_iterator:(triplet_list_iterator + dof^2 - 1)) = B_loc_2_rows(:);
    triplet_list_cols(triplet_list_iterator:(triplet_list_iterator + dof^2 - 1)) = B_loc_2_cols(:);
    triplet_list_entries(triplet_list_iterator:(triplet_list_iterator + dof^2 - 1)) = B_loc_2(:);

    B_penalty = sparse(triplet_list_rows, triplet_list_cols, triplet_list_entries, num_nodes, num_nodes);
end
