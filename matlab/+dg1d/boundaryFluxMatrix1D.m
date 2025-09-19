function B_flux = boundaryFluxMatrix1D(nodes, elements, c_vals)
    % this method assembles the symmetric flux part of the SIP-DG method in 1d
    % for boundary faces, no penalty terms are assembled here.
    arguments (Input)
        nodes               % (num_nodes, 1) node value matrix
        elements            % (num_el, dof) connectivity (element index) matrix 
        c_vals              % (num_nodes, 1) vector with values of c at nodes (c(nodes))  
    end
    arguments (Output)
        B_flux              % (num_nodes, num_nodes) sparse matrix 
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
    [phi_val, dphi_val, ~] = common.getShapeFunctionValueMatrix(dof);

    % lower boundary face contribution     
    lower_boundary_element_idx = 1;
    h = abs(nodes(elements(lower_boundary_element_idx,end))-nodes(elements(lower_boundary_element_idx,1)));
    B_loc_1 = c_vals(elements(lower_boundary_element_idx,1))*(2/h)*(-1)*(phi_val(:,1)*dphi_val(:,1).' + dphi_val(:, 1)*phi_val(:,1).');
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
    B_loc_2 = c_vals(elements(upper_boundary_element_idx,end))*(2/h)*(1)*(phi_val(:,end)*dphi_val(:,end).' + dphi_val(:, end)*phi_val(:,end).');
    B_loc_2_rows = repmat(elements(upper_boundary_element_idx,:).',1,dof);
    B_loc_2_cols = repmat(elements(upper_boundary_element_idx,:), dof, 1);

    % fill local matrices into triplet lists
    triplet_list_rows(triplet_list_iterator:(triplet_list_iterator + dof^2 - 1)) = B_loc_2_rows(:);
    triplet_list_cols(triplet_list_iterator:(triplet_list_iterator + dof^2 - 1)) = B_loc_2_cols(:);
    triplet_list_entries(triplet_list_iterator:(triplet_list_iterator + dof^2 - 1)) = B_loc_2(:);

    B_flux = sparse(triplet_list_rows, triplet_list_cols, triplet_list_entries, num_nodes, num_nodes);
end
