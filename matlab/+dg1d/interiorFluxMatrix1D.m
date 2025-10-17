function B_flux = interiorFluxMatrix1D(nodes, elements, c_vals)
    % this method assembles the symmetric flux part of the SIP-DG method in 1d
    % for interior faces, no penalty terms are assembled here.

    arguments (Input)
        nodes               % (num_nodes, 1) node value matrix
        elements            % (num_el, dof) connectivity (element index) matrix 
        c_vals              % (num_el, num_quad) vector with values of c at nodes (c(nodes))  
    end
    arguments (Output)
        B_flux              % (num_nodes, num_nodes) sparse matrix 
    end

    % initializations 
    num_faces = size(elements, 1)+1; 
    num_nodes = length(nodes);
    dof = size(elements, 2);
    B_flux_max_entries = (num_faces-2)*dof^2*4;
    triplet_list_iterator = 1;

    % preallocation
    triplet_list_rows = zeros(B_flux_max_entries, 1);
    triplet_list_cols = zeros(B_flux_max_entries, 1);
    triplet_list_entries = zeros(B_flux_max_entries, 1);

    % build local matrices (still need to multiply c, (1/h) 
    [phi_val, dphi_val, ~] = common.getShapeFunctionValueMatrix(dof);
    B_loc_1 = 1*(phi_val(:,end)*dphi_val(:,end).' + dphi_val(:, end)*phi_val(:,end).');
    B_loc_21 = 1*phi_val(:,end)*dphi_val(:,1).';
    B_loc_22 = (-1)*dphi_val(:,end)*phi_val(:,1).';
    B_loc_3 = (-1)*(phi_val(:,1)*dphi_val(:,1).' + dphi_val(:, 1)*phi_val(:,1).');


    for k = 2:num_faces-1
        % local matrix
        h_loc_1 = abs(nodes(elements(k-1, 1))-nodes(elements(k-1, end)));
        h_loc_2 = abs(nodes(elements(k, 1))-nodes(elements(k, end)));
        B_loc = [c_vals(k-1, end)*B_loc_1/h_loc_1, (c_vals(k, 1)*B_loc_21/h_loc_2 + c_vals(k-1, end)*B_loc_22/h_loc_1);
                 (c_vals(k, 1)*B_loc_21/h_loc_2 + c_vals(k-1, end)*B_loc_22/h_loc_1).', c_vals(k, 1)*B_loc_3/h_loc_2];
        bordering_elements = [elements(k-1,:), elements(k,:)];
        B_loc_row_idx = repmat(bordering_elements.', 1, 2*dof);
        B_loc_col_idx = repmat(bordering_elements, 2*dof, 1);

        % fill local matrices into triplet lists
        triplet_list_rows(triplet_list_iterator:(triplet_list_iterator + 4*dof^2 - 1)) = B_loc_row_idx(:);
        triplet_list_cols(triplet_list_iterator:(triplet_list_iterator + 4*dof^2 - 1)) = B_loc_col_idx(:);
        triplet_list_entries(triplet_list_iterator:(triplet_list_iterator + 4*dof^2 - 1)) = B_loc(:);
        triplet_list_iterator = triplet_list_iterator + 4*dof^2;
    end

    B_flux = sparse(triplet_list_rows, triplet_list_cols, triplet_list_entries, num_nodes, num_nodes);
end
