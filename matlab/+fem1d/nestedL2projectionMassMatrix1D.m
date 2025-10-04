function M = nestedL2projectionMassMatrix1D(nodes_coarse,elements_coarse,nodes_fine,elements_fine)
        % V1 = span{phi_j^(coarse)}, V2 = span{phi_j^(fine)} are the span of hat basis functions over nodes coarse, fine
        % calculates mass matrix Mij = (phi_i^(2), phi_j^(1))_L2
        % requires for the spaces to be nested, meaning nodes_coarse subset nodes_fine
        % requires nodes to be ordered 
    arguments (Input)
        nodes_coarse            % (num_nodes_coarse,1) node vector
        elements_coarse         % (num_el_coarse, dof) connectivity matrix
        nodes_fine              % (num_nodes_fine,1) node vector
        elements_fine           % (num_el_fine, dof) connectivity matrix
    end

    arguments (Output)
        M               % (num_nodes_fine, num_nodes_coarse) sparse mass matrix
    end

    % Initializations
    num_el_coarse = size(elements_coarse, 1); 
    num_el_fine = size(elements_fine, 1);
    num_nodes_coarse = length(nodes_coarse);
    num_nodes_fine = length(nodes_fine);
    dof_c = size(elements_coarse, 2);
    dof_f = size(elements_fine, 2);
    triplet_list_iterator = 1;
    maxdof = max(dof_f, dof_c);

    % preallocation
    elements_ratio = ceil(num_el_fine/num_el_coarse);
    M_max_entries = num_el_coarse*elements_ratio*dof_f*dof_c;
    triplet_list_rows = zeros(M_max_entries, 1);
    triplet_list_cols = zeros(M_max_entries, 1);
    triplet_list_entries = zeros(M_max_entries, 1);

    % check inputs 
    assert(num_el_fine >= num_el_coarse, "the second input space has to be finer than the first")

    % collect quadrature information
    [phi_val, ~, quad_weights] = common.getShapeFunctionValueMatrix(dof_f, maxdof + 1);
    [quad_nodes, ~] = common.getLobattoQuadrature(maxdof+1);
    barycentric_weights_coarse = common.calculateBarycentricWeights(nodes_coarse, elements_coarse);

    % categorize finer mesh into coarser one
    faces_coarse = [nodes_coarse(elements_coarse(:,1)); nodes_coarse(end)];
    fine_to_coarse_idx_map = discretize(nodes_fine(elements_fine(:,1)), faces_coarse);

    for k = 1:num_el_coarse
        % find fine elements in coarse one
        fine_el_idx = find(fine_to_coarse_idx_map == k);
        fine_el_loc = elements_fine(fine_el_idx, :);

        % find quadrature points in fine elements
        fine_nodes_loc = nodes_fine(fine_el_loc(:, 1));
        fine_nodes_loc = fine_nodes_loc(:);

        % calculate basis function values at coarser nodes 
        weights_loc = barycentric_weights_coarse(elements_coarse(k,:));
        coarse_nodes_loc = nodes_coarse(elements_coarse(k,:));
        Phi_loc = common.evaluateLagrangeBarycentric(fine_nodes_loc, weights_loc, coarse_nodes_loc);

        for s = 1:length(fine_el_idx)
            h = abs(nodes_fine(elements_fine(fine_el_idx(s), end))-nodes_fine(elements_fine(fine_el_idx(s), 1)));
            Phi_loc_s = Phi_loc(:, ((s-1)*dof_f+1):(s*dof_f)).*quad_weights;
            entry_loc = Phi_loc_s*phi_val.';
            row_idx_loc = repmat(elements_coarse(k,:).', 1, dof_f);
            cols_idx_loc = repmat(elements_fine(fine_el_idx(s), :), dof_c, 1);
            entry_length = dof_c*dof_f;
            triplet_list_entries(triplet_list_iterator:triplet_list_iterator+entry_length-1) = entry_loc(:)*(h/2);
            triplet_list_rows(triplet_list_iterator:triplet_list_iterator+entry_length-1) = row_idx_loc(:);
            triplet_list_cols(triplet_list_iterator:triplet_list_iterator+entry_length-1) = cols_idx_loc(:);
            triplet_list_iterator = triplet_list_iterator + entry_length;
        end
    end
    triplet_list_entries = triplet_list_entries(1:triplet_list_iterator-1);
    triplet_list_rows = triplet_list_rows(1:triplet_list_iterator-1);
    triplet_list_cols = triplet_list_cols(1:triplet_list_iterator-1);
    M = sparse(triplet_list_rows, triplet_list_cols, triplet_list_entries, num_nodes_coarse, num_nodes_fine);
end