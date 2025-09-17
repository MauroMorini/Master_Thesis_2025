function B_penalty = boundaryPenaltyMatrix1D(nodes, elements, c_handle, sigma)
    % this method assembles the symmetric penalty part of the SIP-DG method in 1d
    % for boundary faces.
    %
    % Inputs:
    %       nodes:      (num_nodes, 1) node value matrix
    %       elements:   (num_el, dof) connectivity (element index) matrix 
    %       c_handle:   @(x) function handle 
    %
    % Output:   
    %       B_penalty:     (num_nodes, num_nodes) sparse matrix 

    % initializations 
    num_nodes = length(nodes);
    dof = size(elements, 2);
    B_penalty_max_entries = 2*dof^2;
    triplet_list_rows = zeros(B_penalty_max_entries, 1);
    triplet_list_cols = zeros(B_penalty_max_entries, 1);
    triplet_list_entries = zeros(B_penalty_max_entries, 1);
    triplet_list_iterator = 1;
    [phi, dphi] = fem1d.getBasisFun(dof);
    F_ref = @(x, xn_loc, h_loc) (x - xn_loc)/h_loc;

    % lower boundary face contribution      ISSUE: outward normal and index hardcoded!!!
    el_loc = elements(1,:);
    xk = nodes(el_loc(1));
    xn_loc = nodes(el_loc(1));
    h = abs(xk - nodes(el_loc(end)));

    for loc_node_idx_1=1:dof
        for loc_node_idx_2=1:dof
            triplet_list_rows(triplet_list_iterator) = el_loc(loc_node_idx_1);
            triplet_list_cols(triplet_list_iterator) = el_loc(loc_node_idx_2);
            triplet_list_entries(triplet_list_iterator) = c_handle(xk)*sigma/h*phi{loc_node_idx_1}(F_ref(xk,xn_loc,h))*phi{loc_node_idx_2}(F_ref(xk,xn_loc,h));
            triplet_list_iterator = triplet_list_iterator + 1;
        end
    end

    % upper boundary face contribution      ISSUE: outward normal and index hardcoded!!!
    el_loc = elements(end,:);
    xk = nodes(el_loc(end));
    xn_loc = nodes(el_loc(1));
    h = abs(xk - nodes(el_loc(1)));

    for loc_node_idx_1=1:dof
        for loc_node_idx_2=1:dof
            triplet_list_rows(triplet_list_iterator) = el_loc(loc_node_idx_1);
            triplet_list_cols(triplet_list_iterator) = el_loc(loc_node_idx_2);
            triplet_list_entries(triplet_list_iterator) = c_handle(xk)*sigma/h*phi{loc_node_idx_1}(F_ref(xk,xn_loc,h))*phi{loc_node_idx_2}(F_ref(xk,xn_loc,h));
            triplet_list_iterator = triplet_list_iterator + 1;
        end
    end
    
    B_penalty = sparse(triplet_list_rows, triplet_list_cols, triplet_list_entries, num_nodes, num_nodes);
end
