function B_penalty = interiorPenaltyMatrix1D(nodes, elements, c_handle, sigma)
    % this method assembles the symmetric penalty part of the SIP-DG method in 1d
    % for interior faces.
    %
    % Inputs:
    %       nodes:      (num_nodes, 1) node value matrix
    %       elements:   (num_el, dof) connectivity (element index) matrix 
    %       c_handle:   @(x) function handle 
    %       sigma:      scalar penalty value 
    %
    % Output:   
    %       B_penalty:     (num_nodes, num_nodes) sparse matrix 

    % initializations 
    num_faces = size(elements, 1)+1; 
    num_nodes = length(nodes);
    dof = size(elements, 2);
    B_penalty_max_entries = (num_faces-2)*dof^2*4;
    triplet_list_rows = zeros(B_penalty_max_entries, 1);
    triplet_list_cols = zeros(B_penalty_max_entries, 1);
    triplet_list_entries = zeros(B_penalty_max_entries, 1);
    triplet_list_iterator = 1;

    for k = 2:num_faces-1
        bordering_elements = [elements(k-1,:); elements(k,:)];
        bordering_element_nodes = nodes(bordering_elements);
        xk = bordering_element_nodes(1, 2);
        outward_normal = [1, -1];
        
        % comment for clarity:
        % for each fixed face we know only the two bordering elements can
        % have impact due to the local support of the basis functions. So
        % we add up all combinations of basis functions (with jump, avg)
        % evaluated at xk. Note that the basis function only have support
        % in their respective element, such that {phi_i^s(xk)} = 1/2
        % phi_i^s(xk), [phi_i^s(xk)] = outward_normal*phi_i^s(xk)
        for el_idx_1=1:2
            for el_idx_2=1:2
                % collect information of respective element 1,2
                xn_loc_1 = bordering_element_nodes(el_idx_1, 1);        % lower element face node
                h_loc_1 = abs(xn_loc_1-bordering_element_nodes(el_idx_1, end));
                phi_1 = {@(x) 1- (x-xn_loc_1)/h_loc_1, @(x) (x-xn_loc_1)/h_loc_1};

                xn_loc_2 = bordering_element_nodes(el_idx_2, 1);    
                h_loc_2 = abs(xn_loc_2-bordering_element_nodes(el_idx_2, end));
                phi_2 = {@(x) 1- (x-xn_loc_2)/h_loc_2, @(x) (x-xn_loc_2)/h_loc_2};

                h_loc_min = min(h_loc_1, h_loc_2);

                for loc_node_idx_1=1:dof
                    for loc_node_idx_2=1:dof
                        triplet_list_rows(triplet_list_iterator) = bordering_elements(el_idx_1, loc_node_idx_1);
                        triplet_list_cols(triplet_list_iterator) = bordering_elements(el_idx_2, loc_node_idx_2);
                        triplet_list_entries(triplet_list_iterator) =   sigma/h_loc_min*phi_1{loc_node_idx_1}(xk)*phi_2{loc_node_idx_2}(xk)...
                                                                        *outward_normal(loc_node_idx_1)*outward_normal(loc_node_idx_2)...
                                                                        *max(c_handle(xk-eps),c_handle(xk+eps));
                        triplet_list_iterator = triplet_list_iterator + 1;
                    end
                end
            end
        end
    end

    B_penalty = sparse(triplet_list_rows, triplet_list_cols, triplet_list_entries, num_nodes, num_nodes);
end
