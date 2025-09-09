function B_flux = interiorFluxMatrix1D(nodes, elements)

    % initializations 
    num_faces = size(elements, 1)+1; 
    num_nodes = length(nodes);
    dof = size(elements, 2);
    B_flux_max_entries = (num_faces-2)*dof^2*4;
    triplet_list_rows = zeros(B_flux_max_entries, 1);
    triplet_list_cols = zeros(B_flux_max_entries, 1);
    triplet_list_entries = zeros(B_flux_max_entries, 1);
    triplet_list_iterator = 1;

    for k = 2:num_faces-1
        bordering_elements = [elements(k-1,:); elements(k,:)];
        bordering_element_nodes = nodes(bordering_elements);
        xk = bordering_element_nodes(1, 2);
        outward_normal = [1, -1];

        for el_idx_1=1:2
            for el_idx_2=1:2
                xn_loc_1 = bordering_element_nodes(el_idx_1, 1);    
                h_loc_1 = abs(xn_loc_1-bordering_element_nodes(el_idx_1, end));
                phi_1 = {@(x) 1- (x-xn_loc_1)/h_loc_1, @(x) (x-xn_loc_1)/h_loc_1};
                dphi_1 = {@(x) -ones(size(x))/h_loc_1, @(x) ones(size(x))/h_loc_1};

                xn_loc_2 = bordering_element_nodes(el_idx_2, 1);    
                h_loc_2 = abs(xn_loc_2-bordering_element_nodes(el_idx_2, end));
                phi_2 = {@(x) 1- (x-xn_loc_2)/h_loc_2, @(x) (x-xn_loc_2)/h_loc_2};
                dphi_2 = {@(x) -ones(size(x))/h_loc_2, @(x) ones(size(x))/h_loc_2};

                for loc_node_idx_1=1:dof
                    for loc_node_idx_2=1:dof
                        triplet_list_rows(triplet_list_iterator) = bordering_elements(el_idx_1, loc_node_idx_1);
                        triplet_list_cols(triplet_list_iterator) = bordering_elements(el_idx_2, loc_node_idx_2);
                        triplet_list_entries(triplet_list_iterator) =   dphi_2{loc_node_idx_2}(xk)*outward_normal(el_idx_2)*phi_1{loc_node_idx_1}(xk)*(1/2)+...
                                                                        dphi_2{loc_node_idx_2}(xk)*(1/2)*phi_1{loc_node_idx_1}(xk)*outward_normal(el_idx_1);
                        triplet_list_iterator = triplet_list_iterator + 1;
                    end
                end
            end
        end
    end

    B_flux = sparse(triplet_list_rows, triplet_list_cols, triplet_list_entries, num_nodes, num_nodes);
end
