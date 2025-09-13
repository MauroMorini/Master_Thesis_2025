function B_flux = boundaryFluxMatrix1D(nodes, elements, c_handle)
    % this method assembles the symmetric flux part of the SIP-DG method in 1d
    % for boundary faces, no penalty terms are assembled here.
    %
    % Inputs:
    %       nodes:      (num_nodes, 1) node value matrix
    %       elements:   (num_el, dof) connectivity (element index) matrix 
    %       c_handle:   @(x) function handle 
    %
    % Output:   
    %       B_flux:     (num_nodes, num_nodes) sparse matrix 

    % initializations 
    num_nodes = length(nodes);
    dof = size(elements, 2);
    B_flux_max_entries = 2*dof^2;
    triplet_list_rows = zeros(B_flux_max_entries, 1);
    triplet_list_cols = zeros(B_flux_max_entries, 1);
    triplet_list_entries = zeros(B_flux_max_entries, 1);
    triplet_list_iterator = 1;

    % lower boundary face contribution      ISSUE: outward normal and index hardcoded!!!
    el_loc = elements(1,:);
    xk = nodes(el_loc(1));
    xn_loc = nodes(el_loc(1));
    outward_normal = -1;
    h = abs(xk - nodes(el_loc(end)));
    phi = {@(x) 1- (x-xn_loc)/h, @(x) (x-xn_loc)/h};
    dphi = {@(x) -ones(size(x))/h, @(x) ones(size(x))/h};

    for loc_node_idx_1=1:dof
        for loc_node_idx_2=1:dof
            triplet_list_rows(triplet_list_iterator) = el_loc(loc_node_idx_1);
            triplet_list_cols(triplet_list_iterator) = el_loc(loc_node_idx_2);
            triplet_list_entries(triplet_list_iterator) =   c_handle(xk)*dphi{loc_node_idx_2}(xk)*outward_normal*phi{loc_node_idx_1}(xk)*(1/2)+...
                                                            c_handle(xk)*dphi{loc_node_idx_1}(xk)*outward_normal*phi{loc_node_idx_2}(xk)*(1/2);
            triplet_list_iterator = triplet_list_iterator + 1;
        end
    end

    % upper boundary face contribution      ISSUE: outward normal and index hardcoded!!!
    el_loc = elements(end,:);
    xk = nodes(el_loc(end));
    xn_loc = nodes(el_loc(1));
    outward_normal = -1;
    h = abs(xk - nodes(el_loc(1)));
    phi = {@(x) 1- (x-xn_loc)/h, @(x) (x-xn_loc)/h};
    dphi = {@(x) -ones(size(x))/h, @(x) ones(size(x))/h};

    for loc_node_idx_1=1:dof
        for loc_node_idx_2=1:dof
            triplet_list_rows(triplet_list_iterator) = el_loc(loc_node_idx_1);
            triplet_list_cols(triplet_list_iterator) = el_loc(loc_node_idx_2);
            triplet_list_entries(triplet_list_iterator) =   c_handle(xk)*dphi{loc_node_idx_2}(xk)*outward_normal*phi{loc_node_idx_1}(xk)*(1/2)+...
                                                            c_handle(xk)*dphi{loc_node_idx_1}(xk)*outward_normal*phi{loc_node_idx_2}(xk)*(1/2);
            triplet_list_iterator = triplet_list_iterator + 1;
        end
    end



    
    B_flux = sparse(triplet_list_rows, triplet_list_cols, triplet_list_entries, num_nodes, num_nodes);
end
