function neumann_vect = neumannbcVector1D(nodes, elements, c_vals, g_vals)
    % this method assembles the rhs dirichlet vector of the SIP-DG method in 1d
    % (purely boundary conditions, no load).
    % this means all the boundary values of the additionally introduced terms
    % in the bilinear form
    %
    arguments (Input)
        nodes           % (num_nodes, 1) node value matrix
        elements        % (num_el, dof) connectivity (element index) matrix 
        c_vals          % (num_el, num_quad) (elliptic inner function)
        g_vals          % (num_nodes,1)  matrix of the b.c.
    end
    arguments (Output)
        neumann_vect  % (num_nodes, 1) boundary condition rhs vector
    end

    % initializations 
    num_nodes = length(nodes);
    dof = size(elements, 2);
    B_flux_max_entries = 2*dof;

    % preallocation
    triplet_list_rows = zeros(B_flux_max_entries, 1);
    triplet_list_entries = zeros(B_flux_max_entries, 1);

    % get basis function values
    [phi_val, dphi_val, ~] = common.getShapeFunctionValueMatrix(dof);

    % lower boundary face contribution     
    lower_boundary_element_idx = 1;
    h = abs(nodes(elements(lower_boundary_element_idx,end))-nodes(elements(lower_boundary_element_idx,1)));
    triplet_list_rows(1:dof) = elements(lower_boundary_element_idx,:).';
    triplet_list_entries(1:dof) = c_vals(lower_boundary_element_idx,1)*g_vals(elements(lower_boundary_element_idx,1))*(phi_val(:,1)); % note here we have no negative sign because g = n*du

    % upper boundary face contribution     
    upper_boundary_element_idx = size(elements,1);
    h = abs(nodes(elements(upper_boundary_element_idx,end))-nodes(elements(upper_boundary_element_idx,1)));
    triplet_list_rows(dof+1:end) = elements(upper_boundary_element_idx,:).';
    triplet_list_entries(dof+1:end) = c_vals(upper_boundary_element_idx,end)*g_vals(elements(upper_boundary_element_idx,end))*(phi_val(:,end)); 

    neumann_vect = sparse(triplet_list_rows, 1, triplet_list_entries, num_nodes, 1);
end
