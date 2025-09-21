function dirichlet_vect = dirichletbcVector1D(nodes, elements, c_vals, g_vals, sigma)
    % this method assembles the rhs dirichlet vector of the SIP-DG method in 1d
    % (purely boundary conditions, no load).
    % this means all the boundary values of the additionally introduced terms
    % in the bilinear form
    %
    arguments (Input)
        nodes           % (num_nodes, 1) node value matrix
        elements        % (num_el, dof) connectivity (element index) matrix 
        c_vals          % @(x) function handle (elliptic inner function)
        g_vals          % @(x) function handle of the b.c.
        sigma           % scalar penalty coefficient
    end
    arguments (Output)
        dirichlet_vect  % (num_nodes, 1) boundary condition rhs vector
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
    triplet_list_entries(1:dof) = c_vals(elements(lower_boundary_element_idx,1))/h*g_vals(elements(lower_boundary_element_idx,1))*(sigma*phi_val(:,1) + dphi_val(:,1)*(2));

    % upper boundary face contribution     
    upper_boundary_element_idx = size(elements,1);
    h = abs(nodes(elements(upper_boundary_element_idx,end))-nodes(elements(upper_boundary_element_idx,1)));
    triplet_list_rows(dof+1:end) = elements(upper_boundary_element_idx,:).';
    triplet_list_entries(dof+1:end) = c_vals(elements(upper_boundary_element_idx,end))/h*g_vals(elements(upper_boundary_element_idx,end))*(sigma*phi_val(:,end) - dphi_val(:,end)*(2));

    dirichlet_vect = sparse(triplet_list_rows, 1, triplet_list_entries, num_nodes, 1);
end
