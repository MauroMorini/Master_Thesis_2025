function dirichlet_vect = dirichletbcVector1D(nodes, elements, c_handle, g_handle, sigma)
    % this method assembles the rhs dirichlet vector of the SIP-DG method in 1d
    % (purely boundary conditions, no load).
    % this means all the boundary values of the additionally introduced terms
    % in the bilinear form
    %
    % Inputs:
    %       nodes:      (num_nodes, 1) node value matrix
    %       elements:   (num_el, dof) connectivity (element index) matrix 
    %       c_handle:   @(x) function handle (elliptic inner function)
    %       g_handle:   @(x) function handle of the b.c.
    %       sigma:      scalar penalty coefficient
    %
    % Output:   
    %       dirichlet_vect:     (num_nodes, 1) vector with entries at the endpoints

    % initializations 
    num_nodes = length(nodes);
    dof = size(elements, 2);
    B_penalty_max_entries = 2*dof;
    triplet_list_idx = zeros(B_penalty_max_entries, 1);
    triplet_list_entries = zeros(B_penalty_max_entries, 1);
    triplet_list_iterator = 1;
    [phi, dphi] = fem1d.getBasisFun(dof);
    F_ref = @(x, xn_loc, h_loc) (x - xn_loc)/h_loc;

    % lower boundary face contribution      ISSUE: outward normal and index hardcoded!!!
    el_loc = elements(1,:);
    xk = nodes(el_loc(1));
    xn_loc = nodes(el_loc(1));
    h = abs(xk - nodes(el_loc(end)));
    outward_normal = -1;

    for loc_node_idx=1:dof
        triplet_list_idx(triplet_list_iterator) = el_loc(loc_node_idx);
        triplet_list_entries(triplet_list_iterator) =   c_handle(xk)*sigma/h*g_handle(xk)*phi{loc_node_idx}(F_ref(xk,xn_loc,h))-...
                                                        c_handle(xk)*dphi{loc_node_idx}(F_ref(xk,xn_loc,h))/h*g_handle(xk)*outward_normal;
        triplet_list_iterator = triplet_list_iterator + 1;
    end

    % upper boundary face contribution      ISSUE: outward normal and index hardcoded!!!
    el_loc = elements(end,:);
    xk = nodes(el_loc(end));
    xn_loc = nodes(el_loc(1));
    h = abs(xk - nodes(el_loc(1)));
    outward_normal = 1;

    for loc_node_idx=1:dof
        triplet_list_idx(triplet_list_iterator) = el_loc(loc_node_idx);
        triplet_list_entries(triplet_list_iterator) =   c_handle(xk)*sigma/h*g_handle(xk)*phi{loc_node_idx}(F_ref(xk,xn_loc,h))-...
                                                        c_handle(xk)*dphi{loc_node_idx}(F_ref(xk,xn_loc,h))/h*g_handle(xk)*outward_normal;
        triplet_list_iterator = triplet_list_iterator + 1;
    end
    
    dirichlet_vect = sparse(triplet_list_idx, ones(size(triplet_list_idx)), triplet_list_entries, num_nodes, 1);
end
