function [uh, system_matrix] = sip_1d_elliptic_solver(Mesh, boundary_cond, f_vals, c_vals, sigma)
    % general symmetric interior penalty solver for elliptic problems
    % of the form -(cu')' = f 
    %
    arguments (Input)
        Mesh                            % MeshIntervalDG1d object containing mesh information
        boundary_cond struct            % a struct with the following members: values (1,2) double, lower_boundary_type string, upper_boundary_type string
        f_vals                          % (num_el, num_quad) vector with values of load at quadrature nodes
        c_vals                          % (num_el, num_quad) vector with values of c at quadrature nodes (can coincide with nodes, but doesn't have to)
        sigma double                    % scalar penalty constant 
    end
    arguments (Output)
        uh                              % (num_nodes, 1) numerical solution
        system_matrix
    end

    % extract mesh information
    [nodes, ~, elements] = Mesh.getPet();
    [lower_boundary_element_idx, upper_boundary_element_idx] = Mesh.getBoundaryElementIdx();

    % extract boundary condition
    g_vals = boundary_cond.values;

    % collect system matrices
    A = fem1d.stiffnessMatrix1D(nodes, elements, c_vals);
    M = fem1d.massMatrix1D(nodes, elements, ones(size(c_vals)));
    B_flux_int = dg1d.interiorFluxMatrix1D(nodes, elements, c_vals);
    B_flux_bound = dg1d.boundaryFluxMatrix1D(nodes, elements, c_vals);
    B_penalty_int = dg1d.interiorPenaltyMatrix1D(nodes, elements, c_vals, sigma);
    B_penalty_bound = dg1d.boundaryPenaltyMatrix1D(nodes, elements, c_vals, sigma);

    % collect system  vectors
    load_vector = fem1d.loadVector1D(nodes, elements, f_vals);
    dirichlet_vector = dg1d.dirichletbcVector1D(nodes, elements, c_vals, g_vals, sigma);
    neumann_vector = dg1d.neumannbcVector1D(nodes, elements, c_vals, g_vals);

    % implement boundary conditions
    dirichlet_bc_switch = [double(boundary_cond.lower_boundary_type == "dirichlet"), double(boundary_cond.upper_boundary_type == "dirichlet")];
    neumann_bc_switch = [double(boundary_cond.lower_boundary_type == "neumann"), double(boundary_cond.upper_boundary_type == "neumann")];

    B_flux_bound(elements(lower_boundary_element_idx,:), elements(lower_boundary_element_idx,:)) = B_flux_bound(elements(lower_boundary_element_idx,:), elements(lower_boundary_element_idx,:))*dirichlet_bc_switch(1);
    B_flux_bound(elements(upper_boundary_element_idx,:), elements(upper_boundary_element_idx,:)) = B_flux_bound(elements(upper_boundary_element_idx,:), elements(upper_boundary_element_idx,:))*dirichlet_bc_switch(2);
    B_penalty_bound(elements(lower_boundary_element_idx,:), elements(lower_boundary_element_idx,:)) = B_penalty_bound(elements(lower_boundary_element_idx,:), elements(lower_boundary_element_idx,:))*dirichlet_bc_switch(1);
    B_penalty_bound(elements(upper_boundary_element_idx,:), elements(upper_boundary_element_idx,:)) = B_penalty_bound(elements(upper_boundary_element_idx,:), elements(upper_boundary_element_idx,:))*dirichlet_bc_switch(2);
    dirichlet_vector(elements(lower_boundary_element_idx,:)) = dirichlet_vector(elements(lower_boundary_element_idx,:))*dirichlet_bc_switch(1);
    dirichlet_vector(elements(upper_boundary_element_idx,:)) = dirichlet_vector(elements(upper_boundary_element_idx,:))*dirichlet_bc_switch(2);

    neumann_vector(elements(lower_boundary_element_idx,:)) = neumann_vector(elements(lower_boundary_element_idx,:))*neumann_bc_switch(1);
    neumann_vector(elements(upper_boundary_element_idx,:)) = neumann_vector(elements(upper_boundary_element_idx,:))*neumann_bc_switch(2);

    system_matrix = A - B_flux_bound - B_flux_int + B_penalty_int + B_penalty_bound + M;
    system_vector = load_vector + dirichlet_vector + neumann_vector;

    % rescale system 
    scaling_factor = max(diag(system_matrix));
    system_matrix = system_matrix/scaling_factor;
    system_vector = system_vector/scaling_factor;

    % solve system using pcg
    tol = Mesh.h_max^Mesh.dof;
    maxit = numel(system_matrix);
    try
        L = ichol(system_matrix, struct('michol', 'on'));
        [uh, failed_to_converge_flag] = pcg(system_matrix,system_vector,tol,maxit,L,L');
    catch exception
        warning("ichol has failed minres is used")
        [L,U] = ilu(system_matrix);
        [uh, failed_to_converge_flag] = minres(system_matrix, system_vector, tol, maxit, L, U);
    end
    
    if failed_to_converge_flag
        warning("iterative solver has failed to converge normal mldivide is used")
        uh = system_matrix\system_vector;
    end
    % uh = system_matrix\system_vector;
end