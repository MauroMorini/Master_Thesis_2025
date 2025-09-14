function rhs = sipdgDirichletLoadVector1D(nodes, elements, f_handle, c_handle, g_handle, sigma)
    % this method assembles the rhs dirichlet load vector of the SIP-DG method in 1d.
    % this means the complete rhs of the discrete form.
    %
    % Inputs:
    %       nodes:      (num_nodes, 1) node value matrix
    %       elements:   (num_el, dof) connectivity (element index) matrix 
    %       f_handle:   @(x) function handle of load (rhs of pde)
    %       c_handle:   @(x) function handle (elliptic inner function)
    %       g_handle:   @(x) function handle of the b.c.
    %       sigma:      scalar penalty coefficient
    %
    % Output:   
    %       rhs:     (num_nodes, 1) full rhs vector for sipdg

    load_vector = fem1d.loadVectorLinear1D(nodes, elements, f_handle);
    dirichlet_vector = dg1d.dirichletbcVector1D(nodes, elements, c_handle, g_handle, sigma);

    rhs = load_vector + dirichlet_vector;
end
