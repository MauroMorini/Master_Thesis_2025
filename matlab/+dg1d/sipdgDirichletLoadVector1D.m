function rhs = sipdgDirichletLoadVector1D(nodes, elements, f_vals, c_vals, g_vals, sigma)
    % this method assembles the rhs dirichlet load vector of the SIP-DG method in 1d.
    % this means the complete rhs of the discrete form.
    arguments (Input)
        nodes           % (num_nodes, 1) node value matrix
        elements        % (num_el, dof) connectivity (element index) matrix 
        f_vals          % (num_nodes,1) function values of forcing term at nodes
        c_vals          % @(x) function handle (elliptic inner function)
        g_vals          % @(x) function handle of the b.c.
        sigma           % scalar penalty coefficient
    end
    arguments (Output)
        rhs             % (num_nodes, 1) full rhs vector for sipdg

    end

    load_vector = fem1d.loadVector1D(nodes, elements, f_vals);
    dirichlet_vector = dg1d.dirichletbcVector1D(nodes, elements, c_vals, g_vals, sigma);

    rhs = load_vector + dirichlet_vector;
end
