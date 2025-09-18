function load_vector = loadVector1D_v0(nodes, elements, f_handle)
    % assembles load vector depending on degrees of freedom
    %
    % Inputs: 
    %       nodes:          (num_nodes, 1) node matrix
    %       elements:       (num_nodes, dof) connectivity matrix
    %       f_handle:       @(x) function handle of RHS
    % Outputs:
    %       load_vector:    (num_nodes, 1) load vector
    dof = size(elements, 2);
    switch dof
        case 2
            load_vector = fem1d_old_versions.loadVectorLinear1D_v0(nodes.', elements, f_handle);
        case 3      
            load_vector = fem1d_old_versions.loadVectorQuadratic1D_v0(nodes.', elements, f_handle);        
        otherwise
            error("load vector for dof = " + dof + " has not been implemented")
    end
end