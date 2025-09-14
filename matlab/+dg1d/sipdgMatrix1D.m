function B = sipdgMatrix1D(nodes, elements, c_handle, sigma)
    % collects all LHS matrices used for SIP-DG into one sparse 
    % matrix to solve the system (no b.c included)
    %
    % Inputs:
    %       nodes:      (num_nodes, 1) node value matrix
    %       elements:   (num_el, dof) connectivity (element index) matrix 
    %       c_handle:   @(x) function handle 
    %       sigma:      scalar penalty value 
    %
    % Outputs:
    %       B:          (num_nodes, num_nodes) sparse RHS matrix

    A = fem1d.stiffnessMatrix1D(nodes, elements, c_handle);
    B_flux_int = dg1d.interiorFluxMatrix1D(nodes, elements, c_handle);
    B_flux_bound = dg1d.boundaryFluxMatrix1D(nodes, elements, c_handle);
    B_penalty_int = dg1d.interiorPenaltyMatrix1D(nodes, elements, c_handle, sigma);
    B_penalty_bound = dg1d.boundaryPenaltyMatrix1D(nodes, elements, c_handle, sigma);

    B = A - B_flux_int - B_flux_bound + B_penalty_int + B_penalty_bound;
end