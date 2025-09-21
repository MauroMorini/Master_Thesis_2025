function B = sipdgMatrix1D(nodes, elements, c_vals, sigma)
    % collects all LHS matrices used for SIP-DG into one sparse 
    % matrix to solve the system (no b.c included)
    arguments (Input)
        nodes               % (num_nodes, 1) node value matrix
        elements            % (num_el, dof) connectivity (element index) matrix 
        c_vals              % (num_nodes, 1) vector with values of c at nodes (c(nodes))  
        sigma               % scalar penalty constant 
    end
    arguments (Output)
        B                   % (num_nodes, num_nodes) sparse matrix 
    end

    A = fem1d.stiffnessMatrix1D(nodes, elements, c_vals);
    B_flux_int = dg1d.interiorFluxMatrix1D(nodes, elements, c_vals);
    B_flux_bound = dg1d.boundaryFluxMatrix1D(nodes, elements, c_vals);
    B_penalty_int = dg1d.interiorPenaltyMatrix1D(nodes, elements, c_vals, sigma);
    B_penalty_bound = dg1d.boundaryPenaltyMatrix1D(nodes, elements, c_vals, sigma);

    B = A - B_flux_int - B_flux_bound + B_penalty_int + B_penalty_bound;
end