function B_flux = interiorFluxMatrix1D(nodes, elements)

    % initializations 
    nEl = size(elements, 1); 
    n = length(nodes);
    dof = size(elements, 2);
    A_max_entries = nEl*dof^2;
    triplet_list_rows = zeros(A_max_entries, 1);
    triplet_list_cols = zeros(A_max_entries, 1);
    triplet_list_entries = zeros(A_max_entries, 1);
    triplet_list_iterator = 1;


end
