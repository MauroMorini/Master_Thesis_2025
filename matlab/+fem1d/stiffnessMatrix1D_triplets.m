function A = stiffnessMatrix1D_triplets(x, T, c)
    % calculate the nxn stiffness matrix using a connectivity matrix T and a
    % possibly non equidistant grid x of size n for 2 or 3 degrees of freedom
    % A = [(c*phi'_i, phi'_j)_L^2]i,j
    %
    % Inputs :
    % T : (nEl, DoF) connectivity matrix
    % x : (1, n) point vector
    % c : function handle 
    %
    % Outputs:
    % A : (n,n) sparse stiffness matrix
            
    % number of elements 
    nEl = size(T, 1); 
    
    n = length(x);
    dof = size(T, 2);
    A_max_entries = nEl*dof^2;
    triplet_list_rows = zeros(A_max_entries, 1);
    triplet_list_cols = zeros(A_max_entries, 1);
    triplet_list_entries = zeros(A_max_entries, 1);
    triplet_list_iterator = 1;

    
    % iterate over all elements
    for k = 1:nEl
        
        % element
        K = x(T(k,:));
    
        % get element matrix 
        AK = fem1d.stiffnessElementMatrix1D(K, c);
    
        % assembling of stiffness matrix
        for i = 1:dof
            for j = 1:dof
                triplet_list_rows(triplet_list_iterator) = T(k,i);
                triplet_list_cols(triplet_list_iterator) = T(k,j);
                triplet_list_entries(triplet_list_iterator) = AK(i,j);
                triplet_list_iterator = triplet_list_iterator + 1;
            end
        end
    end

    A = sparse(triplet_list_rows, triplet_list_cols, triplet_list_entries, n, n);
end