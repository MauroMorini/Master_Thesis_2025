function A = stiffnessMatrix1D(x, T, c)
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
    A = sparse(n,n);
    
    % iterate over all elements
    for i = 1:nEl
        
        % element
        K = x(T(i,:));
    
        % get element matrix 
        AK = fem1d.stiffnessElementMatrix1D(K, c);
    
        % assembling of stiffness matrix
        A(T(i,:), T(i,:)) = A(T(i,:), T(i,:)) + AK;
    
    end
end