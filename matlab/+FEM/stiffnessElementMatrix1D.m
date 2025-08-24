function AK = stiffnessElementMatrix1D(K, c)
    % calculates element matrix for a given element K on the grid, a function
    % handle c for 3 degrees of freedom (i.e. r=2). Code can be adapted for
    % higher grade polynomials
    %
    % Inputs:
    % K : (1,Dof) element point vector
    % 
    % Output:
    % AK : (Dof,Dof) local element matrix
    
    % Initializations
    Dof = size(K, 2);
    if Dof == 1
        K = K';
        Dof = size(K, 2);
    end
    AK = zeros(Dof);
    
    % initialize variables depending on Dof 
    switch Dof
        case 2
            % first derivative shape functions 
            N = {@(xi) -(1/2)*ones(size(xi)), @(xi) (1/2)*ones(size(xi))};
    
            % quadrature 
            y = [-1, 1];   
            w = [1, 1];     
            B = [N{1}(y); N{2}(y)];
        case 3
            % first derivative shape functions 
            N = {@(xi)1/2*(2*xi-1), @(xi) -2*xi, @(xi) 1/2*(1+2*xi)};
    
            % quadrature 
            y = [-1, 0, 1];   
            w = [1, 4, 1]/3;
            B = [N{1}(y); N{2}(y); N{3}(y)];
        otherwise
            error("local stiffness matrix has not been implemented for " + Dof + " degrees of freedom")
    end
    
    % element length
    h = abs(K(end) - K(1));
    
    % calculate functional values for K of c    
    cVal = c(K);
    
    % calculate integral for (A_K)_pq 
    for p = 1:Dof
        for q = 1:p
            AK(p,q) = (B(p,:).*B(q,:).*cVal)*w.';
            AK(q,p) = AK(p,q);
        end
    end
    AK = AK*(2/h);
end