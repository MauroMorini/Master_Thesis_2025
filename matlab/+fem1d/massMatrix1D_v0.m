function M = massMatrix1D_v0(x, t, c)
        % calculates global mass matrix for linear or quadratic 
        % FE in 1D M(i,j) = int_Omega phi_i*phi_j*c
        % Inputs : 
        % x :  (1, nP) pointvector 
        % t : nEx3 connectivity matrix with elements in rows
        % c : function handle 
        %
        % Outputs : 
        % M : global stiffness matrix nPxnP
        
        % decide if quadratic or linear FE
        DoF = size(t, 2);
        
        nP = length(x);
        nE = size(t, 1);
        M = sparse(nP, nP);
        
        % iterate over elements
        for i = 1:nE
            
            % element
            K = x(t(i,:));
        
            % get element matrix 
            Mloc = fem1d.massElementMatrix1D_v1(K, c);
            
            % Assembling
            M(t(i, :), t(i, :)) = M(t(i, :), t(i, :)) + Mloc;
        end
        end