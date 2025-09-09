function M = massMatrix1D(x, t, c)
        % calculates global mass matrix for linear or quadratic 
        % FE in 1D M(i,j) = int_Omega phi_i*phi_j*c
        % Inputs : 
        % x :  (1, nP) pointvector 
        % t : nEx3 connectivity matrix with elements in rows
        % c : function handle 
        %
        % Outputs : 
        % M : global stiffness matrix nPxnP
        
        % number of elements 
        nEl = size(t, 1); 

        n = length(x);
        dof = size(t, 2);
        M_max_entries = nEl*dof^2;
        triplet_list_rows = zeros(M_max_entries, 1);
        triplet_list_cols = zeros(M_max_entries, 1);
        triplet_list_entries = zeros(M_max_entries, 1);
        triplet_list_iterator = 1;
        
        % iterate over elements
        for k = 1:nEl
            
            % element
            K = x(t(k,:));
        
            % get element matrix 
            Mloc = fem1d.massElementMatrix1D(K, c);
            
            % assembling of stiffness matrix
            for i = 1:dof
                for j = 1:dof
                    triplet_list_rows(triplet_list_iterator) = t(k,i);
                    triplet_list_cols(triplet_list_iterator) = t(k,j);
                    triplet_list_entries(triplet_list_iterator) = Mloc(i,j);
                    triplet_list_iterator = triplet_list_iterator + 1;
                end
            end
        end
        M = sparse(triplet_list_rows, triplet_list_cols, triplet_list_entries, n, n);
        end