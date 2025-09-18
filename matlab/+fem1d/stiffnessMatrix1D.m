function A = stiffnessMatrix1D(nodes, elements, c_vals)
    % calculate the nxn stiffness matrix 
    arguments (Input)
        nodes           % (num_nodes,1) node vector
        elements        % (num_el, dof) connectivity matrix
        c_vals          % (num_el, dof) values of coefficient function in elements
    end

    arguments (Output)
        A               % (num_nodes, num_nodes) sparse stiffness matrix
    end
            
    % initializations
    nEl = size(elements, 1); 
    num_nodes = length(nodes);
    dof = size(elements, 2);
    triplet_list_iterator = 1;

    % preallocation
    A_max_entries = nEl*dof^2;
    triplet_list_rows = zeros(A_max_entries, 1);
    triplet_list_cols = zeros(A_max_entries, 1);
    triplet_list_entries = zeros(A_max_entries, 1);
    
    % collect quadrature information
    [~, dphi_val, quad_weights] = common.getShapeFunctionValueMatrix(dof);

    % iterate over all elements
    for k = 1:nEl
        
        % element
        K = nodes(elements(k,:));
        
        % element length
        h = abs(K(end) - K(1));

        for i = 1:dof
            for j = 1:dof
                triplet_list_rows(triplet_list_iterator) = elements(k,i);
                triplet_list_cols(triplet_list_iterator) = elements(k,j);
                triplet_list_entries(triplet_list_iterator) = (dphi_val(i,:).*dphi_val(j,:).*c_vals(k,:))*quad_weights.'*(2/h);
                triplet_list_iterator = triplet_list_iterator + 1;
            end
        end
    
        % % calculate local stiffness matrix
        % AK = zeros(dof);
        % for p = 1:dof
        %     for q = 1:p
        %         AK(p,q) = (dphi_val(p,:).*dphi_val(q,:).*c_vals(k,:))*quad_weights.';
        %         AK(q,p) = AK(p,q);
        %     end
        % end
        % AK = AK*(2/h);
    
        % % assembling of stiffness matrix    TODO: vectorize
        % for i = 1:dof
        %     for j = 1:dof
        %         triplet_list_rows(triplet_list_iterator) = elements(k,i);
        %         triplet_list_cols(triplet_list_iterator) = elements(k,j);
        %         triplet_list_entries(triplet_list_iterator) = AK(i,j);
        %         triplet_list_iterator = triplet_list_iterator + 1;
        %     end
        % end
    end

    A = sparse(triplet_list_rows, triplet_list_cols, triplet_list_entries, num_nodes, num_nodes);
end
