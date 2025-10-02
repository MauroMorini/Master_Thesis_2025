function M = l2projectionMassMatrix1D(nodes_1,elements_1,nodes_2,elements_2)
        % V1 = span{phi_j^(1)}, V2 = span{phi_j^(2)} are the span of hat basis functions over nodes 1, 2
        % calculates mass matrix Mij = (phi_i^(2), phi_j^(1))_L2
        % note for projection here V2 is the projected space.
        % To project a solution u_2 in V_2 onto the space V_1 solve the system M*P(u_2) = u_2 
        % where P(u_2) is the projection of u_2 onto V_1 
        % only works currently if both spaces have same global dof
    arguments (Input)
        nodes_1           % (num_nodes1,1) node vector
        elements_1        % (num_el1, dof) connectivity matrix
        nodes_2           % (num_nodes2,1) node vector
        elements_2        % (num_el2, dof) connectivity matrix
    end

    arguments (Output)
        M               % (num_nodes, num_nodes) sparse mass matrix
    end

    % Initializations
    num_el_1 = size(elements_1, 1); 
    num_el_2 = size(elements_2, 1);
    num_nodes_1 = length(nodes_1);
    num_nodes_2 = length(nodes_2);
    dof = size(elements_1, 2);
    triplet_list_iterator = 1;
    c_vals_el = c_vals(elements);

    % preallocation
    M_max_entries = num_el_1*dof^2;
    triplet_list_rows = zeros(M_max_entries, 1);
    triplet_list_cols = zeros(M_max_entries, 1);
    triplet_list_entries = zeros(M_max_entries, 1);

    assert(size(elements_2, 2) == dof)

    % collect quadrature information
    [phi_val, ~, quad_weights] = common.getShapeFunctionValueMatrix(dof);
    
    % iterate over elements of V1
    for k = 1:size(elements_1, 1)
        K1 = nodes_1(elements_1(k,:));

        % h1 = abs(K1(1)-K1(end));
        % m1 = (K1(1)+K1(end))/2;

        % find elements in V2 which have a non empty intersection with
        % K
        p2El = [nodes_2(elements_2(:,1)), nodes_2(elements_2(:,end))];
        idxEmpty = p2El(:,2) <= K1(1) | K1(end) <= p2El(:,1);        % elements with empty intersection
        idxEl = find(~idxEmpty);
        for j = 1:length(idxEl)
            K2 = p2(elements_2(idxEl(j),:));

            % h2 = abs(K2(1)-K2(end));
            % m2 = (K2(1) + K2(end))/2;
            % 
            % % Intersection element
            % KInt = [max(K1(1),K2(1)), min(K1(2),K2(end))];
            % hInt = abs(KInt(1)-KInt(2));
            % mInt = (KInt(1) + KInt(2))/2;
            % f = @(x) [N0(FKInv(x,h2,m2))*N0(FKInv(x,h1,m1)), N0(FKInv(x,h2,m2))*N1(FKInv(x,h1,m1));
            %             N1(FKInv(x,h2,m2))*N0(FKInv(x,h1,m1)), N1(FKInv(x,h2,m2))*N1(FKInv(x,h1,m1))];
            % Mloc = hInt/6*(f(KInt(1)) + 4*f(mInt) + f(KInt(2)));

            Mloc = fem1d.projMassElementMatrix1D(K1, K2);
            M(t2(idxEl(j),:), t1(i,:)) = M(t2(idxEl(j),:), t1(i,:)) + Mloc;
        end          
    end
end