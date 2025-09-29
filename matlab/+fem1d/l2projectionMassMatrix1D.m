function M = l2projectionMassMatrix1D(nodes_1,elements_1,nodes_2,elements_2)
        % V1 = span{phi_j^(1)}, V2 = span{phi_j^(2)} are the span of hat basis functions over nodes 1, 2
        % calculates mass matrix Mij = (phi_i^(2), phi_j^(1))_L2
        % note for projection here V2 is the projected space.
        % To project a solution u_2 in V_2 onto the space V_1 solve the system M*P(u_2) = u_2 
        % where P(u_2) is the projection of u_2 onto V_1 
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
    N1 = length(p1);
    N2 = length(p2);
    M = sparse(N2, N1);

    % % functions
    % N0 = @(x) (1-x)/2;
    % N1 = @(x) (1+x)/2;
    % FKInv = @(x,h,m) 2*(x-m)/h;
    
    % iterate over elements of V1
    for i = 1:size(t1, 1)
        K1 = p1(t1(i,:));

        % h1 = abs(K1(1)-K1(end));
        % m1 = (K1(1)+K1(end))/2;

        % find elements in V2 which have a non empty intersection with
        % K
        p2El = [p2(t2(:,1)), p2(t2(:,end))];
        idxEmpty = p2El(:,2) <= K1(1) | K1(end) <= p2El(:,1);        % elements with empty intersection
        idxEl = find(~idxEmpty);
        for j = 1:length(idxEl)
            K2 = p2(t2(idxEl(j),:));

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