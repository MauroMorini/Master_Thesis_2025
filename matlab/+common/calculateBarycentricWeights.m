function barycentric_weights = calculateBarycentricWeights(nodes, elements)
    % calculates barycentric_weights for each node for a fem (node, element) tuple
    arguments (Input)
        nodes                   % (num_nodes, 1) node matrix
        elements                % (num_el, dof) connectivity matrix
    end
    arguments (Output)
        barycentric_weights     % (num_nodes, 1)
    end
    dof = size(elements,2);
    nodes_loc = nodes(elements);

    % in case we have only 1 element nodes_loc should be a row vector 
    if size(elements, 1) == 1
        nodes_loc = nodes_loc.';
    end
    barycentric_weights = ones(size(nodes_loc));
    for i = 1:dof
        for j = 1:dof
            if i == j
                continue
            end
            barycentric_weights(:,i) = barycentric_weights(:, i) * 1./(nodes_loc(:,i) - nodes_loc(:,j));
        end
    end
    barycentric_weights = barycentric_weights';
    barycentric_weights = barycentric_weights(:);
end
