function Phi = evaluateLagrangeBarycentric(x_eval, weights_bary, basis_nodes)
    % evaluates the lagrange basis functions defined through the basis_nodes at the
    % evaluation nodes (nodes_eval) using barycentric interpolation
    arguments (Input)
        x_eval                  % (num_eval, 1) points to be evaluated
        weights_bary            % (num_nodes, 1) barycentric weights externally calculated 
        basis_nodes             % (num_nodes, 1) basis nodes through which lagrangian basis is defined
    end
    arguments (Output)
        Phi                     % (num_nodes, num_eval) [Phi]_{i,j} = L_i(x_eval_j) basis fun i at point j
    end
    
    D = x_eval.' - basis_nodes;

    % take note of points which coincide with basis nodes
    exact = (D==0);
    W = weights_bary ./ D;
    Phi = W ./ sum(W, 1);
    if any(exact(:))
        [exact_row_idx, exact_col_idx] = find(exact);
        for k = 1:length(exact_col_idx)
            Phi(:, exact_col_idx(k)) = 0;
            Phi(exact_row_idx(k),exact_col_idx(k)) = 1;
        end
    end
end