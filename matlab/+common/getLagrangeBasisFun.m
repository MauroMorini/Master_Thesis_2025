function [phi_cell, dphi_cell] = getLagrangeBasisFun(nodes)
    % depending on degrees of freedom provided collects a cell array with function handles 
    % of basis Lagrange functions defined on the reference element K = (-1,1)
    arguments (Input)
        nodes             % (1,dof) node vector 
    end
    arguments (Output)
        phi_cell        % (1,dof) cell array containing function handles of basis functions 
        dphi_cell       % (1,dof) cell array containing function handles of derivatives of basis functions
    end

    % initialization
    dof = length(nodes);
    syms x
    phi_cell = cell(1, dof);
    dphi_cell = cell(1, dof);

    for function_node_idx = 1:dof
        phi = 1;
        for inner_node_idx = 1:dof
            if function_node_idx == inner_node_idx
                continue
            end
            phi = phi*(x - nodes(inner_node_idx))/(nodes(function_node_idx) - nodes(inner_node_idx));
        end
        dphi = diff(phi, 1);
        phi = matlabFunction(phi, 'vars',{x});
        dphi = matlabFunction(dphi, 'vars', {x});
        phi_cell{function_node_idx} = phi;
        dphi_cell{function_node_idx} = dphi;   
    end
end