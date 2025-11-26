function [l2_error, h1_error] = errors_1d(numerical_sol_struct, exact_sol_struct)
    % calculates the error in L^2 and H^1 norm for FEM. Can deal with FEM solution of high h-refinement 
    % as exact solution or just taking the exact solution and it's derivative as function handles 
    % note that for FEM solution as reference the H^1 norm is too good!
    % and FEM solution is L2 projected, so the FEM-space should be a subset of the numerical solution space
    arguments (Input)
        numerical_sol_struct    % struct: {"mesh": MeshIntervalDG1d, "sol": (num_nodes, 1) sol vector, "type": string "numerical_solution"}
        exact_sol_struct        % struct either like numerical sol (for FEM exact sol), or struct: {"u_handle": @(x) u, "du_handle": @(x) du, "type": "exact_solution"}
    end
    arguments (Output)
        l2_error                % scalar value of L^2 error over domain
        h1_error                % scalar value of H^1 error over domain     
    end
    
    % extract numerical solution
    [nodes_num, ~, elements_num] = numerical_sol_struct.mesh.getPet();
    uh = numerical_sol_struct.sol;

    if exact_sol_struct.type == "numerical_solution"
        [nodes_exact, ~, elements_exact] = exact_sol_struct.mesh.getPet();
        u_exact_vals = exact_sol_struct.sol;
        u_exact_vals_projected = fem1d.l2projection1D(nodes_exact, elements_exact, u_exact_vals, nodes_num, elements_num, true);
        [l2_error, h1_error] = fem1d.errors1DBetweenFemSol(nodes_num, elements_num, uh, u_exact_vals_projected);
    elseif exact_sol_struct.type == "exact_solution"
        u_exact_vals = exact_sol_struct.u_handle(nodes_num);
        du_exact_vals = exact_sol_struct.du_handle(nodes_num);
        [l2_error, h1_error] = fem1d.errors1DWithExactSol(nodes_num, elements_num, uh, u_exact_vals, du_exact_vals);
    else
        error("wrong type of exact solution:  " + exact_sol_struct.type)
    end

end