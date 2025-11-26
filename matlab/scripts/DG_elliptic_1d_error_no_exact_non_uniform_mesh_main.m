% Solves elliptic problem without exact solution and with non uniform mesh using SIPG 
% approximates errors using another numerical solution
clc;clear;close all;

% Imports
import mesh.*
import fem1d.*

% Settings
c_handle_idx = 2;
u_exact_handle_idx = 10;
dof = 2;
sigma = 10*dof^2;
num_refinement_iterations = 10;
overwrite_functions_bool = false;

% define function handles (real solution)   
% Cell array of 10 C^2 functions on [0,1]
syms x
cell_exact_fun = {
    x.^2;                   % polynomial
    x.^9;                   % polynomial
    sin(pi*x);              % sine
    cos(pi*x);              % cosine
    x.^2 .* (1-x);          % cubic-like
    exp(x);                 % exponential
    log(x+1);               % smooth on [0,1]
    sin(2*pi*x) + x;        % sine + linear
    x.^4 - x.^2;            % quartic
    exp(-x).*sin(5*x)       % damped oscillation
};
cell_c_fun = {
    1;
    (sin(10*x) + 2)*(x+1)
};
c_handle = cell_c_fun{c_handle_idx};
u_exact_handle = cell_exact_fun{u_exact_handle_idx};
f_exact_handle = diff(c_handle*diff(-u_exact_handle, 1), 1) + u_exact_handle;
du_exact_handle = diff(u_exact_handle, 1);
u_exact_handle = matlabFunction(u_exact_handle, 'vars', {x});
du_exact_handle = matlabFunction(du_exact_handle, 'vars', {x});
f_exact_handle = matlabFunction(f_exact_handle, 'vars', {x});
if isnumeric(c_handle)
    c_handle = @(x) c_handle + zeros(size(x));
else
    c_handle = matlabFunction(c_handle, 'vars', {x});
end

% initialize parameters and preallocate
refine_factor = 2;
boundary_nodes = [0,1];
initial_meshsize = abs(boundary_nodes(1) - boundary_nodes(2))/3;
errors = zeros(1, num_refinement_iterations);
condition_B = zeros(1,num_refinement_iterations);
numerical_solutions = cell(1, num_refinement_iterations);
exact_solution_struct = struct("u_handle", u_exact_handle, "du_handle", du_exact_handle, "type", "exact_solution");

% set boundary conditions
boundary_cond = struct("values", [u_exact_handle(boundary_nodes(1)), du_exact_handle(boundary_nodes(2))], "lower_boundary_type", "dirichlet", "upper_boundary_type", "neumann");

% create initial mesh
resonators_mat = [0.5, 0.6; 0.6, 0.9; 0.9, 1];
H_meshsizes = zeros(1,num_refinement_iterations); 
H_meshsizes(1) = initial_meshsize;
Mesh = mesh.MeshIntervalDG1d(boundary_nodes, [initial_meshsize, initial_meshsize/100]);
Mesh.dof = dof;
Mesh.buildResonatorMesh(resonators_mat, [initial_meshsize/2, initial_meshsize/20]);

% overwrite function handles
if overwrite_functions_bool
    c_handle = @(x) 2*(x>=boundary_nodes(1) & x< resonators_mat(1)) + ...
                    0.5*(x<=boundary_nodes(2) & x>resonators_mat(2)) + ...
                    (11+10*sin(300*x)).*(x>=resonators_mat(1) & x<=resonators_mat(2));
    % c_handle = @(x) ones(size(x));
    f_exact_handle = @(x) zeros(size(x));
    f_exact_handle = @(x) 100*sin(3*x).*(x < 0.4 & x >= 0) + x.^2.*(x >= 0.4 & x <= 1);
    boundary_cond = struct("values", [0, 1], "lower_boundary_type", "dirichlet", "upper_boundary_type", "neumann");
end

for i = 1:num_refinement_iterations

    % update mesh
    Mesh.updatePet();
    [nodes, ~, elements] = Mesh.getPet();

    % set values from handles
    c_vals = c_handle(nodes(elements));
    f_vals = f_exact_handle(nodes(elements));

    % solve system
    [uh, B] = dg1d.sip_1d_elliptic_solver(Mesh, boundary_cond, f_vals, c_vals, sigma);
    condition_B(i) = condest(B);
    disp("calculated uh for h = " + Mesh.h_max + "   ...   i = " + i)

    % collect solution
    numerical_solutions{i} = struct("mesh", copy(Mesh), "sol", uh, "type", "numerical_solution");

    % refine mesh
    H_meshsizes(i) = Mesh.h_max;
    if i < num_refinement_iterations
        Mesh.h_min = Mesh.h_min/refine_factor;
        refine_idx = true(size(elements, 1),1);
        Mesh.refineElementsByFact(refine_idx, refine_factor);
        Mesh.h_max = Mesh.h_max/refine_factor;
    end
end

% calculate high h-refined FEM sol
    Mesh.h_min = Mesh.h_min/8;
    refine_idx = true(size(elements, 1),1);
    Mesh.refineElementsByFact(refine_idx, 8);
    Mesh.h_max = Mesh.h_max/8;
    % update mesh
    Mesh.updatePet();
    [nodes, boundary_nodes_idx, elements] = Mesh.getPet();

    % % set values from handles
    % c_vals = c_handle(nodes);
    % f_vals = f_exact_handle(nodes);
    % [uh_ref, ~] = dg1d.sip_1d_elliptic_solver(Mesh, boundary_cond, f_vals, c_vals, sigma);
    % reference_sol_struct = struct("mesh", copy(Mesh), "sol", uh_ref, "type", "numerical_solution");
    % error_ref = errors_1d(reference_sol_struct, exact_solution_struct);
    % exact_solution_struct = reference_sol_struct;
    

% calculate errors
reference_struct = exact_solution_struct;
for i = 1:num_refinement_iterations-1
    if overwrite_functions_bool
        reference_struct = numerical_solutions{i+1};
    end
    [errors(1,i),errors(2,i)] = fem1d.errors_1d(numerical_solutions{i}, reference_struct);
    errors_1d_obj = fem1d.Errors1D(exact_solution_struct.u_handle, exact_solution_struct.du_handle, numerical_solutions{i}.sol, numerical_solutions{i}.mesh);
    errors_1d_obj.run();
    [errors(1,i),errors(2,i)] = errors_1d_obj.getErrors();
end

% plot solution
solution_idx = 1;
f = numerical_solutions{solution_idx}.mesh.plotDGsol(numerical_solutions{solution_idx}.sol);

% plot condition
figure;
loglog(H_meshsizes, H_meshsizes.^(-2), '--', H_meshsizes, condition_B);
xlabel('Step Size (H)');
ylabel('Condition of B');
legend("h^{-2}", "cond(B)")

% scale errors for better plot
errors_plot = errors;

% plot errors
figure;
loglog(H_meshsizes, H_meshsizes.^(dof-1), '--', H_meshsizes, H_meshsizes.^(dof),'--', H_meshsizes, errors_plot(1,:), H_meshsizes, errors_plot(2,:));
xlabel('Step Size (H)');
ylabel('Error');
legend("h^"+{dof-1}, "h^"+{dof}, "L2", "H1")
title("Convergence of Errors for P^"+(dof-1)+ "elements");