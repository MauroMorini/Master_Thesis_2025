% Solves elliptic problem without exact solution and with non uniform mesh using SIPG 
% approximates errors using another numerical solution
clc;clear;close all;

% Imports
import mesh.*
import fem1d.*

% Settings
c_handle_idx = 1;
u_exact_handle_idx = 1;
dof = 3;
num_refinement_iterations = 9;

% define function handles (real solution)   
% Cell array of 10 C^2 functions on [0,1]
syms x
cell_exact_fun = {
    x;                   % polynomial
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
    sin(10*x) + 2
};
c_handle = cell_c_fun{c_handle_idx};
u_exact_handle = cell_exact_fun{u_exact_handle_idx};
f_exact_handle = diff(c_handle*diff(-u_exact_handle, 1), 1);
du_exact_handle = diff(u_exact_handle, 1);
u_exact_handle = matlabFunction(u_exact_handle, 'vars', {x});
du_exact_handle = matlabFunction(du_exact_handle, 'vars', {x});
f_exact_handle = matlabFunction(f_exact_handle, 'vars', {x});
if isnumeric(c_handle)
    c_handle = @(x) c_handle + zeros(size(x));
else
    c_handle = matlabFunction(c_handle, 'vars', {x});
end

u_exact_handle = @(x) ones(size(x));
du_exact_handle = @(x) zeros(size(x));

% initialize parameters and preallocate
refine_factor = 2;
boundary_nodes = [0,1];
initial_meshsize = abs(boundary_nodes(1) - boundary_nodes(2))/2;
errors = zeros(1, num_refinement_iterations);
exact_solution_struct = struct("u_handle", u_exact_handle, "du_handle", du_exact_handle, "type", "exact_solution");

% create initial mesh
H_meshsizes = zeros(1,num_refinement_iterations); 
H_meshsizes(1) = initial_meshsize;
Mesh = mesh.MeshIntervalDG1d(boundary_nodes, [initial_meshsize, initial_meshsize/100]);
Mesh.dof = dof;
% Mesh.buildResonatorMesh([0.7,0.75; 0.3, 0.4], [0.2, 0.01]);

for i = 1:num_refinement_iterations

    % update mesh
    Mesh.updatePet();
    [nodes, ~, elements] = Mesh.getPet();

    u_exact_vals = u_exact_handle(nodes);
    
    % collect solution
    numerical_solutions{i} = struct("mesh", copy(Mesh), "sol", u_exact_vals, "type", "numerical_solution");

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

    uh_ref = u_exact_handle(nodes);
    reference_sol_struct = struct("mesh", copy(Mesh), "sol", uh_ref, "type", "numerical_solution");

% project and calculate errors
for i = 1:num_refinement_iterations-1
    [nodes_origin,~,elements_origin] = numerical_solutions{i}.mesh.getPet();
    [nodes_target,~,elements_target] = numerical_solutions{i+1}.mesh.getPet();
    values_origin = numerical_solutions{i}.sol;
    u_projected = l2projection1D(nodes_origin, elements_origin, values_origin, nodes_target, elements_target, true);
    projection_solution = struct("mesh", numerical_solutions{i+1}.mesh, "sol", u_projected, "type", "numerical_solution");
    [errors(1,i),~] = fem1d.errors1D(projection_solution, exact_solution_struct);
end

% plot errors
figure;
loglog(H_meshsizes, H_meshsizes.^(dof-1), '--', H_meshsizes, H_meshsizes.^(dof),'--', H_meshsizes, errors(1,:));
xlabel('Step Size (H)');
ylabel('Error');
legend("h^"+{dof-1}, "h^"+{dof}, "L2")
title("Convergence of Errors for P^"+(dof-1)+ "elements");