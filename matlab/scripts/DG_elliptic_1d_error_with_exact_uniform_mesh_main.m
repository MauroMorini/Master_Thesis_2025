% Solves elliptic problem with given exact solution in 1d using SIPDG and plots errors
clc;clear;close all;

% Imports
import mesh.*
import fem1d.*

% Settings
c_handle_idx = 3;
u_exact_handle_idx = 10;
dof = 3;
sigma = 10*dof^2;

% define function handles (real solution)   
% Cell array of 10 C^2 functions on [0,1]
syms x
cell_exact_fun = {
    x.^2;                   % polynomial
    x.^3;                   % polynomial
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
    sin(10*x) + 2;
    x.^(dof-1) + 1
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
H_meshsizes = 2.^-(2:10);
errors = zeros(3, length(H_meshsizes));
condition_B = zeros(1,length(H_meshsizes));
exact_solution_struct = struct("u_handle", u_exact_handle, "du_handle", du_exact_handle, "type", "exact_solution");
boundary_nodes = [0,1];
numerical_solution = cell(1, length(H_meshsizes));

% set boundary conditions
boundary_cond = struct("values", [u_exact_handle(boundary_nodes(1)), u_exact_handle(boundary_nodes(2))], "lower_boundary_type", "dirichlet", "upper_boundary_type", "dirichlet");


for i = 1:length(H_meshsizes)
    % initialize mesh
    h = H_meshsizes(i);
    Mesh = mesh.MeshIntervalDG1d([0,1], [2*h, 2*h/100]);
    Mesh.createUniformMesh(h);
    Mesh.dof = dof;
    Mesh.updatePet();
    [nodes, boundary_nodes_idx, elements] = Mesh.getPet();

    % set values from handles
    c_vals = c_handle(nodes(elements));
    f_vals = f_exact_handle(nodes(elements));
    u_exact_vals = u_exact_handle(nodes);
    du_exact_vals = du_exact_handle(nodes);
    
    % solve system
    [uh, B] = dg1d.sip_1d_elliptic_solver(Mesh, boundary_cond, f_vals, c_vals, sigma);
    condition_B(i) = condest(B);
    numerical_solution{i} = struct("mesh", copy(Mesh), "sol", uh, "type", "numerical_solution");
    disp("calculated uh for h = " + Mesh.h_max + "   ...   i = " + i)

    % calculate errors 
    [errors(1,i),errors(2,i)] = fem1d.errors1D(numerical_solution{i}, exact_solution_struct);
    % errors(3,i) = dg1d.energyNormError1D(nodes, elements, uh, c_vals, sigma, u_exact_vals, du_exact_vals);
end

% plot solution
solution_idx = 2;
f = numerical_solution{solution_idx}.mesh.plotDGsol(numerical_solution{solution_idx}.sol);

% plot condition
figure;
loglog(H_meshsizes, H_meshsizes.^(-2), '--', H_meshsizes, condition_B);
xlabel('Step Size (H)');
ylabel('Condition of B');
legend("h^{-2}", "cond(B)")

% plot errors
line_width = 2;
figure;
loglog(H_meshsizes, H_meshsizes.^(dof-1), '--','LineWidth', line_width);
hold on
loglog(H_meshsizes, H_meshsizes.^(dof), '--', 'LineWidth', line_width);
loglog(H_meshsizes, errors(1,:), 'LineWidth', line_width);
loglog(H_meshsizes, errors(2,:), 'LineWidth', line_width);
loglog(H_meshsizes, errors(3,:), 'LineWidth', line_width);
hold off
xlabel('Step Size (H)');
ylabel('Error');
legend("h^"+(dof-1), "h^"+(dof), "L2", "H1", "energy")
title("Convergence of Errors for P^"+(dof-1)+ "elements");

