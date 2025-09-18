% Solves simple elliptic problem in 1d using FEM and plots errors
clc;clear;close all;

% Imports
import common.*
import mesh.*
import fem1d.*

% Settings
c_handle_idx = 1;
u_exact_handle_idx = 3;
dof = 2;

% define function handles (real solution)   
% Cell array of 10 C^2 functions on [0,1]
syms x
cell_exact_fun = {
    x;                   % polynomial
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
    sin(10*x) + 2
};
c_handle = cell_c_fun{c_handle_idx};
u_exact_handle = cell_exact_fun{u_exact_handle_idx};
f_exact_handle = diff(c_handle*diff(-u_exact_handle, 1), 1)+ u_exact_handle;
du_exact_handle = diff(u_exact_handle, 1);
u_exact_handle = matlabFunction(u_exact_handle, 'vars', {x});
du_exact_handle = matlabFunction(du_exact_handle, 'vars', {x});
f_exact_handle = matlabFunction(f_exact_handle, 'vars', {x});
if isnumeric(c_handle)
    c_handle = @(x) c_handle + zeros(size(x));
else
    c_handle = matlabFunction(c_handle, 'vars', {x});
end

% initialize 
H_stepsizes = 2.^-(1:10);
errors = zeros(1, length(H_stepsizes));
for i = 1:length(H_stepsizes)
    % initialize mesh
    h = H_stepsizes(i);
    Mesh = mesh.MeshIntervalFEM1d([0,1], [h, h/100]);
    Mesh.dof = dof;
    Mesh.updatePet();
    [nodes, boundary_nodes_idx, elements] = Mesh.getPet();

    % set values from handles
    c_vals = c_handle(nodes(elements));
    f_values = f_exact_handle(nodes(elements));

    % assemble matrices
    num_nodes = length(nodes);
    A = fem1d.stiffnessMatrix1D(nodes, elements, c_vals);
    M = fem1d.massMatrix1D_v1(nodes, elements, c_handle);
    LHS = A + M;
    load_vec = fem1d.loadVector1D(nodes, elements, f_values);
    uh = zeros(num_nodes, 1);
    
    % boundary conditions
    uh(boundary_nodes_idx) = u_exact_handle(nodes(boundary_nodes_idx));
    
    % solve interior problem:
    interior_nodes_idx = 1:num_nodes;
    interior_nodes_idx(boundary_nodes_idx) = [];
    uh(interior_nodes_idx) = LHS(interior_nodes_idx, interior_nodes_idx)\(load_vec(interior_nodes_idx) - LHS(interior_nodes_idx, boundary_nodes_idx)*uh(boundary_nodes_idx));

    % calculate errors 
    errors(i) = fem1d.errors1D(elements, nodes, uh, du_exact_handle, u_exact_handle);
    disp("calculated h = "+ h + "  i = " + i + " cond(A) = " + condest(A(interior_nodes_idx, interior_nodes_idx)))
end

% plot solution
figure;
plot(nodes, uh, nodes, u_exact_handle(nodes))
xlabel("x")
ylabel("y")
legend("uh", "u\_exact")

% plot errors
figure;
loglog(H_stepsizes, H_stepsizes.^2, '--', H_stepsizes, H_stepsizes.^3, '--', H_stepsizes, errors);
xlabel('Step Size (H)');
ylabel('Error');
legend("h^2", "h^3", "L2")
title('Convergence of Errors');

