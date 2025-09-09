% Solves simple elliptic problem in 1d using DG and plots errors
clc;clear;close all;

% Imports
import mesh.*
import fem1d.*

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
u_exact_handle = cell_exact_fun{6};
f_exact_handle = diff(-u_exact_handle, 2);
du_exact_handle = diff(u_exact_handle, 1);
u_exact_handle = matlabFunction(u_exact_handle, 'vars', {x});
du_exact_handle = matlabFunction(du_exact_handle, 'vars', {x});
f_exact_handle = matlabFunction(f_exact_handle, 'vars', {x});
c_handle = @(x) ones(size(x));

% initialize 
H_stepsizes = 2.^-(1:10);
errors = zeros(1, length(H_stepsizes));
for i = 1:length(H_stepsizes)
    % initialize mesh
    h = H_stepsizes(i);
    Mesh = Mesh1d([0,1], [h, h/100]);
    [nodes, boundary_nodes_idx, elements] = Mesh.getPet();
    
    % assemble matrices
    num_nodes = length(nodes);
    A = fem1d.stiffnessMatrix1D(nodes', elements, c_handle);
    load_vec = fem1d.loadVectorLinear1D(nodes', elements, f_exact_handle);
    uh = zeros(num_nodes, 1);
    
    % boundary conditions
    uh(boundary_nodes_idx) = u_exact_handle(nodes(boundary_nodes_idx));
    
    % solve interior problem:
    interior_nodes_idx = 1:num_nodes;
    interior_nodes_idx(boundary_nodes_idx) = [];
    uh(interior_nodes_idx) = A(interior_nodes_idx, interior_nodes_idx)\(load_vec(interior_nodes_idx) - A(interior_nodes_idx, boundary_nodes_idx)*uh(boundary_nodes_idx));

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
loglog(H_stepsizes, H_stepsizes.^2, '--', H_stepsizes, errors);
xlabel('Step Size (H)');
ylabel('Error');
legend("h²", "L2")
title('Convergence of Errors');

% interior nodes
num_nodes = length(nodes);
interior_nodes_idx = 1:num_nodes;
interior_nodes_idx(boundary_nodes_idx) = [];

u_h_zero_bc = A(interior_nodes_idx, interior_nodes_idx)\load_v(interior_nodes_idx);

