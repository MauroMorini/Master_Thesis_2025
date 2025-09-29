% Solves elliptic problem with given exact solution in 1d using SIPDG and plots errors
clc;clear;close all;

% Imports
import mesh.*
import fem1d.*

% Settings
c_handle_idx = 2;
u_exact_handle_idx = 9;
sigma = 10;
dof = 6;
overwrite_functions_bool = true;

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

% initialize parameters and preallocate
refine_factor = 2;
boundary_nodes = [0,1];
resonators_mat = [0.5, 0.6];
initial_meshsize = abs(boundary_nodes(1) - boundary_nodes(2))/20;

% set boundary conditions
boundary_cond = struct("values", [1, 10], "lower_boundary_type", "dirichlet", "upper_boundary_type", "neumann");

% overwrite function handles
if overwrite_functions_bool
    c_handle = @(x) 2*(x>=boundary_nodes(1) & x< resonators_mat(1)) + ...
                    0.5*(x<=boundary_nodes(2) & x>resonators_mat(2)) + ...
                    (11+10*sin(300*x)).*(x>=resonators_mat(1) & x<=resonators_mat(2));
    f_exact_handle = @(x) zeros(size(x));
    % f_exact_handle = @(x) 100*sin(3*x).*(x < 0.4 & x >= 0) + x.^2.*(x >= 0.4 & x <= 1);
end

% create initial mesh
Mesh = mesh.MeshIntervalDG1d(boundary_nodes, [initial_meshsize, initial_meshsize/100]);
Mesh.dof = dof;
Mesh.buildResonatorMesh(resonators_mat, [initial_meshsize/2, initial_meshsize/20]);

% update mesh
Mesh.updatePet();
[nodes, boundary_nodes_idx, elements] = Mesh.getPet();

% set values from handles
c_vals = c_handle(nodes);
f_vals = f_exact_handle(nodes);
g_vals = zeros(size(nodes));g_vals(1) = 1; g_vals(end) = 5;

[uh, B] = dg1d.sip_1d_elliptic_solver(Mesh, boundary_cond, f_vals, c_vals, sigma);

% refine mesh
Mesh.h_min = Mesh.h_min/refine_factor;
refine_idx = true(size(elements, 1),1);
Mesh.refineElementsByFact(refine_idx, refine_factor);
Mesh.h_max = Mesh.h_max/refine_factor;

% plot solution
figure;
plot(nodes, uh, nodes, c_vals)
xlabel("x")
ylabel("y")
legend("uh", "c")


