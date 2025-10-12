% Solves simple elliptic problem in 1d using FEM and plots errors
clc;clear;close all;

% Imports
import common.*
import mesh.*
import fem1d.*

% Settings
c_handle_idx = 1;
u_exact_handle_idx = 10;
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
H_meshsizes = 2.^-(1:10);
errors = zeros(2, length(H_meshsizes));
for i = 1:length(H_meshsizes)
    % initialize mesh
    h = H_meshsizes(i);
    Mesh = mesh.MeshIntervalFEM1d([0,1], [h, h/100]);
    Mesh.dof = dof;
    Mesh.updatePet();
    [nodes, boundary_nodes_idx, elements] = Mesh.getPet();

    % set values from handles
    c_vals = c_handle(nodes);
    f_values = f_exact_handle(nodes);
    u_exact_vals = u_exact_handle(nodes);
    du_exact_vals = du_exact_handle(nodes);

    % assemble matrices
    num_nodes = length(nodes);
    A = fem1d.stiffnessMatrix1D(nodes, elements, c_vals);
    M = fem1d.massMatrix1D(nodes, elements, c_vals);
    LHS = A + M;
    load_vec = fem1d.loadVector1D(nodes, elements, f_values);
    uh = zeros(num_nodes, 1);
    
    % boundary conditions
    uh(boundary_nodes_idx) = u_exact_handle(nodes(boundary_nodes_idx));
    
    % solve interior problem:
    interior_nodes_idx = 1:num_nodes;
    interior_nodes_idx(boundary_nodes_idx) = [];
    uh_loc = LHS(interior_nodes_idx, interior_nodes_idx)\(load_vec(interior_nodes_idx) - LHS(interior_nodes_idx, boundary_nodes_idx)*uh(boundary_nodes_idx));
    system_matrix = LHS(interior_nodes_idx, interior_nodes_idx);
    system_vector = (load_vec(interior_nodes_idx) - LHS(interior_nodes_idx, boundary_nodes_idx)*uh(boundary_nodes_idx));

        % rescale system 
    scaling_factor = max(diag(system_matrix));
    system_matrix = system_matrix/scaling_factor;
    system_vector = system_vector/scaling_factor;

     % solve system using pcg
    tol = Mesh.h_max^Mesh.dof;
    maxit = numel(system_matrix);
    try
        L = ichol(system_matrix, struct('michol', 'on'));
    catch exception
        warning("ichol has failed shifted ichol is used")
        alpha = max(sum(abs(system_matrix),2)./diag(system_matrix));
        L = ichol(system_matrix, struct('type','ict','droptol',1e-3,'diagcomp',alpha));
    end
    [uh_loc, failed_to_converge_flag] = pcg(system_matrix,system_vector,tol,maxit,L,L');
    if failed_to_converge_flag
        warning("pcg has failed to converge normal mldivide is used")
        uh_loc = system_matrix\system_vector;
    end
    % uh_loc = LHS(interior_nodes_idx, interior_nodes_idx)\(load_vec(interior_nodes_idx) - LHS(interior_nodes_idx, boundary_nodes_idx)*uh(boundary_nodes_idx));
    uh(interior_nodes_idx) = uh_loc;

    % calculate errors 
    [errors(1,i), errors(2,i)] = fem1d.errors1DWithExactSol(nodes, elements, uh, u_exact_vals, du_exact_vals);
    disp("calculated h = "+ h + "  i = " + i + " cond(A) = " + condest(A(interior_nodes_idx, interior_nodes_idx)))
end

% plot solution
figure;
plot(nodes, uh, nodes, u_exact_vals)
xlabel("x")
ylabel("y")
legend("uh", "u\_exact")

% scale errors for better plot
errors_plot = errors;
if errors(1,1) > 1e-10
    errors_plot(1,:) = errors(1,:)/(errors(1,1))*H_meshsizes(1)^(dof);
end
if errors(2,1) > 1e-10
    errors_plot(2,:) = errors(2,:)/(errors(2,1))*H_meshsizes(1)^(dof-1);
end

% plot errors
figure;
loglog(H_meshsizes, H_meshsizes.^(dof-1), '--', H_meshsizes, H_meshsizes.^(dof),'--', H_meshsizes, errors_plot(1,:), H_meshsizes, errors_plot(2,:));
xlabel('Step Size (H)');
ylabel('Error');
legend("h^"+{dof-1}, "h^"+{dof}, "L2", "H1")
title("Convergence of Errors for P^"+(dof-1)+ "elements");

