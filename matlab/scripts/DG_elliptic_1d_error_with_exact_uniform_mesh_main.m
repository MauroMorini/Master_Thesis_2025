% Solves elliptic problem with given exact solution in 1d using SIPDG and plots errors
clc;clear;close all;

% Imports
import mesh.*
import fem1d.*

% Settings
c_handle_idx = 1;
u_exact_handle_idx = 10;
dof = 3;
sigma = 100*dof^2;

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
H_meshsizes = 2.^-(2:11);
errors = zeros(1, length(H_meshsizes));
condition_B = zeros(1,length(H_meshsizes));

for i = 1:length(H_meshsizes)
    % initialize mesh
    h = H_meshsizes(i);
    Mesh = mesh.MeshIntervalDG1d([0,1], [h, h/100]);
    Mesh.dof = dof;
    Mesh.updatePet();
    [nodes, boundary_nodes_idx, elements] = Mesh.getPet();

    % set values from handles
    c_vals = c_handle(nodes);
    f_vals = f_exact_handle(nodes);
    u_exact_vals = u_exact_handle(nodes);
    du_exact_vals = du_exact_handle(nodes);
    
    % assemble matrices
    num_nodes = length(nodes);
    B = dg1d.sipdgMatrix1D(nodes, elements, c_vals, sigma);
    rhs_vector = dg1d.sipdgDirichletLoadVector1D(nodes, elements, f_vals, c_vals, u_exact_vals, sigma);
    
    % solve system
    uh = B\rhs_vector;
    condition_B(i) = condest(B);

    % calculate errors 
    [errors(1,i),errors(2,i)] = fem1d.errors1DWithExactSol(nodes, elements, uh, u_exact_vals, du_exact_vals);
    disp("calculated uh for h = " + h)
end

% plot solution
figure;
plot(nodes, uh, nodes, u_exact_handle(nodes))
xlabel("x")
ylabel("y")
legend("uh", "u\_exact")

% plot condition
figure;
loglog(H_meshsizes, H_meshsizes.^(-2), '--', H_meshsizes, condition_B);
xlabel('Step Size (H)');
ylabel('Condition of B');
legend("h^{-2}", "cond(B)")

% plot errors
figure;
loglog(H_meshsizes, H_meshsizes.^(dof-1), '--', H_meshsizes, H_meshsizes.^(dof),'--', H_meshsizes, errors(1,:),H_meshsizes, errors(2,:));
xlabel('Step Size (H)');
ylabel('Error');
legend("h^"+(dof-1), "h^"+(dof), "L2", "H1")
title("Convergence of Errors for P^"+(dof-1)+ "elements");

