% Solves elliptic problem without exact solution and with non uniform mesh using SIPG 
% approximates errors using another numerical solution
clc;clear;close all;

% Imports
import mesh.*
import fem1d.*

% Settings
c_handle_idx = 1;
u_exact_handle_idx = 6;
sigma = 10;
dof = 2;
num_refinement_iterations = 10;

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
initial_meshsize = abs(boundary_nodes(1) - boundary_nodes(2))/3;
errors = zeros(1, num_refinement_iterations);
condition_B = zeros(1,num_refinement_iterations);

% set boundary conditions
boundary_cond = struct("values", [u_exact_handle(boundary_nodes(1)), du_exact_handle(boundary_nodes(2))], "lower_boundary_type", "dirichlet", "upper_boundary_type", "neumann");

% create initial mesh
H_meshsizes = zeros(1,num_refinement_iterations); 
H_meshsizes(1) = initial_meshsize;
Mesh = mesh.MeshIntervalDG1d(boundary_nodes, [initial_meshsize, initial_meshsize/100]);
Mesh.dof = dof;
Mesh.buildResonatorMesh([0.7,0.75; 0.3, 0.4], [0.2, 0.01]);
for i = 1:num_refinement_iterations

    % update mesh
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

    [uh, B] = dg1d.sip_1d_elliptic_solver(Mesh, boundary_cond, f_vals, c_vals, sigma);
    condition_B(i) = condest(B);

    % calculate errors 
    [errors(1,i),errors(2,i)] = fem1d.errors1DBetweenFemSol(nodes, elements, uh, u_exact_vals);
    disp("calculated uh for h = " + Mesh.h_max)

    % refine mesh
    H_meshsizes(i) = Mesh.h_max;
    Mesh.h_min = Mesh.h_min/refine_factor;
    refine_idx = true(size(elements, 1),1);
    Mesh.refineElementsByFact(refine_idx, refine_factor);
    Mesh.h_max = Mesh.h_max/refine_factor;
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

