% This is the main script of chapter 1 numerical results. It is subdivided
% into sections just like the Numerical Results section in chapter 1.
% The sections are self-contained and can be ran individually, therefore
% they also share a lot of repeating code.
clc; clear; close all;

%% Rate of Convergence
% This section replicates the desired rates of convergences play with the
% settings to get different variations of the experiment

% Settings
c_index = 2;            % index of the coefficient c choose int 1-3
u_exact_handle_idx = 2; % index of the exact sol u choose int 1-10
is_resonator = false;   % choose if initial mesh is uniform or not
dof = 2;                % dof = r+1 the polynomial degree (for P1 elements do dof = 2)
sigma = 10*dof^2;       % penalization parameter sigma 

% additional settings 
num_ref = 11;
initial_meshsize = 0.25;
refine_factor = 2;
boundary_nodes = [0, 1];

% functions
syms x
cell_exact_fun = {
    x.^(dof-1);
    exp(-x).*sin(5*x);
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
    x*0 + 1;
    sin(x) + 2;
    x.^(dof-1) + 1
};
c_handle = cell_c_fun{c_index};
u_exact_handle = cell_exact_fun{u_exact_handle_idx};
f_exact_handle = -diff(c_handle*diff(u_exact_handle, 1), 1);
du_exact_handle = diff(u_exact_handle, 1);
u_exact_handle = matlabFunction(u_exact_handle, 'vars', {x});
du_exact_handle = matlabFunction(du_exact_handle, 'vars', {x});
f_exact_handle = matlabFunction(f_exact_handle, 'vars', {x});
c_handle = matlabFunction(c_handle, 'vars', {x});

c_handle = @(x) c_handle(x) + zeros(size(x));
u_exact_handle = @(x) u_exact_handle(x) + zeros(size(x)); 
du_exact_handle = @(x) du_exact_handle(x) + zeros(size(x)); 
f_exact_handle = @(x) f_exact_handle(x) + zeros(size(x));

boundary_cond = struct("values", [u_exact_handle(boundary_nodes(1)), u_exact_handle(boundary_nodes(2))], "lower_boundary_type", "dirichlet", "upper_boundary_type", "dirichlet");

% initialization
errors = zeros(2, num_ref);
meshsizes = zeros(1, num_ref);
waveguide = mesh.MeshIntervalDG1d(boundary_nodes, [2*initial_meshsize, initial_meshsize/50]);
waveguide.createUniformMesh(initial_meshsize);
if is_resonator
    waveguide.buildResonatorMesh([0.3, 0.7], [initial_meshsize, initial_meshsize/5]);
end
waveguide.dof = dof;
waveguide.updatePet();

for i = 1:num_ref
    if i > 1
        waveguide.refineAll(refine_factor);
        waveguide.updatePet();
    end

    % set values from handles
    quad_mesh = copy(waveguide);
    quad_mesh.dof = quad_mesh.dof + 1;
    quad_mesh.updatePet();
    [quad_nodes, ~, quad_elements] = quad_mesh.getPet();
    c_vals = c_handle(quad_nodes(quad_elements));
    f_vals = f_exact_handle(quad_nodes(quad_elements));

    % solve system
    [uh, B] = dg1d.sip_1d_elliptic_solver(waveguide, boundary_cond, f_vals, c_vals, sigma);
    condition_B(i) = condest(B);

    % calculate errors
    error_obj = fem1d.Errors1D(u_exact_handle, du_exact_handle, uh, waveguide);
    error_obj.initialize_dg_settings(c_handle, sigma);
    error_obj.run();

    [errors(1,i), errors(2,i), errors(3,i)] = error_obj.getErrors();
    meshsizes(i) = waveguide.h_max;
end

% plot errors
dof_plot = dof;
line_width = 2;
figure;
loglog(meshsizes, meshsizes.^(dof_plot-1), '--','LineWidth', line_width);
hold on
loglog(meshsizes, meshsizes.^(dof_plot), '--', 'LineWidth', line_width);
loglog(meshsizes, errors(1,:),'^-' ,'LineWidth', line_width);
loglog(meshsizes, errors(2,:),'o-' ,'LineWidth', line_width);
loglog(meshsizes, errors(3,:), 'x-','LineWidth', line_width);
hold off
xlabel('meshsize h');
ylabel('error');
legend("h^"+(dof_plot-1), "h^"+(dof_plot), "L2", "H1", "energy")
title("Convergence of Errors for P^"+(dof-1)+ "elements");
grid on;
cr = fem1d.Errors1D.interpolate_convergence_rates(meshsizes, errors);
fprintf('Interpolated Convergence Rates ------------- \nl2: %g,  h1: %g,  energy: %g \n\n', cr(2,1),cr(2,2),cr(2,3));

%% Visualization of the SIPG Solution
% Here we visualize the numerical solution after a certain set of
% refinements

% Settings
c_index = 2;            % index of the coefficient c choose int 1-3
u_exact_handle_idx = 2; % index of the exact sol u choose int 1-10
is_resonator = false;   % choose if initial mesh is uniform or not
dof = 3;                % dof = r+1 the polynomial degree (for P1 elements do dof = 2)
sigma = 1.1;       % penalization parameter sigma 
num_ref = 1;            % number of refinement cycles before plotting

% additional settings 
initial_meshsize = 0.25;
refine_factor = 2;
boundary_nodes = [0, 1];

% functions
syms x
cell_exact_fun = {
    x.^(dof-1);
    exp(-x).*sin(5*x);
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
    x*0 + 1;
    sin(x) + 2;
    x.^(dof-1) + 1
};
c_handle = cell_c_fun{c_index};
u_exact_handle = cell_exact_fun{u_exact_handle_idx};
f_exact_handle = -diff(c_handle*diff(u_exact_handle, 1), 1);
du_exact_handle = diff(u_exact_handle, 1);
u_exact_handle = matlabFunction(u_exact_handle, 'vars', {x});
du_exact_handle = matlabFunction(du_exact_handle, 'vars', {x});
f_exact_handle = matlabFunction(f_exact_handle, 'vars', {x});
c_handle = matlabFunction(c_handle, 'vars', {x});

c_handle = @(x) c_handle(x) + zeros(size(x));
u_exact_handle = @(x) u_exact_handle(x) + zeros(size(x)); 
du_exact_handle = @(x) du_exact_handle(x) + zeros(size(x)); 
f_exact_handle = @(x) f_exact_handle(x) + zeros(size(x));

boundary_cond = struct("values", [u_exact_handle(boundary_nodes(1)), u_exact_handle(boundary_nodes(2))], "lower_boundary_type", "dirichlet", "upper_boundary_type", "dirichlet");

% initialization
waveguide = mesh.MeshIntervalDG1d(boundary_nodes, [2*initial_meshsize, initial_meshsize/50]);
waveguide.createUniformMesh(initial_meshsize);
if is_resonator
    waveguide.buildResonatorMesh([0.3, 0.7], [initial_meshsize, initial_meshsize/5]);
end
waveguide.dof = dof;
waveguide.updatePet();

% refine mesh
for i = 1:num_ref
    waveguide.refineAll(2);
    waveguide.updatePet();
end

% compute numerical solution
[nodes, ~, elements] = waveguide.getPet();

% set values from handles
quad_mesh = copy(waveguide);
quad_mesh.dof = quad_mesh.dof + 1;
quad_mesh.updatePet();
[quad_nodes, ~, quad_elements] = quad_mesh.getPet();
c_vals = c_handle(quad_nodes(quad_elements));
f_vals = f_exact_handle(quad_nodes(quad_elements));

% solve system
[uh, B] = dg1d.sip_1d_elliptic_solver(waveguide, boundary_cond, f_vals, c_vals, sigma);
condition_B(i) = condest(B);

f = figure;
waveguide.plotDGsol(uh, f);
ylim([-0.5, 1])


%% Influence of the Quadrature Rule on the Convergence Rate
% here we compute again two sets of error rates, one with a 
% higher order quadrature rule and one where the quadrature nodes coincide 
% with the basis nodes

% Settings
c_index = 4;            % index of the coefficient c choose int 1-3
u_exact_handle_idx = 2; % index of the exact sol u choose int 1-10
is_resonator = false;   % choose if initial mesh is uniform or not
dof = 2;                % dof = r+1 the polynomial degree (for P1 elements do dof = 2)
sigma = 10*dof^2;       % penalization parameter sigma 

% additional settings 
num_ref = 11;
initial_meshsize = 0.25;
refine_factor = 2;
boundary_nodes = [0, 1];

% functions
syms x
cell_exact_fun = {
    x.^(dof-1);
    exp(-x).*sin(5*x);
    sin(pi*x);              % sine
    cos(pi*x);              % cosine
    x.^2 .* (1-x);          % cubic-like
    exp(x);                 % exponential
    log(x+1);               % smooth on [0,1]
    sin(2*pi*x) + x;        % sine + linear
    x.^4 - x.^2;            % quartic
    exp(-x).*sin(50*x)       % damped oscillation
};
cell_c_fun = {
    x*0 + 1;
    sin(x) + 2;
    x.^(dof-1) + 1;
    sin(20*x) + 2
};
c_handle = cell_c_fun{c_index};
u_exact_handle = cell_exact_fun{u_exact_handle_idx};
f_exact_handle = -diff(c_handle*diff(u_exact_handle, 1), 1) + u_exact_handle;
du_exact_handle = diff(u_exact_handle, 1);
u_exact_handle = matlabFunction(u_exact_handle, 'vars', {x});
du_exact_handle = matlabFunction(du_exact_handle, 'vars', {x});
f_exact_handle = matlabFunction(f_exact_handle, 'vars', {x});
c_handle = matlabFunction(c_handle, 'vars', {x});

c_handle = @(x) c_handle(x) + zeros(size(x));
u_exact_handle = @(x) u_exact_handle(x) + zeros(size(x)); 
du_exact_handle = @(x) du_exact_handle(x) + zeros(size(x)); 
f_exact_handle = @(x) f_exact_handle(x) + zeros(size(x));

boundary_cond = struct("values", [u_exact_handle(boundary_nodes(1)), u_exact_handle(boundary_nodes(2))], "lower_boundary_type", "dirichlet", "upper_boundary_type", "dirichlet");

% initialization
errors_mass_lump  = zeros(2, num_ref);
errors_exact = zeros(2, num_ref);
meshsizes = zeros(1, num_ref);
waveguide = mesh.MeshIntervalDG1d(boundary_nodes, [2*initial_meshsize, initial_meshsize/50]);
waveguide.createUniformMesh(initial_meshsize);
if is_resonator
    waveguide.buildResonatorMesh([0.3, 0.7], [initial_meshsize, initial_meshsize/5]);
end
waveguide.dof = dof;
waveguide.updatePet();

for i = 1:num_ref
    if i > 1
        waveguide.refineAll(2);
        waveguide.updatePet();
    end

    [nodes, ~, elements] = waveguide.getPet();

    % set values from handles
    quad_mesh = copy(waveguide);
    quad_mesh.dof = quad_mesh.dof + 1;
    quad_mesh.updatePet();
    [quad_nodes, ~, quad_elements] = quad_mesh.getPet();
    c_vals_quad = c_handle(quad_nodes(quad_elements));
    f_vals_quad = f_exact_handle(quad_nodes(quad_elements));
    c_vals = c_handle(nodes(elements));
    f_vals = f_exact_handle(nodes(elements));

    % solve system
    [uh_mass_lump, ~] = dg1d.sip_1d_elliptic_solver(waveguide, boundary_cond, f_vals, c_vals, sigma, true);
    [uh_exact, ~] = dg1d.sip_1d_elliptic_solver(waveguide, boundary_cond, f_vals_quad, c_vals_quad, sigma, true);
    
    % calculate errors
    error_obj_mass_lump = fem1d.Errors1D(u_exact_handle, du_exact_handle, uh_mass_lump , waveguide);
    error_obj_mass_lump.initialize_dg_settings(c_handle, sigma);
    error_obj_mass_lump.additional_quadrature_dof = 3;
    error_obj_mass_lump.run();

    error_obj_exact = fem1d.Errors1D(u_exact_handle, du_exact_handle, uh_exact , waveguide);
    error_obj_exact.initialize_dg_settings(c_handle, sigma);
    error_obj_exact.additional_quadrature_dof = 3;
    error_obj_exact.run();

    [errors_mass_lump(1,i), errors_mass_lump(2,i), errors_mass_lump(3,i)] = error_obj_mass_lump.getErrors();
    [errors_exact(1,i), errors_exact(2,i), errors_exact(3,i)] = error_obj_exact.getErrors();
    meshsizes(i) = waveguide.h_max;
end

% calculate relative distance
errors_difference = (errors_mass_lump - errors_exact);

% plot errors
dof_plot = dof;
line_width = 2;
figure;
loglog(meshsizes, errors_mass_lump(1,:),'^-' ,'LineWidth', line_width);
hold on
loglog(meshsizes, errors_exact(1,:),'^-' ,'LineWidth', line_width);
loglog(meshsizes, errors_mass_lump(2,:),'o-' ,'LineWidth', line_width);
loglog(meshsizes, errors_exact(2,:),'o-' ,'LineWidth', line_width);
loglog(meshsizes, errors_mass_lump(3,:), 'x-','LineWidth', line_width);
loglog(meshsizes, errors_exact(3,:), 'x-','LineWidth', line_width);
hold off
xlabel('meshsize h');
ylabel('error');
legend("L2 mass lump", "L2 exact","H1 mass lump", "H1 exact","energy mass lump", "energy exact")
title("Error rates with higher and lower quadrature rule for P^"+(dof-1)+ "elements");
grid on

figure;
loglog(meshsizes, errors_difference(1,:),'^-' ,'LineWidth', line_width)
hold on
loglog(meshsizes, errors_difference(2,:),'o-' ,'LineWidth',line_width)
loglog(meshsizes, errors_difference(3,:), 'x-','LineWidth', line_width)
hold off
xlabel('meshsize h');
ylabel('error');
legend("L2","H1","energy")
title("Comparison of error rates, lower quad error - higher quad error for P^"+(dof-1)+ "elements")
grid on

%% Modeling an Inhomogeneous Membrane
% here we now model a problem only knowing the external forces and not the exact solution on just one 
% well chosen mesh and plot the result

% Settings
dof = 3;                % dof = r+1 the polynomial degree (for P1 elements do dof = 2)
sigma = 10*dof^2;       % penalization parameter sigma 


% additional settings
resonators_matrix = [0.3, 0.7];
boundary_nodes = [0, 1];
h_background = 0.2;
h_res = h_background/5;

% functions
c_handle = @(x) 1*(x < resonators_matrix(1) | x > resonators_matrix(2)) + ...
                20*(x>= resonators_matrix(1) & x <= resonators_matrix(2));
f_handle = @(x) ones(size(x));

boundary_cond = struct("values", [0, 0], "lower_boundary_type", "dirichlet", "upper_boundary_type", "dirichlet");

% mesh
waveguide = mesh.MeshIntervalDG1d(boundary_nodes, [2*h_background,h_res/10]);
waveguide.dof = dof;
waveguide.buildResonatorMesh(resonators_matrix, [h_background, h_res]);
waveguide.updatePet();

% solve problem
quad_mesh = copy(waveguide);
quad_mesh.dof = quad_mesh.dof + 1;
quad_mesh.updatePet();
[quad_nodes, ~, quad_elements] = quad_mesh.getPet();
eval_nodes = quad_nodes(quad_elements);
eval_nodes(:, 1) = eval_nodes(:, 1) + 1e-14;
eval_nodes(:, end) = eval_nodes(:, end) - 1e-14;
c_vals = c_handle(eval_nodes);
f_vals = f_handle(eval_nodes);

% solve system
uh = dg1d.sip_1d_elliptic_solver(waveguide, boundary_cond, f_vals, c_vals, sigma);

% plot solution
f = figure;
waveguide.plotDGsol(uh, f);
grid on;