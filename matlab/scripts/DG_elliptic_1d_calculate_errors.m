% script to calculate and visualize errors of leapfrog sipg solution of the wave equation with time
% dependent coefficients
clc; clear; close all;

plot_errors = true;

% Settings
c_index = 1;
u_exact_handle_idx = 3;
is_resonator = true;
dof = 5;
num_ref = 6;
initial_meshsize = 0.25;
refine_factor = 2;
sigma = 10*dof^2;
boundary_nodes = [0, 1];

% functions
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
    x*0 + 1;
    sin(10*x) + 2;
    x.^(dof-1) + 1
};
c_handle = cell_c_fun{c_index};
u_exact_handle = cell_exact_fun{u_exact_handle_idx};
f_exact_handle = diff(c_handle*diff(-u_exact_handle, 1), 1) + u_exact_handle;
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
    waveguide.buildResonatorMesh([0.5, 0.6], [initial_meshsize, initial_meshsize/5]);
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
if plot_errors
    dof_plot = dof;
    line_width = 2;
    figure;
    loglog(meshsizes, meshsizes.^(dof_plot-1), '--','LineWidth', line_width);
    hold on
    loglog(meshsizes, meshsizes.^(dof_plot), '--', 'LineWidth', line_width);
    loglog(meshsizes, errors(1,:), 'LineWidth', line_width);
    loglog(meshsizes, errors(2,:), 'LineWidth', line_width);
    loglog(meshsizes, errors(3,:), 'LineWidth', line_width);
    hold off
    xlabel('Step Size (H)');
    ylabel('Error');
    legend("h^"+(dof_plot-1), "h^"+(dof_plot), "L2", "H1", "energy")
    title("Convergence of Errors for P^"+(dof-1)+ "elements");
    cr = fem1d.Errors1D.interpolate_convergence_rates(meshsizes, errors);
    fprintf('Interpolated Convergence Rates ------------- \nl2: %g,  h1: %g,  energy: %g \n\n', cr(2,1),cr(2,2),cr(2,3));
end
