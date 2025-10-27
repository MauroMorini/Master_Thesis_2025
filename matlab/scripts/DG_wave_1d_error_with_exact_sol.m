% script to calculate and visualize errors of leapfrog sipg solution of the wave equation with time
% dependent coefficients
clc; clear; close all;

% Settings
h = 1;
dof = 2;
plot_time = 2;
num_plot_nodes = 1000;

pde_data = fem1d.PDEData.generate_gaussian_puls_data_on_waveguide();
waveguide = mesh.MeshIntervalDG1d(pde_data.boundary_points, [2*h, h/50]);
waveguide.createUniformMesh(h);
waveguide.dof = dof;
waveguide.updatePet();
sipg_solver = dg1d.SIPGWaveSolver1D(waveguide, pde_data);
sipg_solver.run();

[~, time_idx] = min(abs(sipg_solver.times - plot_time));
uh = sipg_solver.solution(:,time_idx);
f = figure;
waveguide.plotDGsol(uh, f);

figure;
plot_nodes = linspace(pde_data.boundary_points(1),pde_data.boundary_points(2),num_plot_nodes)';
plot(plot_nodes, pde_data.u_exact_fun(plot_nodes, sipg_solver.times(time_idx)))

f = figure;
waveguide.plotMesh(f);

%% errors

% Settings
initial_meshsize = 2;
dof = 2;
num_ref = 8;
refine_factor = 2;

% initialization
errors = zeros(2, num_ref);
meshsizes = zeros(1, num_ref);
pde_data = fem1d.PDEData.generate_gaussian_puls_data_on_waveguide();
waveguide = mesh.MeshIntervalDG1d(pde_data.boundary_points, [2*initial_meshsize, initial_meshsize/50]);
waveguide.createUniformMesh(initial_meshsize);
waveguide.dof = dof;
waveguide.updatePet();

for i = 1:num_ref
    if i > 1
        waveguide.refineAll(2);
        waveguide.updatePet();
    end
    sipg_solver = dg1d.SIPGWaveSolver1D(waveguide, pde_data);
    sipg_solver.run();
    numerical_sol_struct = struct("mesh", waveguide, "sol", sipg_solver.solution(:, end), "type","numerical_solution");
    final_time_num = sipg_solver.times(end);
    exact_sol_struct = struct("u_handle", @(x) pde_data.u_exact_fun(x,final_time_num),...
                              "du_handle", @(x) pde_data.grad_u_exact_fun(x,final_time_num),"type","exact_solution");
    [errors(1,i), errors(2,i)] = fem1d.errors1D(numerical_sol_struct, exact_sol_struct);
    meshsizes(i) = waveguide.h_max;
end

% plot errors
line_width = 2;
figure;
loglog(meshsizes, meshsizes.^(dof-1), '--','LineWidth', line_width);
hold on
loglog(meshsizes, meshsizes.^(dof), '--', 'LineWidth', line_width);
loglog(meshsizes, errors(1,:), 'LineWidth', line_width);
loglog(meshsizes, errors(2,:), 'LineWidth', line_width);
% loglog(meshsizes, errors(3,:), 'LineWidth', line_width);
hold off
xlabel('Step Size (H)');
ylabel('Error');
legend("h^"+(dof-1), "h^"+(dof), "L2", "H1", "energy")
title("Convergence of Errors for P^"+(dof-1)+ "elements");

