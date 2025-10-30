% script to calculate and visualize errors of leapfrog sipg solution of the wave equation with time
% dependent coefficients
clc; clear; close all;

% Settings
h = 1;
dof = 3;
plot_time = 2;
num_plot_nodes = 1000;

pde_data = fem1d.PDEData.generate_gaussian_puls_data_on_waveguide();
waveguide = mesh.MeshIntervalDG1d(pde_data.boundary_points, [2*h, h/50]);
waveguide.createUniformMesh(h);
waveguide.dof = dof;
waveguide.updatePet();
sipg_solver = dg1d.SIPGWaveSolver1D(waveguide, pde_data);
%sipg_solver.sigma = 20;
sipg_solver.run();

wave_postprocessor = dg1d.WavePostprocessor1D(sipg_solver);
[uh, t] = wave_postprocessor.get_solution_at_time(plot_time);
f = figure;
hold on
plot_nodes = linspace(pde_data.boundary_points(1),pde_data.boundary_points(2),num_plot_nodes)';
plot(plot_nodes, pde_data.u_exact_fun(plot_nodes, t), 'LineWidth', 3)
waveguide.plotDGsol(uh, f);
hold off



%% errors

% Settings
initial_meshsize = 2;
dof = 2;
num_ref = 6;
refine_factor = 2;
c_index = 2;
dt_scaling_factor = 0.5;

% initialization
errors = zeros(2, num_ref);
meshsizes = zeros(1, num_ref);
pde_data = fem1d.PDEData.generate_gaussian_puls_data_on_waveguide(c_index);
waveguide = mesh.MeshIntervalDG1d(pde_data.boundary_points, [2*initial_meshsize, initial_meshsize/50]);
waveguide.createUniformMesh(initial_meshsize);
waveguide.buildResonatorMesh([4, 6], [initial_meshsize, initial_meshsize/5]);
waveguide.dof = dof;
waveguide.updatePet();

for i = 1:num_ref
    if i > 1
        waveguide.refineAll(2);
        waveguide.updatePet();
    end
    sipg_solver = dg1d.SIPGWaveSolver1D(waveguide, pde_data);
    % sipg_solver.dt = waveguide.h_min*dt_scaling_factor;
    sipg_solver.run();
    wave_postprocessor = dg1d.WavePostprocessor1D(sipg_solver);
    wave_postprocessor.calculate_errors();
    [errors(1,i), errors(2,i), errors(3,i)] = wave_postprocessor.errors_obj.getErrors();
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
loglog(meshsizes, errors(3,:), 'LineWidth', line_width);
hold off
xlabel('Step Size (H)');
ylabel('Error');
dof_plot = 2;
legend("h^"+(dof_plot-1), "h^"+(dof_plot), "L2", "H1", "energy")
title("Convergence of Errors for P^"+(dof-1)+ "elements");

