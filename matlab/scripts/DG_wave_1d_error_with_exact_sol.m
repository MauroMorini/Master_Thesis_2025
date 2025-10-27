% script to calculate and visualize errors of leapfrog sipg solution of the wave equation with time
% dependent coefficients
clc; clear; close all;

% Settings
h = 0.1;
dof = 2;
plot_time = 2;
num_plot_nodes = 1000;

pde_data = fem1d.PDEData.generate_gaussian_puls_data_on_waveguide();
waveguide = mesh.MeshIntervalDG1d(pde_data.boundary_points, [2*h, h/50]);
waveguide.createRngMesh();
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

