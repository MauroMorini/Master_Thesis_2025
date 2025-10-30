% script to simulate and observe a wave solved with sipg in space and leapfrog in time

clc;clear;close all;

% Settings
h = 1;
dof = 2;
plot_times = 0:0.3:2.5;
num_plot_nodes = 1000;
c_index = 3;

pde_data = fem1d.PDEData.generate_gaussian_puls_data_on_waveguide(c_index);
waveguide = mesh.MeshIntervalDG1d(pde_data.boundary_points, [2*h, h/50]);
waveguide.dof = dof;
waveguide.buildResonatorMesh([4, 6; 8, 9], [h, h/5]);
waveguide.updatePet();
sipg_solver = dg1d.SIPGWaveSolver1D(waveguide, pde_data);
%sipg_solver.sigma = 20;
sipg_solver.matrix_update_type = "piecewise-const-coefficient-in-space";
sipg_solver.run();

for plot_time = plot_times
    wave_postprocessor = dg1d.WavePostprocessor1D(sipg_solver);
    [uh, t] = wave_postprocessor.get_solution_at_time(plot_time);
    f = figure;
    hold on
    plot_nodes = linspace(pde_data.boundary_points(1),pde_data.boundary_points(2),num_plot_nodes)';
    plot(plot_nodes, pde_data.u_exact_fun(plot_nodes, t), 'LineWidth', 3)
    waveguide.plotDGsol(uh, f);
    hold off
end