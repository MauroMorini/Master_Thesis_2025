% script to simulate and observe a wave solved with sipg in space and leapfrog in time

clc;clear;close all;

% Settings
% h = 0.1;
% dof = 2;
% plot_times = 2;
% plot_times = 0:0.3:2.5;
% num_plot_nodes = 1000;
% c_index = 3;
% dt_scaling = 0.5/10;

% pde_data = fem1d.PDEData.generate_gaussian_puls_data_on_waveguide(c_index);
% waveguide = mesh.MeshIntervalDG1d(pde_data.boundary_points, [2*h, h/50]);
% waveguide.dof = dof;
% waveguide.buildResonatorMesh([4, 6; 8, 9], [h, h/5]);
% waveguide.updatePet();
% sipg_solver = dg1d.SIPGWaveSolver1D(waveguide, pde_data);
% sipg_solver.dt = waveguide.h_min*dt_scaling/waveguide.dof;
% % sipg_solver.sigma = 20;
% %c sipg_solver.matrix_update_type = "brute-force";
% sipg_solver.run();

%% 
% Settings
h = 0.3;
dof = 3;
dt_scaling = 0.5/10;
c_idx = 3;

[pde_data, res_mat] = fem1d.PDEData.generate_gaussian_puls_data_on_waveguide_with_resonators(c_idx);
waveguide = mesh.MeshIntervalDG1d(pde_data.boundary_points, [2*h, h/50]);
waveguide.dof = dof;
waveguide.buildResonatorMesh(res_mat, [h, h/5]);
waveguide.updatePet();
sipg_solver = dg1d.SIPGWaveSolver1D(waveguide, pde_data);
sipg_solver.dt = waveguide.h_min*dt_scaling/waveguide.dof;
% sipg_solver.sigma = 20;
% sipg_solver.matrix_update_type = "brute-force";
sipg_solver.run();

%% plot
plot_times = 0:5;
num_plot_nodes = 1000;
for plot_time = plot_times
    wave_postprocessor = dg1d.WavePostprocessor1D(sipg_solver);
    [uh, t] = wave_postprocessor.get_solution_at_time(plot_time);
    f = figure;
    hold on
    % plot_nodes = linspace(pde_data.boundary_points(1),pde_data.boundary_points(2),num_plot_nodes)';
    % plot(plot_nodes, pde_data.u_exact_fun(plot_nodes, t), 'LineWidth', 3)
    sipg_solver.initial_mesh.plotDGsol(uh, f);
    hold off
end

%% plot animation
f = figure;
plot_times = 0:0.1:20;
num_plot_nodes = 1000;
wave_postprocessor = dg1d.WavePostprocessor1D(sipg_solver);
F(length(plot_times)) = struct('cdata',[],'colormap',[]);

for i = 1:length(plot_times)
    [uh, t] = wave_postprocessor.get_solution_at_time(plot_times(i));
    f_exact = plot(plot_nodes, sipg_solver.pde_data.u_exact_fun(plot_nodes, t), 'LineWidth', 3);
    f_num = sipg_solver.initial_mesh.plotDGsol(uh, f);
    drawnow limitrate;
    F(i) = getframe(f);
    pause(0.1)
    cla;
end

v = VideoWriter('wave_simulation.avi');
v.FrameRate = 30;
open(v);
writeVideo(v, F);
close(v);
disp('Video saved as wave_simulation.avi');