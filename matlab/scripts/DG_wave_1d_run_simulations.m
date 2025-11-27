% script to simulate and observe a wave solved with sipg in space and leapfrog in time

clc;clear;close all;
save_to_h5 = false;
plot_images = true;

% Settings
h = 0.2;
filename_index = 310;
u_exact_index = 4;
wave_speed_index = 2;
dof = 3;
dt_scaling = 0.5/10;
res_mat = [5, 6];
plot_times = 0:5;

filepath = "data/matlab/wave/hdf5/";
filename = filepath + "wave-" + filename_index + ".h5";
[pde_data, res_mat] = fem1d.PDEData.generate_pde_data_with_resonators(u_exact_index, wave_speed_index, res_mat);
% pde_data = fem1d.PDEData.generate_smooth_pde_data(u_exact_index, wave_speed_index);
waveguide = mesh.MeshIntervalDG1d(pde_data.boundary_points, [2*h, h/50]);
waveguide.dof = dof;
% waveguide.createUniformMesh(h);
waveguide.buildResonatorMesh(res_mat, [h, h/5]);
waveguide.updatePet();
sipg_solver = dg1d.SIPGWaveSolver1D(waveguide, pde_data);
% sipg_solver.dt = waveguide.h_min*dt_scaling/waveguide.dof;
% sipg_solver.sigma = 100;
% sipg_solver.matrix_update_type = "piecewise-const-coefficient-in-space";
% sipg_solver.matrix_update_type = "brute-force";
sipg_solver.run();

wave_postprocessor = dg1d.WavePostprocessor1D(sipg_solver);

if save_to_h5
    wave_postprocessor.write_to_hdf5(filename);
end

if plot_images
    wave_postprocessor.plot_solutions(plot_times);
    drawnow;
end