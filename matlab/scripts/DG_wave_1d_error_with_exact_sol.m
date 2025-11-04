% script to calculate and visualize errors of leapfrog sipg solution of the wave equation with time
% dependent coefficients
clc; clear; close all;

% Settings
initial_meshsize = 1;
dof = 3;
num_ref = 5;
refine_factor = 2;
c_index = 2;
dt_scaling_factor = 0.5/10;

% initialization
errors = zeros(2, num_ref);
meshsizes = zeros(1, num_ref);
pde_data = fem1d.PDEData.generate_gaussian_puls_data_on_waveguide(c_index);
waveguide = mesh.MeshIntervalDG1d(pde_data.boundary_points, [2*initial_meshsize, initial_meshsize/50]);
waveguide.createUniformMesh(initial_meshsize);
waveguide.buildResonatorMesh([4, 4.5; 8.5, 9], [initial_meshsize, initial_meshsize/5]);
waveguide.dof = dof;
waveguide.updatePet();

for i = 1:num_ref
    if i > 1
        waveguide.refineAll(2);
        waveguide.updatePet();
    end
    sipg_solver = dg1d.SIPGWaveSolver1D(waveguide, pde_data);
    sipg_solver.matrix_update_type = "piecewise-const-coefficient-in-space";
    % sipg_solver.dt = waveguide.h_min*dt_scaling_factor/dof;
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

