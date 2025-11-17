% script to calculate and visualize errors of leapfrog sipg solution of the wave equation with time
% dependent coefficients
clc; clear; close all;

save_csv = true;
plot_errors = true;

% Settings
filename_index = 201;
c_index = 2;
u_exact_index = 1;
is_resonator = false;
dof = 2;
dt_scaling_factor = 0.5/10;
num_ref = 6;
initial_meshsize = 1;
refine_factor = 2;
tic 

% initialization
errors = zeros(2, num_ref);
meshsizes = zeros(1, num_ref);
pde_data = fem1d.PDEData.generate_smooth_pde_data(u_exact_index, c_index);
waveguide = mesh.MeshIntervalDG1d(pde_data.boundary_points, [2*initial_meshsize, initial_meshsize/50]);
waveguide.createUniformMesh(initial_meshsize);
if is_resonator
    waveguide.buildResonatorMesh([4, 4.5; 8.5, 9], [initial_meshsize, initial_meshsize/5]);
end
waveguide.dof = dof;
waveguide.updatePet();

for i = 1:num_ref
    if i > 1
        waveguide.refineAll(2);
        waveguide.updatePet();
    end
    sipg_solver = dg1d.SIPGWaveSolver1D(waveguide, pde_data);
    % sipg_solver.dt = waveguide.h_min*dt_scaling_factor/dof;
    sipg_solver.run();
    wave_postprocessor = dg1d.WavePostprocessor1D(sipg_solver);
    wave_postprocessor.calculate_errors();
    [errors(1,i), errors(2,i), errors(3,i)] = wave_postprocessor.errors_obj.getErrors();
    meshsizes(i) = waveguide.h_max;
end
time_elapsed = toc;

% save to csv
if save_csv
    left_bc = pde_data.boundary_conditions{1}.bc_type;
    right_bc = pde_data.boundary_conditions{2}.bc_type;
    syms x t
    exact_sol = string(pde_data.generateSmoothExactSol(u_exact_index));
    wave_speed = string(pde_data.generateSmoothWaveSpeed(c_index));
    metadata = sprintf(['# Boundary Conditions: [%s, %s] \n'...
                        '# Exact Solution: %s \n' ...
                        '# Wave Speed: %s \n'...
                        '# Time-interval: [%d, %d]\n'...
                        '# Computation time: %d \n'],...
                        left_bc,right_bc, exact_sol, wave_speed,...
                        pde_data.initial_time, pde_data.final_time, time_elapsed);
    filename_temp = "wave_error-" + string(filename_index) + ".csv";
    filename = fullfile('data','matlab', 'wave', 'csv', filename_temp);
    wave_postprocessor.errors_obj.write_errors_to_csv(filename, errors, meshsizes, metadata);
end 

% plot errors
if plot_errors
    dof_plot = 2;
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
    cr = wave_postprocessor.errors_obj.interpolate_convergence_rates(meshsizes, errors);
    fprintf('Interpolated Convergence Rates ------------- \nl2: %g,  h1: %g,  energy: %g \n\n', cr(2,1),cr(2,2),cr(2,3));
end
