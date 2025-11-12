% script to calculate and visualize errors of leapfrog sipg solution of the wave equation with time
% dependent coefficients
clc; clear; close all;

save_csv = true;
plot_errors = false;

% Settings
filename_index = 1;
c_index = 1;
is_resonator = true;
dof = 2;
dt_scaling_factor = 0.5/10;
num_ref = 5;
initial_meshsize = 1;
refine_factor = 2;

% metadata
exact_sol = "@(x,t) exp(-(x-t+1).^2)";

% initialization
errors = zeros(2, num_ref);
meshsizes = zeros(1, num_ref);
pde_data = fem1d.PDEData.generate_exact_gaussian_puls_data_on_waveguide(c_index);
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
    % sipg_solver.matrix_update_type = "piecewise-const-coefficient-in-space";
    % sipg_solver.dt = waveguide.h_min*dt_scaling_factor/dof;
    sipg_solver.run();
    wave_postprocessor = dg1d.WavePostprocessor1D(sipg_solver);
    wave_postprocessor.calculate_errors();
    [errors(1,i), errors(2,i), errors(3,i)] = wave_postprocessor.errors_obj.getErrors();
    meshsizes(i) = waveguide.h_max;
end

% save to csv
if save_csv
    left_bc = pde_data.boundary_conditions{1}.bc_type;
    right_bc = pde_data.boundary_conditions{2}.bc_type;
    syms x t
    wave_speed = string(func2str(matlabFunction(pde_data.generateSmoothWaveSpeed(c_index), 'Vars', {x,t})));
    metadata = sprintf(['# Boundary Conditions: [%s, %s] \n'...
                        '# Exact Solution: %s \n' ...
                        '# Wave Speed: %s \n'...
                        '# Time-interval: [%d, %d]\n'],...
                        left_bc,right_bc, exact_sol, wave_speed,...
                        pde_data.initial_time, pde_data.final_time);
    filename_temp = "wave_error-" + filename_index + ".csv";
    filename = fullfile('data','matlab', 'wave', 'csv', filename_temp);
    wave_postprocessor.errors_obj.write_errors_to_csv(filename, errors, meshsizes, metadata);
end 

% plot errors
if plot_errors
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
end
