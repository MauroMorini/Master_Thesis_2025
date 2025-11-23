% This is the main script of chapter 2 numerical results. It is subdivided
% into sections corresponding to certain numerical experiments presented in chapter 2.
% The sections are self-contained and can be ran individually, therefore
% they also share a lot of repeating code.

clc; clear; close all;
%% Convergence Rates
% script to calculate and visualize errors of leapfrog sipg solution of the wave equation with time
% dependent coefficients on a uniform mesh.
% If the computation takes too long consider choosing a smaller final time or reducing the max 
% refinement cycles "num_ref".

% Settings
c_index = 1;            % coefficient index choose int 1-4
u_exact_index = 2;      % exact solution index choose int 1-2
dof = 3;                % number of basis nodes per element, r = dof-1, P1 elements are dof = 2 


% additional Settings
final_time = 10;        % T = 10 endtime
is_resonator = false;
dt_scaling_factor = 0.5/2;
num_ref = 9;
initial_meshsize = 2;
refine_factor = 2;


% initialization
errors = zeros(2, num_ref);
meshsizes = zeros(1, num_ref);
pde_data = fem1d.PDEData.generate_smooth_pde_data(u_exact_index, c_index);
pde_data.final_time = final_time;
waveguide = mesh.MeshIntervalDG1d(pde_data.boundary_points, [2*initial_meshsize, initial_meshsize/50]);
waveguide.createUniformMesh(initial_meshsize);
if is_resonator
    waveguide.buildResonatorMesh([4, 4.5; 8.5, 9], [initial_meshsize, initial_meshsize/5]);
end
waveguide.dof = dof;
waveguide.updatePet();

for i = 1:num_ref
    if i > 1
        waveguide.refineAll(refine_factor);
        waveguide.updatePet();
    end
    sipg_solver = dg1d.SIPGWaveSolver1D(waveguide, pde_data);
    if ~is_resonator
        sipg_solver.dt = waveguide.h_min*dt_scaling_factor/(dof-1);
    end
    sipg_solver.run();
    wave_postprocessor = dg1d.WavePostprocessor1D(sipg_solver);
    wave_postprocessor.calculate_errors();
    [errors(1,i), errors(2,i), errors(3,i)] = wave_postprocessor.errors_obj.getErrors();
    meshsizes(i) = waveguide.h_max;
end

dof_plot = dof;
line_width = 2;
figure;
loglog(meshsizes, meshsizes.^(dof_plot-1), '--','LineWidth', line_width);
hold on
loglog(meshsizes, meshsizes.^(dof_plot), '--', 'LineWidth', line_width);
loglog(meshsizes, errors(1,:),'^-' ,'LineWidth', line_width);
loglog(meshsizes, errors(2,:),'o-' ,'LineWidth', line_width);
loglog(meshsizes, errors(3,:), 'x-','LineWidth', line_width);
hold off
xlabel('meshsize h');
ylabel('error');
legend("h^"+(dof_plot-1), "h^"+(dof_plot), "L2", "H1", "energy")
title("Error rates for P^"+(dof-1)+ "elements");
grid on;
cr = fem1d.Errors1D.interpolate_convergence_rates(meshsizes, errors);
fprintf('Interpolated Convergence Rates ------------- \nl2: %g,  h1: %g,  energy: %g \n\n', cr(2,1),cr(2,2),cr(2,3));

