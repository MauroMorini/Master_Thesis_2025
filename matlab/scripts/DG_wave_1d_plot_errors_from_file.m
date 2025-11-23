%  extracts csv file containing errors and plots them for visualization
% clc;clear;close all;

% Settings
filename_index = 2;
dof = 3;

% read csv
file_path = fullfile("data", "matlab", "wave", "csv","wave_error-");
filename_temp = filename_index + ".csv";
filename = file_path + filename_temp;
[meshsizes, errors, metadata] = fem1d.Errors1D.read_errors_from_csv(filename);

% plot 
disp(metadata);
line_width = 2;
dof_plot = dof;
figure;
loglog(meshsizes, meshsizes.^(dof_plot-1), '--','LineWidth', line_width);
hold on
grid on;
loglog(meshsizes, meshsizes.^(dof_plot), '--', 'LineWidth', line_width);
loglog(meshsizes, errors(1,:),'^-' ,'LineWidth', line_width);
loglog(meshsizes, errors(2,:),'o-' ,'LineWidth', line_width);
loglog(meshsizes, errors(3,:), 'x-','LineWidth', line_width);
hold off
grid on;
xlabel('meshsize h');
ylabel('Error');
legend("h^"+(dof_plot-1), "h^"+(dof_plot), "L2", "H1", "energy")
title("Errors of the SIPG wave solution");