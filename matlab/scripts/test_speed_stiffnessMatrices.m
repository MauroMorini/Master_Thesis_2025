% test different stiffnessMatrix algos against each other to compare speeds
clc;clear;close all;

% Imports
import mesh.*
import fem1d.*
import fem1d_old_versions.*

num_nodes_list = [200, 2000, 10000, 20000, 40000];
stepsizes = 1./(num_nodes_list + 1);
times_mat = zeros(3, length(num_nodes_list));
c = @(x) ones(size(x));


for i = 1:length(num_nodes_list)
    % create mesh
    h = stepsizes(i);
    Mesh = mesh.Mesh1dBroken([0,1], [h, h/100]);
    [nodes, boundary_nodes_idx, elements] = Mesh.getPet();

    % standard matlab stiffness original
    tic;
    A_matlab_stand = fem1d_old_versions.stiffnessMatrix1D_v0(nodes, elements, c);
    times_mat(1, i) = toc;

    % matlab stiffness triplets with function calls
    tic;
    A_matlab_triplet = fem1d_old_versions.stiffnessMatrix1D_v1(nodes, elements, c);
    times_mat(2, i) = toc;

    % matlab stiffness triplets without function calls
    tic;
    c_vals = c(nodes(elements));
    A_matlab_triplet_no_fun = fem1d.stiffnessMatrix1D(nodes, elements, c_vals);
    times_mat(3, i) = toc;
end

% plot results  
figure;
plot(num_nodes_list, times_mat(1,:), num_nodes_list, times_mat(2,:),num_nodes_list, times_mat(3,:))
title("comparison time needed to assemble stiffness matrices")
legend("matlab v0", "matlab v1", "matlab v2")
xlabel("number of nodes")
ylabel("time required")
