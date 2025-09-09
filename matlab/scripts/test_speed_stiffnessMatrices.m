% test different stiffnessMatrix algos against each other to compare speeds
clc;clear;close all;

% Imports
import mesh.*
import fem1d.*

num_nodes_list = [200, 2000, 10000, 20000, 40000];
stepsizes = 1./(num_nodes_list + 1);
times_mat = zeros(2, length(num_nodes_list));
c = @(x) ones(size(x));

for i = 1:length(num_nodes_list)
    % create mesh
    h = stepsizes(i);
    Mesh = Mesh1d([0,1], [h, h/100]);
    [nodes, boundary_nodes_idx, elements] = Mesh.getPet();

    % standard matlab stiffness
    tic;
    A_matlab_stand = fem1d.stiffnessMatrix1D(nodes, elements, c);
    times_mat(1, i) = toc;

    % standard matlab stiffness
    tic;
    A_matlab_triplet = fem1d.stiffnessMatrix1D_triplets(nodes, elements, c);
    times_mat(2, i) = toc;
end

% plot results
figure;
plot(num_nodes_list, times_mat(1,:), num_nodes_list, times_mat(2,:))
title("comparison time needed to assemble stiffness matrices")
legend("standard matlab", "matlab with triplets")
xlabel("number of nodes")
ylabel("time required")