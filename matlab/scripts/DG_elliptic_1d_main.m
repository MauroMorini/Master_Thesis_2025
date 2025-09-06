% Solves simple elliptic problem in 1d using DG and plots errors
clc;clear;close all;

% Imports
import mesh.*
import fem1d.*

% initialize mesh
Mesh = Mesh1d([0,1], [0.2, 0.002]);
[nodes, boundary_nodes_idx, elements] = Mesh.getPet();

% define function handles (real solution)   
syms x
u_exact_handle = sin(x);
f_exact_handle = diff(u_exact_handle, 2);
u_exact_handle = matlabFunction(u_exact_handle);
f_exact_handle = matlabFunction(f_exact_handle);
c_handle = @(x) ones(size(x));

% assemble matrices
A = fem1d.stiffnessMatrix1D(nodes', elements, c_handle);
load_v = fem1d.loadVectorLinear1D(nodes', elements, f_exact_handle);

