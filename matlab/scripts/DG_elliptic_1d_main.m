% Solves simple elliptic problem in 1d using DG and plots errors
clc;clear;close all;

% Imports
import mesh.*

% initialize mesh
Mesh = Mesh1d();

% define function handles (real solution)   
syms x
u_exact_handle = sin(x);
f_exact_handle = diff(u_exact_handle, 2);
u_exact_handle = matlabFunction(u_exact_handle);
f_exact_handle = matlabFunction(f_exact_handle);
