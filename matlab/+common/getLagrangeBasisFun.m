function [phi_cell, dphi_cell] = getLagrangeBasisFun(dof)
    % depending on degrees of freedom provided collects a cell array with function handles 
    % of basis Lagrange functions defined on the reference element K = (-1,1)
    % 
    % Inputs:
    %       dof:        scalar degrees of freedom parameter
    %       
    % Output:   
    %       phi_cell:   (1,dof) cell array containing function handles of basis functions 
    %       dphi_cell:  (1,dof) cell array containing function handles of derivatives of basis functions
    switch dof
        case 2
            phi_cell = {@(xi) (1-xi)/2, @(xi) (1+xi)/2};
            dphi_cell = {@(xi) -(1/2)*ones(size(xi)), @(xi) (1/2)*ones(size(xi))};
        case 3
            phi_cell = {@(xi)1/2*(xi.^2 - xi), @(xi) 1 - xi.^2, @(xi) 1/2*(xi + xi.^2)};
            dphi_cell = {@(xi)1/2*(2*xi-1), @(xi) -2*xi, @(xi) 1/2*(1+2*xi)};
        otherwise     
    end
end