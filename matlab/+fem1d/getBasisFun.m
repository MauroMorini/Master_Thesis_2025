function [phi_cell, dphi_cell] = getBasisFun(dof)
    % depending on degrees of freedom provided collects a cell array with function handles 
    % of basis functions defined on the reference element K = (-1,1)
    % 
    % Inputs:
    %       dof:        scalar degrees of freedom parameter
    %       
    % Output:   
    %       phi_cell:   (1,dof) cell array containing function handles of basis functions 
    %       dphi_cell:  (1,dof) cell array containing function handles of derivatives of basis functions
    switch dof
        case 2
            phi_cell = {@(x) 1-x, @(x) x};
            dphi_cell = {@(x) -ones(size(x)), @(x) ones(size(x))};
        case 3
            phi_cell = {@(x) 2*x.^2-3*x+1, @(x) -4*x.^2+4*x, @(x) 2*x.^2-x};
            dphi_cell = {@(x) 4*x-3, @(x) -8*x+4, @(x) 4*x-1};
        otherwise     
    end
end