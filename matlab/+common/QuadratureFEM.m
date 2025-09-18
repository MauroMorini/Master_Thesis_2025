classdef QuadratureFEM
    % This class encapsulates the quadrature functionality as well as 
    % the basis functions for fem
    methods (Static)
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
        function [quad_nodes, quad_weights] = getLobattoQuadrature(dof)
            % depending on dof parameter yields quadrature points and weights, currently only yields 
            % uniformly distributed nodal basis
            %
            % Inputs:
            %       dof:            scalar degrees of freedom
            % Outputs:
            %       quad_nodes:     (1, dof) node vector
            %       quad_weights:   (1, dof) weight vector
            switch dof
                case 2
                    quad_nodes = [-1, 1];
                    quad_weights = [1, 1];
                case 3
                    quad_nodes = [-1, 0, 1];
                    quad_weights = [1, 4, 1]/3;    
                case 4
                    quad_nodes = [-1, -1/sqrt(5), 1/sqrt(5), 1];
                    quad_weights = [1, 5, 5, 1]/6;            
                otherwise
                    error("method has not been implemented for dof = " + dof)
            end
        end
        function [phi_val, dphi_val, quad_weights] = getShapeFunctionValueMatrix(dof)
            % calculates matrix which contains the values of the shape functions
            % and their derivatives at the quadrature nodes
            %
            % Inputs:      
            %       dof:            scalar degree of freedom
            % Outputs:
            %       phi_val:        (dof, dof) matrix with entry (i,j) phi_i(x_j)
            %       dphi_val:       (dof, dof) matrix with entry (i,j) dphi_i(x_j)
            %       quad_weights:   (1, dof) weight vector
            phi_val = zeros(dof,dof);
            dphi_val = zeros(dof, dof);
            [quad_nodes, quad_weights] = common.QuadratureFEM.getLobattoQuadrature(dof);
            [phi_cell, dphi_cell] = common.QuadratureFEM.getLagrangeBasisFun(dof);
            for i = 1:dof
                phi_val(i,:) = phi_cell{i}(quad_nodes);
                dphi_val(i,:) = dphi_cell{i}(quad_nodes);
            end
        end
    end
end