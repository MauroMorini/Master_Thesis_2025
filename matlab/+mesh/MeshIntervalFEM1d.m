classdef MeshIntervalFEM1d < handle
    % MeshIntervalFEM1d is a class for creating and managing a 1D mesh
    % for continuous Galerkin methods.
    %
    % Properties:
    %   nodes - Array of node coordinates in the mesh.
    %   elements - Connectivity matrix defining the mesh elements.
    %   boundary_interface_node_idx - Indices of boundary nodes in the mesh.
    %   dof - Number of degrees of freedom per element.
    %   h_max - Maximum mesh size.
    %   h_min - Minimum mesh size.
    %   lower_interval_bound - Lower bound of the mesh interval.
    %   upper_interval_bound - Upper bound of the mesh interval.
    %   element_interface_nodes - Coordinates of the element interface nodes.
    %   face_node_to_element_map - Mapping of internal face nodes to elements.
    %
    % Methods:
    %   MeshIntervalDG1d - Constructor for the MeshIntervalDG1d class.
    %   updatePet - Updates the mesh properties and connectivity based on the
    %               current element interface nodes.
    %   getPet - Returns the current nodes, elements, and boundary node indices.
    %   plotMesh - Plots the mesh and highlights the boundary nodes.
    properties
        nodes
        elements
        boundary_interface_node_idx
        dof
        h_max
        h_min
        lower_interval_bound
        upper_interval_bound
        element_interface_nodes
        face_node_to_element_map
    end

    methods
        function obj = MeshIntervalFEM1d(varargin)
            switch (nargin)
                case 0
                    % default: equidistant mesh on [0,1] with mesh size 0.1
                    obj.lower_interval_bound = 0;
                    obj.upper_interval_bound = 1;
                    obj.h_max = 0.1;
                    obj.h_min = obj.h_max/100;    
                    obj.dof = 2;                
                case 2
                    % inputs: ([lower_interval_bound,upper_interval_bound],[h_max, h_min])
                    intBounds = varargin{1};
                    Hbounds = varargin{2};
                    obj.lower_interval_bound = intBounds(1);
                    obj.upper_interval_bound = intBounds(2);
                    obj.h_max = Hbounds(1);
                    obj.h_min = Hbounds(2);
                    obj.dof = 2;
                otherwise
                    error("Constructor has wrong number of inputs")
            end
            assert(obj.h_max >= 5*obj.h_min, "h_min should be small enough, maximally h_max/5")
            obj.element_interface_nodes = (obj.lower_interval_bound:obj.h_max:(obj.upper_interval_bound-obj.h_min))'; 
            obj.element_interface_nodes(end+1) = obj.upper_interval_bound;
            while abs(obj.element_interface_nodes(end) - obj.element_interface_nodes(end-1)) > obj.h_max
                obj.element_interface_nodes = [obj.element_interface_nodes(1:end-1); (obj.element_interface_nodes(end-1)+obj.upper_interval_bound)/2; obj.upper_interval_bound];
            end
            obj = updatePet(obj);
        end

        function obj = updatePet(obj)
            % updates connectivity matrix and edge matrix 
            obj.element_interface_nodes = sort(obj.element_interface_nodes);        % possibly remove sorting
            if obj.element_interface_nodes(1) ~= obj.lower_interval_bound || obj.element_interface_nodes(end) ~= obj.upper_interval_bound
                error("boundary points of element_interface_nodes don't coincide with properties of class")
            end
            % check that Hmin and Hmax are maintained
            D = abs(diff(obj.element_interface_nodes));
            corrSize = all(obj.h_min-10*eps <= D & D <= obj.h_max+10*eps);
            if ~corrSize
                warning("There are elements which are too small or too big")
            end

            obj.setNodes();

            N = length(obj.element_interface_nodes);
            % connectivity matrix
            elements_temp = zeros(N-1, obj.dof);
            for i = 1:obj.dof
                elements_temp(:,i) = (i:(obj.dof-1):((N-2)*(obj.dof-1)+i))';
            end
            obj.elements = elements_temp;

            % face_node_to_element_map
            obj.face_node_to_element_map = [NaN, 1; (1:N-2).', (2:N-1).'; N-1, NaN];
        end

        function obj = setNodes(obj)
            % creates interior nodes based on Gauss-Lobatto
            % and updates object
            N = length(obj.element_interface_nodes);
            [quad_nodes, ~] = obj.getLobattoQuadrature(obj.dof);
            element_meshsizes = diff(obj.element_interface_nodes);
            element_midpoints = obj.element_interface_nodes(1:end-1) + element_meshsizes/2;
            nodes_matrix_extended = zeros(N-1, obj.dof-1);
            nodes_matrix_extended(:, 1) = obj.element_interface_nodes(1:end-1);
            nodes_matrix_extended(:, 2:end) = element_midpoints + quad_nodes(2:end-1).*element_meshsizes/2;
            nodes_matrix_extended = nodes_matrix_extended.';
            obj.nodes = [nodes_matrix_extended(:);obj.upper_interval_bound];
            num_nodes = length(obj.nodes);
            obj.boundary_interface_node_idx = [1, num_nodes];
            
        end
        function [nodes, boundary_interface_node_idx, elements] = getPet(obj)
            nodes = obj.nodes;
            elements = obj.elements;
            boundary_interface_node_idx = obj.boundary_interface_node_idx;
            return 
        end

        function f = plotMesh(obj, f)
            if nargin == 1
                f = figure;
            end
            figure(f);            
            plot(obj.element_interface_nodes,0,'b.','MarkerSize',10)
            hold on
            plot(obj.element_interface_nodes(obj.boundary_interface_node_idx), [0,0], 'rx','MarkerSize',10)
            hold off
        end
    end
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
    end
end