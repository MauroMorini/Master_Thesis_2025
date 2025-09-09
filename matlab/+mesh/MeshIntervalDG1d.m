classdef MeshIntervalDG1d < handle

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
    end

    methods
        function obj = MeshIntervalDG1d(varargin)
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
            N = length(obj.element_interface_nodes);

            % add additional nodes for higher dofs to obj.nodes and add
            % double nodes for DG
            H = diff(obj.element_interface_nodes)/(obj.dof-1);
            nodes_matrix_extended = zeros(N-1, obj.dof);
            nodes_matrix_extended(:, 1) = obj.element_interface_nodes(1:end-1);
            for i = 1:(obj.dof-2)
                local_nodes = obj.element_interface_nodes(1:end-1) + H*i;
                nodes_matrix_extended(:, i+1) = local_nodes;
            end
            nodes_matrix_extended(:, end) = obj.element_interface_nodes(2:end);
            nodes_matrix_extended = nodes_matrix_extended.';
            obj.nodes = nodes_matrix_extended(:);
            obj.boundary_interface_node_idx = [1, length(obj.nodes)];
            
            % connectivity matrix
            elements_temp = 1:length(obj.nodes);
            elements_temp = reshape(elements_temp, obj.dof, []);
            obj.elements = elements_temp.';
        end

        function [nodes, elements, boundary_interface_node_idx] = getPet(obj)
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
end