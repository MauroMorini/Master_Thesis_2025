classdef MeshInterval1d < handle

    properties
        nodes
        elements
        boundary_element_idx
        dof
        h_max
        h_min
        lower_interval_bound
        upper_interval_bound
        element_boundary_nodes
    end

    methods
        function obj = MeshInterval1d(varargin)
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
            obj.element_boundary_nodes = (obj.lower_interval_bound:obj.h_max:(obj.upper_interval_bound-obj.h_min))'; 
            obj.element_boundary_nodes(end+1) = obj.upper_interval_bound;
            while abs(obj.element_boundary_nodes(end) - obj.element_boundary_nodes(end-1)) > obj.h_max
                obj.element_boundary_nodes = [obj.element_boundary_nodes(1:end-1); (obj.element_boundary_nodes(end-1)+obj.upper_interval_bound)/2; obj.upper_interval_bound];
            end
            obj = updatePet(obj);
        end

        function obj = updatePet(obj)
            % updates connectivity matrix and edge matrix 
            obj.element_boundary_nodes = sort(obj.element_boundary_nodes);
            if obj.element_boundary_nodes(1) ~= obj.lower_interval_bound || obj.element_boundary_nodes(end) ~= obj.upper_interval_bound
                error("boundary points of element_boundary_nodes don't coincide with properties of class")
            end
            % check that Hmin and Hmax are maintained
            D = abs(diff(obj.element_boundary_nodes));
            corrSize = all(obj.h_min-10*eps <= D & D <= obj.h_max+10*eps);
            if ~corrSize
                warning("There are elements which are too small or too big")
            end
            N = length(obj.element_boundary_nodes);
            obj.boundary_element_idx = [1, N];

            % add additional nodes for higher dofs to obj.nodes
            p = obj.element_boundary_nodes;
            H = diff(p)/(obj.dof-1);
            for i = 1:(obj.dof-2)
                local_nodes = obj.element_boundary_nodes(1:end-1) + H*i;
                p = [p;local_nodes];
            end
            obj.nodes = sort(p);
            
            obj.elements = zeros(size(obj.element_boundary_nodes-1,1), obj.dof);
            % connectivity matrix
            for i = 1:obj.dof
                obj.elements(:,i) = (i:(obj.dof-1):(length(p)-(obj.dof-i)))';
            end
        end

        function [nodes, elements, boundary_element_idx] = getPet(obj)
            nodes = obj.nodes;
            elements = obj.elements;
            boundary_element_idx = obj.boundary_element_idx;
            return 
        end

        function f = plotMesh(obj, f)
            if nargin == 1
                f = figure;
            end
            figure(f);            
            plot(obj.element_boundary_nodes,0,'b.','MarkerSize',10)
            hold on
            plot(obj.element_boundary_nodes(obj.boundary_element_idx), [0,0], 'rx','MarkerSize',10)
            hold off
        end
    end
end