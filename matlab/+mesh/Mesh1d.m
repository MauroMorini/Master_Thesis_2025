classdef Mesh1d < handle
    % Mesh object for 1d FEM/DG methods, currently only supports connected intervals 
    %
    % properties:
    %   nodes:                  (N,1) matrix of node values contained in mesh
    %   connectivity_matrix:    (N,dof) index matrix of nodes, each line corresponding to one element
    %   boundary_nodes_idx:     (2,1) index matrix of boundary values
    %   dof:                    int, global degree of freedom per element (for now all elements have the same dof)
    %   number_of_nodes         int, corresponds to length(nodes) = N
    %   
    % methods:
    %   Mesh1d:                 default constructor, possible input 
    %   
    %   

    properties
        nodes
        connectivity_matrix 
        boundary_nodes_idx
        dof
        number_of_nodes
    end

    methods
        function obj = Mesh1d(stepsize, lower_bound, upper_bound, dof)
            % Default constructor, has 3 mode:
            %   mode 1: no inputs
            %   mode 2: only input stepsize
            %   mode 3: input everything
            %
            % Inputs: (optional)
            %   stepsize:       positive scalar uniform stepsize for initial construction
            %   lower_bound:    lower bound of interval
            %   upper_bound:    upper bound of interval
            %   dof:            degrees of freedom as property
            %
            % Output:
            %   Mesh1d object
            if nargin < 1
                lower_bound = 0;
                upper_bound = 1;
                stepsize = 0.1;
                obj.dof = 2; 
            elseif nargin == 1
                lower_bound = 0;
                upper_bound = 1;
                obj.dof = 2;                 
            end
            obj.number_of_nodes = length(obj.nodes);
            obj.boundary_nodes_idx = [1; obj.number_of_nodes];
            obj.connectivity_matrix = [1:obj.dof:(obj.number_of_nodes-1), 2:obj.dof:obj.number_of_nodes];
        end
        function obj = Mesh1d(lower_bound, upper_bound, h)
            %UNTITLED15 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end