classdef BoundaryCondition1D
    % describes a boundary condition at a certain boundary 
    % could be time dependent or not 
    % and Neumann, Dirichlet or even Absorbing
    properties
        bc_type                 % string, could be "neumann", "dirichlet", "transparent"
        bc_location             % scalar point where b.c. is
        bc_fun                  % @(x,t) function handle of rhs value 
        outward_normal          % scalar in {1, -1} 
    end

    methods
        function obj = BoundaryCondition1D(bc_type, bc_location, bc_fun)
            % constructor
            arguments (Input)
                bc_type string 
                bc_location double 
                bc_fun function_handle
            end
            obj.bc_type = bc_type;
            obj.bc_location = bc_location;
            obj.bc_fun = bc_fun;
        end

        function bc_val = get_bc_val(obj, time)
            bc_val = obj.bc_fun(obj.bc_location, time);
        end
    end
end