classdef PDEData
    % class with function handles containing coefficents, boundary-/initial conditions
    % and rhs load function of the PDE as well as domain information.
    % for now only for the wave equation u_tt - (c(x,t) u_x(x,t))_x = f(x,t)
    % and only for a interval domain
    properties
        initial_time                % scalar initial time 
        final_time                  % scalar final time
        initial_displacement        % @(x) function handle of u(x,0)
        initial_velocity            % @(x) u_t(x,0)
        boundary_points             % (1, 2) vector with lower boundary point in first and upper boundary point in second entry
        boundary_conditions         % (1, 2) cell with lower b_c in first and upper bc in second entry
        rhs_fun                     % @(x,t) function handle of load (rhs)
        wave_speed_coeff_fun        % @(x,t) function handle of coefficient c in space-time
    end 

    methods
        function obj = PDEData(initial_time, final_time, initial_displacement, initial_velocity, boundary_points,...
                                boundary_conditions, rhs_fun, wave_speed_coeff_fun)
            % constructor 
            arguments (Input)
                initial_time double
                final_time double
                initial_displacement function_handle
                initial_velocity function_handle
                boundary_points (1,2) double
                boundary_conditions cell
                rhs_fun function_handle
                wave_speed_coeff_fun function_handle
            end
            obj.initial_time = initial_time;
            obj.final_time = final_time;
            obj.initial_displacement= initial_displacement;
            obj.initial_velocity = initial_velocity;
            obj.boundary_points = boundary_points;
            obj.boundary_conditions = boundary_conditions;
            obj.rhs_fun = rhs_fun;
            obj.wave_speed_coeff_fun = wave_speed_coeff_fun;
        end
    end
end