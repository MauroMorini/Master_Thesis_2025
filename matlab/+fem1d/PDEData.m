classdef PDEData < handle
    % class with function handles containing coefficents, boundary-/initial conditions
    % and rhs load function of the PDE as well as domain information.
    % for now only for the wave equation u_tt - (c(x,t) u_x(x,t))_x = f(x,t)
    % and only for a interval domain
    properties
        initial_time                                    % scalar initial time 
        final_time                                      % scalar final time
        initial_displacement                            % @(x) function handle of u(x,0)
        initial_velocity                                % @(x) u_t(x,0)
        boundary_points                                 % (1, 2) vector with lower boundary point in first and upper boundary point in second entry
        boundary_conditions                             % (1, 2) cell with lower b_c in first and upper bc in second entry
        rhs_fun                                         % @(x,t) function handle of load (rhs)
        wave_speed_coeff_fun                            % @(x,t) function handle of coefficient c in space-time
        u_exact_fun                                     % @(x,t) exact solution function handle
        grad_u_exact_fun                                % @(x,t) spacial derivative function handle of exact solution
        has_exact_solution logical                      % true if u_exact_fun and grad_u_exact_fun are initialized 
        wave_speed_is_time_dependent = NaN              % helps deciding matrix update scheme
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

    methods (Static)
        function pde_data = generate_gaussian_puls_data_on_waveguide(wave_speed_index)
            arguments (Input)
                wave_speed_index = 1;
            end
            boundary_points = [0, 10];
            initial_time = 0;
            final_time = 2;
            has_exact_solution = true;

            % symbolic calculations
            syms x t
            [wave_speed_coeff_sym, wave_speed_is_time_dependent] = fem1d.PDEData.generateWaveSpeed(wave_speed_index);
            u_exact_sym = exp(-(x-5*t+1)^2);
            u_t_exact_sym = diff(u_exact_sym, t);
            grad_u_exact_sym = diff(u_exact_sym, x);
            rhs_sym = diff(u_t_exact_sym,t) - diff(wave_speed_coeff_sym*grad_u_exact_sym, x);

            % transform into matlab function
            u_exact_fun = matlabFunction(u_exact_sym, 'Vars', {x,t});
            grad_u_exact_fun = matlabFunction(grad_u_exact_sym, 'Vars', {x,t});
            u_t_exact_fun = matlabFunction(u_t_exact_sym, 'Vars', {x,t});
            rhs_fun = matlabFunction(rhs_sym, 'Vars', {x,t});
            wave_speed_coeff_fun = matlabFunction(wave_speed_coeff_sym, 'Vars', {x,t});

            % check for scalar functions
            if isempty(symvar(u_exact_sym))
                u_exact_fun = @(x,t) u_exact_fun(x,t) + zeros(size(x));
            end
            if isempty(symvar(grad_u_exact_sym))
                grad_u_exact_fun = @(x,t) grad_u_exact_fun(x,t) + zeros(size(x));
            end
            if isempty(symvar(rhs_sym))
                rhs_fun = @(x,t) rhs_fun(x,t) + zeros(size(x));
            end
            if isempty(symvar(wave_speed_coeff_sym))
                wave_speed_coeff_fun = @(x,t) wave_speed_coeff_fun(x,t) + zeros(size(x));
            end
            if isempty(symvar(u_t_exact_sym))
                u_t_exact_fun = @(x,t) u_t_exact_fun(x,t) + zeros(size(x));
            end

            initial_displacement = @(x) u_exact_fun(x, initial_time);
            initial_velocity = @(x) u_t_exact_fun(x, initial_time);

            boundary_conditions = cell(1,2);
            boundary_conditions{1} = fem1d.BoundaryCondition1D("dirichlet", boundary_points(1), u_exact_fun);
            boundary_conditions{2} = fem1d.BoundaryCondition1D("dirichlet", boundary_points(2), u_exact_fun);
            pde_data = fem1d.PDEData(initial_time, final_time, initial_displacement, initial_velocity, boundary_points,...
                                boundary_conditions, rhs_fun, wave_speed_coeff_fun);
            pde_data.u_exact_fun = u_exact_fun;
            pde_data.grad_u_exact_fun = grad_u_exact_fun;
            pde_data.has_exact_solution = has_exact_solution;
            pde_data.wave_speed_is_time_dependent = wave_speed_is_time_dependent;
        end

        function [c_sym, isTimeDependent] = generateWaveSpeed(index)
            syms x t
            c_cell = {  struct("c_sym", x*0+1, "isTimeDependent", false);
                        struct("c_sym", (sin(x) + 2)*(cos(t)+2), "isTimeDependent", true);
                        struct("c_sym", (sin(x) + 2), "isTimeDependent", false);
            };
            c_sym = c_cell{index}.c_sym;
            isTimeDependent = c_cell{index}.isTimeDependent;
        end
    end
end