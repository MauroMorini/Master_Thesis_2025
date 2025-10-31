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
        wave_speed_type = "brute-force"                 % helps deciding matrix update scheme
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
            final_time = 5;
            has_exact_solution = true;

            % symbolic calculations
            syms x t
            [wave_speed_coeff_sym, wave_speed_type] = fem1d.PDEData.generateSmoothWaveSpeed(wave_speed_index);
            u_exact_sym = exp(-(x-t+1)^2);
            % u_exact_sym = cos(2*pi*(x - t));
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
            u_exact_fun = @(x,t) u_exact_fun(x,t) + zeros(size(x));
            grad_u_exact_fun = @(x,t) grad_u_exact_fun(x,t) + zeros(size(x));
            rhs_fun = @(x,t) rhs_fun(x,t) + zeros(size(x));
            wave_speed_coeff_fun = @(x,t) wave_speed_coeff_fun(x,t) + zeros(size(x));
            u_t_exact_fun = @(x,t) u_t_exact_fun(x,t) + zeros(size(x));

            initial_displacement = @(x) u_exact_fun(x, initial_time);
            initial_velocity = @(x) u_t_exact_fun(x, initial_time);

            boundary_conditions = cell(1,2);
            boundary_conditions{1} = fem1d.BoundaryCondition1D("neumann", boundary_points(1), @(x,t) -grad_u_exact_fun(x,t) );
            boundary_conditions{2} = fem1d.BoundaryCondition1D("transparent", boundary_points(2), u_exact_fun);
            pde_data = fem1d.PDEData(initial_time, final_time, initial_displacement, initial_velocity, boundary_points,...
                                boundary_conditions, rhs_fun, wave_speed_coeff_fun);
            pde_data.u_exact_fun = u_exact_fun;
            pde_data.grad_u_exact_fun = grad_u_exact_fun;
            pde_data.has_exact_solution = has_exact_solution;
            pde_data.wave_speed_type = wave_speed_type;
        end

        function [pde_data, resonator_matrix] = generate_gaussian_puls_data_on_waveguide_with_resonators()

            boundary_points = [0, 10];
            initial_time = 0;
            final_time = 20;
            has_exact_solution = false;
            resonator_matrix = [6, 7; 9, 9.5];
            wave_speed_type = "piecewise-const-coefficient-in-space";

            % symbolic calculations
            syms x t
            u_exact_sym = exp(-(x-t+2)^2);
            % u_exact_sym = cos(2*pi*(x - t));
            u_t_exact_sym = diff(u_exact_sym, t);
            grad_u_exact_sym = diff(u_exact_sym, x);
            rhs_sym = diff(u_t_exact_sym,t) - diff(grad_u_exact_sym, x);

            % transform into matlab function
            u_exact_fun = matlabFunction(u_exact_sym, 'Vars', {x,t});
            grad_u_exact_fun = matlabFunction(grad_u_exact_sym, 'Vars', {x,t});
            u_t_exact_fun = matlabFunction(u_t_exact_sym, 'Vars', {x,t});
            rhs_fun = matlabFunction(rhs_sym, 'Vars', {x,t});
            wave_speed_coeff_fun = @(x,t) 1*(x < resonator_matrix(1,1) | (x > resonator_matrix(1,2) & x < resonator_matrix(2,1)) | x > resonator_matrix(2,2)) +... 
                (1+0.2*cos(pi/2*t))/0.1.*( (x >= resonator_matrix(1,1) & x <= resonator_matrix(1,2)) | (x >= resonator_matrix(2,1) & x <= resonator_matrix(2,2)) );
            wave_speed_coeff_fun = @(x,t) ones(size(x));
            % check for scalar functions
            u_exact_fun = @(x,t) u_exact_fun(x,t) + zeros(size(x));
            grad_u_exact_fun = @(x,t) grad_u_exact_fun(x,t) + zeros(size(x));
            rhs_fun = @(x,t) rhs_fun(x,t) + zeros(size(x));
            wave_speed_coeff_fun = @(x,t) wave_speed_coeff_fun(x,t) + zeros(size(x));
            u_t_exact_fun = @(x,t) u_t_exact_fun(x,t) + zeros(size(x));

            initial_displacement = @(x) u_exact_fun(x, initial_time);
            initial_velocity = @(x) u_t_exact_fun(x, initial_time);

            initial_displacement = @(x) zeros(size(x));
            initial_velocity = @(x) zeros(size(x));
            rhs_fun = @(x, t) zeros(size(x));

            boundary_conditions = cell(1,2);
            boundary_conditions{1} = fem1d.BoundaryCondition1D("neumann", boundary_points(1), @(x,t) -grad_u_exact_fun(x,t) );
            boundary_conditions{2} = fem1d.BoundaryCondition1D("transparent", boundary_points(2), u_exact_fun);
            pde_data = fem1d.PDEData(initial_time, final_time, initial_displacement, initial_velocity, boundary_points,...
                                boundary_conditions, rhs_fun, wave_speed_coeff_fun);
            pde_data.u_exact_fun = u_exact_fun;
            pde_data.grad_u_exact_fun = grad_u_exact_fun;
            pde_data.has_exact_solution = has_exact_solution;
            pde_data.wave_speed_type = wave_speed_type;
        end

        function [c_sym, type] = generateSmoothWaveSpeed(index)
            syms x t
            c_cell = {  struct("c_sym", x*0+1, "type", "time-independent");
                        struct("c_sym", (sin(x) + 2)*(cos(t)+2), "type", "brute-force");
                        struct("c_sym", (sin(x) + 2), "type", "time-independent");
                        struct("c_sym", (sin(1*t) + 2), "type", "piecewise-const-coefficient-in-space");
            };
            c_sym = c_cell{index}.c_sym;
            type = c_cell{index}.type;
        end

        function [c_sym, type] = generateDiscontinuousWaveSpeed(resIdx)
            syms x t
            c_cell = {  struct("c_sym", x*0+1, "type", "time-independent");
                        struct("c_sym", (sin(x) + 2)*(cos(t)+2), "type", "brute-force");
                        struct("c_sym", (sin(x) + 2), "type", "time-independent");
                        struct("c_sym", (sin(1*t) + 2), "type", "piecewise-const-coefficient-in-space");
            };
            c_sym = c_cell{index}.c_sym;
            type = c_cell{index}.type;
        end
    end
end