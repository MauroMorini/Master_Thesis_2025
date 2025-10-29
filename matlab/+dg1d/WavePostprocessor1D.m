classdef WavePostprocessor1D < handle
    % postprocessor for a 1d dg-leapfrog wave solution with a static (time-independent) mesh 
    properties
        wave_solver
        errors_obj
    end

    methods
        function obj = WavePostprocessor1D(wave_solver)
            obj.wave_solver = wave_solver;
        end

        function obj = calculate_errors(obj)
            if ~obj.wave_solver.pde_data.has_exact_solution
                obj.l2_error = NaN;
                obj.h1_error = NaN;
                obj.energy_error = NaN;
                warning("no error could be calculated since there is no exact solution")
                return
            end
            [uh, T] = obj.get_solution_at_time(obj.wave_solver.time_vector(end));
            obj.errors_obj = fem1d.Errors1D(@(x) obj.wave_solver.pde_data.u_exact_fun(x,T), @(x) obj.wave_solver.pde_data.grad_u_exact_fun(x,T), uh, obj.wave_solver.initial_mesh);
            obj.errors_obj.initialize_dg_settings(@(x) obj.wave_solver.pde_data.wave_speed_coeff_fun(x,T), obj.wave_solver.sigma);
            obj.errors_obj.run();
        end

        function [u_h, t_h] = get_solution_at_time(obj, t)
            % for a given time t extracts the time t_h in the time_vector which is closest to t and 
            % returns the solution u_h(:, t_h)
            [~, t_h_index] = min(abs(obj.wave_solver.time_vector - t));
            u_h = obj.wave_solver.solution(:, t_h_index);
            t_h = obj.wave_solver.time_vector(t_h_index);
        end
    end
end