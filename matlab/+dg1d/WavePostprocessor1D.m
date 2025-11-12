classdef WavePostprocessor1D < handle
    % postprocessor for a 1d dg-leapfrog wave solution with a static (time-independent) mesh 
    properties
        solution
        mesh
        time_vector
        sigma
        pde_data
        metadata
        errors_obj
    end

    methods
        function obj = WavePostprocessor1D(wave_solver)
            arguments (Input)
                wave_solver = [];
            end
            if isempty(wave_solver)
                return;
            end
            obj.solution = wave_solver.solution;
            obj.mesh = wave_solver.initial_mesh;
            obj.time_vector = wave_solver.time_vector;
            obj.sigma = wave_solver.sigma;
            obj.pde_data = wave_solver.pde_data;
        end

        function obj = calculate_errors(obj)
            if ~obj.pde_data.has_exact_solution
                warning("no error could be calculated since there is no exact solution")
                return
            end
            [uh, T] = obj.get_solution_at_time(obj.time_vector(end));
            obj.errors_obj = fem1d.Errors1D(@(x) obj.pde_data.u_exact_fun(x,T), @(x) obj.pde_data.grad_u_exact_fun(x,T), uh, obj.mesh);
            obj.errors_obj.initialize_dg_settings(@(x) obj.pde_data.wave_speed_coeff_fun(x,T), obj.sigma);
            obj.errors_obj.run();
        end

        function [u_h, t_h] = get_solution_at_time(obj, t)
            % for a given time t extracts the time t_h in the time_vector which is closest to t and 
            % returns the solution u_h(:, t_h)
            [~, t_h_index] = min(abs(obj.time_vector - t));
            u_h = obj.solution(:, t_h_index);
            t_h = obj.time_vector(t_h_index);
        end

        function obj = write_to_hdf5(obj, filename)
            arguments (Input)
                obj
                filename string
            end
            if isfile(filename)
                warning("a file with the name: " + filename + " already exists..... will be saved with the ending -temp")
                filename = erase(filename, ".h5");
                filename = filename + "-temp.h5";
                if isfile(filename)
                    delete(filename);
                    fprintf('overwriting %s ... \n', filename)
                end
            end

            % write mesh
            dof = obj.mesh.dof;
            h5create(filename, '/Mesh/dof', 1);
            h5write(filename, '/Mesh/dof', dof);
            element_interface_nodes = obj.mesh.element_interface_nodes;
            h5create(filename, '/Mesh/element_interface_nodes', size(element_interface_nodes));
            h5write(filename, '/Mesh/element_interface_nodes', element_interface_nodes);
            resonators_matrix = obj.mesh.resonators_matrix;
            h5create(filename, '/Mesh/resonators_matrix', size(resonators_matrix));
            h5write(filename, '/Mesh/resonators_matrix', resonators_matrix);
            H = [obj.mesh.h_max, obj.mesh.h_min];
            h5create(filename, '/Mesh/meshsizes', size(H));
            h5write(filename, '/Mesh/meshsizes', H);
            
            % write time 
            h5create(filename, '/time_vector', size(obj.time_vector));
            h5write(filename, '/time_vector', obj.time_vector);

            % write solution
            h5create(filename, '/solution', size(obj.solution));
            h5write(filename, '/solution', obj.solution);

            % write metadata
            h5create(filename, '/Metadata/sigma', size(obj.sigma));
            h5write(filename, '/Metadata/sigma', obj.sigma);
            wave_speed_string = string(func2str(obj.pde_data.wave_speed_coeff_fun));
            exact_solution_string = string(func2str(obj.pde_data.u_exact_fun));
            h5create(filename, '/Metadata/wave_speed', [1 1], 'Datatype', 'string');
            h5write(filename, '/Metadata/wave_speed', wave_speed_string);
            h5create(filename, '/Metadata/exact_solution', [1 1], 'Datatype', 'string');
            h5write(filename, '/Metadata/exact_solution', exact_solution_string);
        end  

        function obj = read_from_hdf5(obj, filename)
            if ~isfile(filename)
                error("a file with the name: " + filename + " doesn't exist")
            end
            obj.solution = h5read(filename,'/solution');
            obj.sigma = h5read(filename, '/Metadata/sigma');
            obj.time_vector = h5read(filename, '/time_vector');

            % recreate mesh
            element_interface_nodes = h5read(filename, '/Mesh/element_interface_nodes');
            H = h5read(filename, '/Mesh/meshsizes');
            obj.mesh = mesh.MeshIntervalDG1d([element_interface_nodes(1), element_interface_nodes(end)], H);
            obj.mesh.resonators_matrix = h5read(filename, '/Mesh/resonators_matrix');
            obj.mesh.element_interface_nodes = element_interface_nodes;
            obj.mesh.dof = h5read(filename, '/Mesh/dof');
            if ~isempty(obj.mesh.resonators_matrix)
                obj.mesh.is_resonator_mesh = true;
            end 
            obj.mesh.updatePet();

            % metadata string
            wave_speed = h5read(filename, '/Metadata/wave_speed');
            exact_solution = h5read(filename, '/Metadata/exact_solution');
            obj.metadata = sprintf(['Time interval: [%g, %g]\n' ...
                                    'Domain: (%g,%g)\n' ...
                                    'sigma: %g \n'...
                                    'Exact Solution: %s \n'...
                                    'Wave Speed: %s'], ...
                                    obj.time_vector(1),obj.time_vector(end),...
                                    obj.mesh.lower_interval_bound, obj.mesh.upper_interval_bound,...
                                    obj.sigma, exact_solution, wave_speed);

        end

        function [obj, figure_cell] = plot_solutions(obj, plot_times, fast_plot)
            arguments (Input)
                obj
                plot_times
                fast_plot = false;
            end
            num_fig = length(plot_times);
            if num_fig > 50
                error("Plotting  " + num_fig + " figures is too much choose plot_times smaller")
            end
            figure_cell = cell(1, num_fig);
            for i = 1:num_fig
                [uh, t] = obj.get_solution_at_time(plot_times(i));
                f = figure('Visible', 'on');
                figure_cell{i} = f;
                obj.mesh.plotDGsol(uh, f, fast_plot);
                title("t = " + t)
                ylim([-2,2])
            end
        end

        function obj = create_and_save_wave_animation(obj, filename, number_of_frames, ylim_vector)
            % creates a video of the solution 
            arguments (Input)
                obj
                filename
                number_of_frames = 1000
                ylim_vector = [-2,2]
            end
            v = VideoWriter(filename);
            v.FrameRate = 60;
            open(v);
            
            f = figure("Visible","off");
            f.Position = get(0, 'ScreenSize');
            ax = axes('Parent', f, 'Position',[0 0 1 1]);
            ylim(ax, ylim_vector)
            xlim(ax, [obj.mesh.lower_interval_bound, obj.mesh.upper_interval_bound])
            axis(ax, 'manual');
            
            t0 = obj.time_vector(1);
            T = obj.time_vector(end);
            plot_times = linspace(t0, T, number_of_frames);
            
            progress_percent = 0:10:100;
            progress_frames = floor(progress_percent/100 * number_of_frames);
            disp("starting to draw the frames")
            for i = 1:number_of_frames
                if any(i == progress_frames)
                    pct = progress_percent(progress_frames == i);
                    disp("Progress:  " + pct + "% -------------------------")
                end
                cla(ax);
                [uh, ~] = obj.get_solution_at_time(plot_times(i));
                f = obj.mesh.plotDGsol(uh, f);
                ax = findall(f, 'type', 'axes');
        
                frame = getframe(f);
                frame_resized = imresize(frame.cdata, [1080, 1920]);
                writeVideo(v, frame);
            end 
            close(v);
            disp("Video saved as " + filename);
        end
    end
end