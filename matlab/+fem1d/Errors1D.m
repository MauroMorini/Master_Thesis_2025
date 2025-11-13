classdef Errors1D < handle
    % calculates L2 and H1 error in space for a 1d FEM solution 
    properties
        u_exact_fun
        grad_u_exact_fun
        c_fun
        uh
        mesh
        quadrature_mesh
        uh_quad
        l2_error
        h1_error
        energy_error = NaN
        sigma = NaN
    end

    methods
        function obj = Errors1D(u_exact_fun, grad_u_exact_fun, uh, mesh)
            if nargin == 0
                return
            end
            obj.u_exact_fun = u_exact_fun;
            obj.grad_u_exact_fun = grad_u_exact_fun;
            obj.uh = uh;
            obj.mesh = mesh;
        end

        function [l2_error, h1_error, energy_error] = getErrors(obj)
            l2_error = obj.l2_error;
            h1_error = obj.h1_error;
            energy_error = obj.energy_error;
        end

        function obj = initialize_dg_settings(obj, c_fun, sigma)
            obj.sigma = sigma;
            obj.c_fun = c_fun;
        end

        function obj = generate_quadrature_mesh(obj, additional_quadrature_dof)
            % creates quadrature mesh with more nodes than original mesh
            arguments (Input)
                obj
                additional_quadrature_dof = 1;
            end
            obj.quadrature_mesh = copy(obj.mesh);
            obj.quadrature_mesh.dof = obj.quadrature_mesh.dof + additional_quadrature_dof;
            obj.quadrature_mesh.updatePet();
        end

        function obj = interpolate_sol_at_quad_mesh(obj)
            % interpolates solution at higher dof quadrature mesh
            bary_weights = common.calculateBarycentricWeights(obj.mesh.nodes, obj.mesh.elements);
            uh_quad_temp = zeros(size(obj.quadrature_mesh.nodes));
            el_quad_temp = obj.quadrature_mesh.elements;
            num_el = size(el_quad_temp, 1);
            for k = 1:num_el
                quadrature_nodes_loc = obj.quadrature_mesh.nodes(el_quad_temp(k,:));
                basis_nodes_loc = obj.mesh.nodes(obj.mesh.elements(k,:));
                Phi = common.evaluateLagrangeBarycentric(quadrature_nodes_loc, bary_weights(obj.mesh.elements(k,:)), basis_nodes_loc);
                uh_quad_temp(el_quad_temp(k,:)) = Phi.'*obj.uh(obj.mesh.elements(k, :));
            end
            obj.uh_quad = uh_quad_temp;
        end

        function obj = calculate_errors(obj)
            u_exact_vals = obj.u_exact_fun(obj.quadrature_mesh.nodes);
            grad_u_exact_vals = obj.grad_u_exact_fun(obj.quadrature_mesh.nodes);
            [obj.l2_error, obj.h1_error] = fem1d.errors1DWithExactSol(obj.quadrature_mesh.nodes, obj.quadrature_mesh.elements, obj.uh_quad, u_exact_vals, grad_u_exact_vals);
            
            if ~isnan(obj.sigma)
                c_vals = obj.c_fun(obj.quadrature_mesh.nodes(obj.quadrature_mesh.elements));
                obj.energy_error = dg1d.energyNormError1D(obj.quadrature_mesh.nodes, obj.quadrature_mesh.elements, obj.uh_quad, c_vals, obj.sigma, u_exact_vals, grad_u_exact_vals);
            end
        end

        function obj = write_errors_to_csv(obj, filename, errors, meshsizes, metadata_input)
            % takes a matrix of multiple errors and meshsizes and creates a csv file 
            % writing into it the errors and meshsizes and commented lines for metadata
            arguments (Input)
                obj
                filename string
                errors                      % (num_errors, num_meshsizes)
                meshsizes                   % (1, num_meshsizes)
                metadata_input string       % string containing more metadata given from the outside
            end
            if isfile(filename)
                warning("a file with the name: " + filename + " already exists..... will be saved with the ending -temp")
                filename = erase(filename, ".csv");
                filename = filename + "-temp.csv";
                if isfile(filename)
                    delete(filename);
                    fprintf('overwriting %s ... \n', filename)
                end
            end

            assert(size(errors, 2) == size(meshsizes, 2), "there have to be as many meshsizes as error columns");

            % collect metadata
            dof = obj.mesh.dof;
            quad_dof = obj.quadrature_mesh.dof;
            cr = obj.interpolate_convergence_rates(meshsizes, errors);

            file_id = fopen(filename, 'w');
            
            % write metadata 
            fprintf(file_id, '# Errors of SIPG numerical solution \n');
            fprintf(file_id, '# Domain: (%g, %g) \n', obj.mesh.lower_interval_bound, obj.mesh.upper_interval_bound);
            fprintf(file_id, '# DoF of the solution (per cell): %g \n', dof);
            fprintf(file_id, '# DoF of the quadrature mesh: %g \n', quad_dof);
            fprintf(file_id, '# sigma: %g \n', obj.sigma);
            fprintf(file_id, '# Resonators: %s \n', mat2str(obj.mesh.resonators_matrix));
            fprintf(file_id, '\n%s \n', metadata_input);
            fprintf(file_id, '# Interpolated Convergence Rates ------------- \n# l2: %g,  h1: %g,  energy: %g \n\n', cr(2,1),cr(2,2),cr(2,3));

            % write errors
            fprintf(file_id, 'meshsize,L2-error,H1-error,energy-error \n');
            fclose(file_id);
            M = [meshsizes', errors'];
            writematrix(M, filename, 'Delimiter', ',', 'WriteMode', 'append');

        end

        function obj = run(obj)
            obj.generate_quadrature_mesh();
            obj.interpolate_sol_at_quad_mesh();
            obj.calculate_errors();
        end
    end

    methods (Static)
        function convergence_rates = interpolate_convergence_rates(meshsizes, errors)
            % calculate convergence rates using linear interpolation (least squares)
            arguments (Input)
                meshsizes (1,:)             % (1, num_meshsizes)
                errors                      % (num_error_types, num_meshsizes)
            end
            assert(size(errors, 2) == size(meshsizes, 2), "there must be as many error measurements as meshsizes")
            X = [ones(length(meshsizes),1), log(meshsizes')];
            convergence_rates = X\log(errors');
        end

        function [meshsizes, errors, metadata] = read_errors_from_csv(filename)
            % reads a csv file created by write_errors function 
            arguments (Input)
                filename string
            end
            arguments (Output)
                meshsizes (1,:)     % (1, num_meshsizes) matrix with corresponding meshsizes
                errors              % (num_error_types, num_meshsizes) error matrix
                metadata string     % a string containing all the metadata for output (not processed)
            end

            % get data
            data = readmatrix(filename, "CommentStyle", "#");
            errors = data(:, 2:end)';
            meshsizes = data(:, 1)';

            % get metadata
            f_id = fopen(filename, 'r');
            comments = textscan(f_id, '#%[^\n]');
            fclose(f_id);
            metadata = strjoin(strcat('#', comments{1}), newline);
        end
    end
end