classdef SIPGWaveSolver1D < handle
    % solver class for the full problem of solving the wave
    % equation with time-dependent coefficients using 
    % SIPG in space and leapfrog in time for constant meshes.
    % (method of lines)
    %
    properties
        sigma double = 0                    % scalar penalty constant "sigma" for SIPG 
        initial_mesh                        % MeshIntervalDG1d object
        quadrature_mesh                     % MeshIntervalDG1d object for quadrature 
        pde_data                            % PDEData1D object
        solution                            % (num_nodes, num_steps) solution matrix
        initial_matrix_struct               % struct as a collection of assembled initial matrices
        dt double = 0                       % global stepsize for time integration    
        times (1,:) double                  % time vector containing all steps made                             
    end

    methods
        function obj = SIPGWaveSolver1D(initial_mesh, pde_data)
            % constructor 
            obj.initial_mesh = initial_mesh;
            obj.pde_data = pde_data;
        end

        function obj = initialize_quadrature(obj)
            obj.quadrature_mesh = copy(obj.initial_mesh);
            obj.quadrature_mesh.dof = obj.quadrature_mesh.dof + 1;
        end

        function obj = setup_settings(obj)
            if obj.sigma <= 0
                obj.sigma = 10*obj.initial_mesh.dof^2;
            end
        end

        function obj = assemble_initial_matrices(obj)
            % collects system matrices into a struct for later time adaptation

            % initialization
            t0 = obj.pde_data.initial_time;
            c_fun_initial = @(x) obj.pde_data.wave_speed_coeff_fun(x, t0);

            % initialize c_vals
            [nodes_quad, ~, elements_quad] = obj.quadrature_mesh.getPet();
            c_vals = c_fun_initial(nodes_quad(elements_quad));

            % get matrices
            [nodes, ~, elements] = obj.initial_mesh.getPet();
            A = fem1d.stiffnessMatrix1D(nodes, elements, c_vals);
            M = fem1d.massMatrix1D(nodes, elements, ones(size(c_vals)));
            B_flux_int = dg1d.interiorFluxMatrix1D(nodes, elements, c_vals);
            B_flux_bound = dg1d.boundaryFluxMatrix1D(nodes, elements, c_vals);
            B_penalty_int = dg1d.interiorPenaltyMatrix1D(nodes, elements, c_vals, obj.sigma);
            B_penalty_bound = dg1d.boundaryPenaltyMatrix1D(nodes, elements, c_vals, obj.sigma);

            obj.initial_matrix_struct = struct("A", A, "M", M, "B_flux_int", B_flux_int, "B_flux_bound", B_flux_bound, ...
                                         "B_penalty_int", B_penalty_int, "B_penalty_bound", B_penalty_bound);
            
        end

        function obj = calculate_stable_stepsize(obj)
            % calculates a stable stepsize dt if one has not already been chosen using the initial mass and stiffness matrices
            if obj.dt > 0
                return
            end
            B_loc = obj.initial_matrix_struct.A - obj.initial_matrix_struct.B_flux_bound - obj.initial_matrix_struct.B_flux_int + obj.initial_matrix_struct.B_penalty_int + obj.initial_matrix_struct.B_penalty_bound;
            M_loc = obj.initial_matrix_struct.M;
            lMax = eigs(B_loc,1);
            lMin = eigs(M_loc,1,0);
            obj.dt = sqrt(lMin/lMax);
        end

        function system_struct = setup_system(obj, current_time)
            % sets up system by introducing time dependent components and boundary conditions
            % TODO: update stiffness matrix in time
            system_struct = obj.initial_matrix_struct;
            system_struct.time = current_time;
            
            % collect load vector
            f_values = obj.pde_data.rhs_fun(obj.quadrature_mesh.nodes(obj.quadrature_mesh.elements), current_time);
            system_struct.load_vector = fem1d.loadVector1D(obj.initial_mesh.nodes, obj.initial_mesh.elements, f_values);
            system_struct = obj.impose_boundary_conditions(system_struct);
        end

        function x = solve_system(obj, A, b)
            % solves system Ax = b
            x = A\b;
        end

        function system_struct = impose_boundary_conditions(obj, system_struct)
            % imposes boundary condition at a fixed time to a already prepared system
            arguments (Input)
                obj
                system_struct struct            % contains time dependent matrices and vectors for the system
            end
            current_time = system_struct.time;

            % initializations
            [nodes_quad, ~, elements_quad] = obj.quadrature_mesh.getPet();
            [nodes, ~, elements] = obj.initial_mesh.getPet();
            [lower_boundary_element_idx, upper_boundary_element_idx] = obj.initial_mesh.getBoundaryElementIdx();
            boundary_element_idx = [lower_boundary_element_idx, upper_boundary_element_idx];
            c_vals = obj.pde_data.wave_speed_coeff_fun(nodes_quad(elements_quad), current_time);
            g_vals = [obj.pde_data.boundary_conditions{1}.get_bc_val(current_time), ...
                        obj.pde_data.boundary_conditions{2}.get_bc_val(current_time)];
            
            % collect b.c vectors
            dirichlet_vector = dg1d.dirichletbcVector1D(nodes, elements, c_vals, g_vals,obj. sigma);
            neumann_vector = dg1d.neumannbcVector1D(nodes, elements, c_vals, g_vals);

            for i = 1:2
                bc_type = obj.pde_data.boundary_conditions{i}.bc_type;
                switch bc_type
                    case "dirichlet"
                        neumann_vector(elements(boundary_element_idx(i),:)) = 0;
                    case "neumann"
                        dirichlet_vector(elements(boundary_element_idx(i),:)) = 0;
                        system_struct.B_flux_bound(elements(boundary_element_idx(i),:), elements(boundary_element_idx(i),:)) = 0;
                        system_struct.B_penalty_bound(elements(boundary_element_idx(i),:), elements(boundary_element_idx(i),:)) = 0;
                    otherwise
                        error("the boundary condition type << " + bc_type + " >> has not been implemented")
                end
            end
            system_struct.load_vector = system_struct.load_vector + dirichlet_vector + neumann_vector;
            system_struct.B = system_struct.A - system_struct.B_flux_bound - system_struct.B_flux_int + system_struct.B_penalty_int + system_struct.B_penalty_bound;
        end

        function obj = implement_initial_conditions(obj)
            % introduces the already known initial displacement into solution and calculates the solution at the time t0 + dt
            % using a taylor-expansion
            u0 = obj.pde_data.initial_displacement(obj.initial_mesh.nodes);
            v0 = obj.pde_data.initial_velocity(obj.initial_mesh.nodes); 
            system_struct = obj.setup_system(obj.times(1));
            obj.solution(:, 1) = u0;
            system_matrix = system_struct.M;
            system_rhs = system_struct.M*u0 + obj.dt*system_struct.M*v0 + obj.dt^2/2*(system_struct.load_vector - system_struct.B*u0);
            obj.solution(:, 2) = obj.solve_system(system_matrix, system_rhs);
        end

        function obj = setup_iteration(obj)
            % main preparation function for the iteration, initializes and preallocates everything needed for the actual 
            % iteration
            obj.initialize_quadrature();
            obj.setup_settings();
            obj.assemble_initial_matrices();
            obj.calculate_stable_stepsize();

            obj.times = obj.pde_data.initial_time:obj.dt:obj.pde_data.final_time;
            obj.solution = zeros(size(obj.initial_mesh.nodes,1), size(obj.times, 2));
            obj.implement_initial_conditions();
        end

        function obj = leap_frog_leap(obj)
            % main iteration, applies leapfrog time integration
            
            for i = 2:length(obj.times)
                system_struct = obj.setup_system(obj.times(i));
                system_matrix = system_struct.M;
                system_rhs = obj.dt^2*system_struct.load_vector + ...
                            (2*system_struct.M - obj.dt^2*system_struct.B)*obj.solution(:, i) - ...
                             system_struct.M*obj.solution(:, i-1);
                obj.solution(:, i+1) = obj.solve_system(system_matrix, system_rhs);
            end
        end

        function obj = run(obj)
            obj.setup_iteration();
            obj.leap_frog_leap();
        end
    end
end