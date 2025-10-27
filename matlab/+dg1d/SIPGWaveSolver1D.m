classdef SIPGWaveSolver1D < handle
    % solver class for the full problem of solving the wave
    % equation with time-dependent coefficients using 
    % SIPG in space and leapfrog in time for constant meshes.
    % (method of lines)
    %
    properties
        sigma               % scalar penalty constant "sigma" for SIPG 
        initial_mesh        % MeshIntervalDG1d object
        quadrature_mesh     % MeshIntervalDG1d object for quadrature 
        pde_data            % PDEData1D object
        solution            % (num_nodes, num_steps) solution matrix
        matrix_struct       % struct as a collection of assembled initial matrices
    end

    methods
        function obj = SIPGWaveSolver(initial_mesh, pde_data)
            % constructor 
            obj.initial_mesh = initial_mesh;
            obj.pde_data = pde_data;
        end

        function obj = initialize_quadrature(obj)
            obj.quadrature_mesh = copy(obj.initial_mesh);
            obj.quadrature_mesh.dof = obj.quadrature_mesh.dof + 1;
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

            obj.matrix_struct = struct("A", A, "M", M, "B_flux_int", B_flux_int, "B_flux_bound", B_flux_bound, ...
                                         "B_penalty_int", B_penalty_int, "B_penalty_bound", B_penalty_bound);
        end

        function obj = setup_current_system(input)
            
        end

        function system_struct = impose_boundary_conditions(obj, system_struct)
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
            system_struct.system_rhs = system_struct.system_rhs + dirichlet_vector + neumann_vector;
        end
    end
end