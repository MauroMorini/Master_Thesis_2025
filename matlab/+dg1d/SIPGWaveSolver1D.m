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
        neutered_matrix_struct              % struct containing initial matrices with removed values at resonator boundaries and divided by c locally for piecewise-const process 
        dt double = 0                       % global stepsize for time integration    
        time_vector (1,:) double            % time vector containing all steps made 
        matrix_update_type = []             % dictates how matrices are updated in each time step
                                            % in {"time-independent", "brute-force", "piecewise-const-coefficient-in-space"} 
        matrix_resonator_masks              % struct, see recompute_matrix_resonator_masks() for definition       
        wave_speed_cell                     % (1,2) cell with the c_vals calculated at the initial time int the first entry and for the current time in the second entry                  
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
            obj.quadrature_mesh.updatePet();
        end

        function obj = setup_settings(obj)
            if obj.sigma <= 0
                obj.sigma = 10*obj.initial_mesh.dof^2;
            end
            if isempty(obj.matrix_update_type)
                obj.matrix_update_type = obj.pde_data.wave_speed_type;
            end
        end

        function obj = precompute_matrix_resonator_masks(obj, system_struct)
            % precomputes the indices of the resonators for the respective triplets of A, B_flux, B_penalty
            % to reuse in each piecewise_const_coeff matrix update

            % initialization
            num_res = size(obj.initial_mesh.resonators_matrix, 1);
            A_masks = cell(1, num_res+1);
            B_flux_masks = cell(1, num_res+1);
            B_penalty_masks = cell(1, num_res+1);
            
            % dissect sparse matrices
            [A_rows, A_cols, ~] = find(system_struct.A);
            [B_flux_rows, B_flux_cols, ~] = find(system_struct.B_flux_int);
            [B_penalty_rows, B_penalty_cols, ~] = find(system_struct.B_penalty_int);

            % iterate over resonators
            for i = 0:num_res

                % find global dof indices which are in the current resonator
                res_idx = find(obj.initial_mesh.element_idx_to_resonator_idx_map == i);
                if isempty(res_idx)
                    continue
                end
                res_idx = obj.initial_mesh.elements(res_idx,:);
                res_idx = res_idx(:);      
                
                % find triplet indices, such that both rows and cols are in the resonator
                A_masks{i+1} = find(ismember(A_rows,res_idx) & ismember(A_cols,res_idx));
                B_flux_masks{i+1} = find(ismember(B_flux_rows,res_idx) & ismember(B_flux_cols,res_idx));
                B_penalty_masks{i+1} = find(ismember(B_penalty_rows,res_idx) & ismember(B_penalty_cols,res_idx));           
            end
            obj.matrix_resonator_masks = struct("A_masks", A_masks, "B_flux_masks", B_flux_masks, "B_penalty_masks", B_penalty_masks);
        end

        function obj = neuter_initial_matrix_struct(obj)
            % is called to initialize neutered_matrix_struct by subtracting out the boundary face contribution of the resonator boundary
            % and then dividing by c such that we get matrices prepared for piecewise-const wave speed update

            system_struct = obj.initial_matrix_struct;

            % dissect sparse matrices
            [A_rows, A_cols, A_vals] = find(system_struct.A);
            [B_flux_rows, B_flux_cols, B_flux_vals] = find(system_struct.B_flux_int);
            [B_penalty_rows, B_penalty_cols, B_penalty_vals] = find(system_struct.B_penalty_int);

            % recalculate interior resonator interface values for flux and penalty matrix (updating connecting information between resonators)
            res_boundary_upper_element_idx = obj.initial_mesh.getInteriorResonatorBoundaryFace();

            % extend triplets
            dof = obj.initial_mesh.dof;
            empty_triplet_slots = zeros(dof^2*4*length(res_boundary_upper_element_idx), 1);
            flux_triplet_iterator = length(B_flux_rows) + 1;
            penalty_triplet_iterator = length(B_penalty_rows) + 1;
            B_flux_rows = [B_flux_rows; empty_triplet_slots];
            B_flux_cols = [B_flux_cols; empty_triplet_slots];
            B_flux_vals = [B_flux_vals; empty_triplet_slots];
            B_penalty_rows = [B_penalty_rows; empty_triplet_slots];
            B_penalty_cols = [B_penalty_cols; empty_triplet_slots];
            B_penalty_vals = [B_penalty_vals; empty_triplet_slots];

            % build local matrices for B_flux (still need to multiply c, (1/h) 
            [phi_val, dphi_val, ~] = common.getShapeFunctionValueMatrix(dof);
            B_flux_ref_1 = 1*(phi_val(:,end)*dphi_val(:,end).' + dphi_val(:, end)*phi_val(:,end).');
            B_flux_ref_21 = 1*phi_val(:,end)*dphi_val(:,1).';
            B_flux_ref_22 = (-1)*dphi_val(:,end)*phi_val(:,1).';
            B_flux_ref_3 = (-1)*(phi_val(:,1)*dphi_val(:,1).' + dphi_val(:, 1)*phi_val(:,1).');

            B_penalty_ref_1 = (1)*phi_val(:, end)*phi_val(:,end).';
            B_penalty_ref_2 = (-1)*phi_val(:, end)*phi_val(:,1).';
            B_penalty_ref_3 = (1)*phi_val(:,1)*phi_val(:,1).';
            for k = res_boundary_upper_element_idx.'
            
                % initializations for both 
                elements_loc = obj.initial_mesh.elements(k-1:k,:);
                nodes_loc = obj.initial_mesh.nodes(elements_loc);
                h_loc = [abs(nodes_loc(1,1) - nodes_loc(1,end)); abs(nodes_loc(2,1) - nodes_loc(2,end))];
                c_vals_loc = [obj.wave_speed_cell{2}(k-1,end), obj.wave_speed_cell{1}(k-1,end);
                              obj.wave_speed_cell{2}(k,1), obj.wave_speed_cell{1}(k,1)];
                bordering_elements = [elements_loc(1,:), elements_loc(2,:)];
                B_loc_row_idx = repmat(bordering_elements.', 1, 2*dof);
                B_loc_col_idx = repmat(bordering_elements, 2*dof, 1);
                
                % flux update
                B_flux_loc_minus = [c_vals_loc(1,2)*B_flux_ref_1/h_loc(1), (c_vals_loc(2,2)*B_flux_ref_21/h_loc(2) + c_vals_loc(1,2)*B_flux_ref_22/h_loc(1));
                                  (c_vals_loc(2, 2)*B_flux_ref_21/h_loc(2) + c_vals_loc(1, 2)*B_flux_ref_22/h_loc(1)).', c_vals_loc(2, 2)*B_flux_ref_3/h_loc(2)];
                B_flux_loc = -B_flux_loc_minus;
                B_flux_vals(flux_triplet_iterator:(flux_triplet_iterator + 4*dof^2 - 1)) = B_flux_loc(:);
                B_flux_rows(flux_triplet_iterator:(flux_triplet_iterator + 4*dof^2 - 1)) = B_loc_row_idx(:);
                B_flux_cols(flux_triplet_iterator:(flux_triplet_iterator + 4*dof^2 - 1)) = B_loc_col_idx(:);
                flux_triplet_iterator = flux_triplet_iterator + 4*dof^2;

                % penalty update
                h_min_loc = min(h_loc);
                c_max_loc = max(c_vals_loc, [], 1);
                B_penalty_loc_minus = obj.sigma*c_max_loc(2)/h_min_loc*[B_penalty_ref_1, B_penalty_ref_2;
                                                                B_penalty_ref_2.', B_penalty_ref_3];
                B_penalty_loc = -B_penalty_loc_minus;
                B_penalty_vals(penalty_triplet_iterator:(penalty_triplet_iterator + 4*dof^2 - 1)) = B_penalty_loc(:);
                B_penalty_rows(penalty_triplet_iterator:(penalty_triplet_iterator + 4*dof^2 - 1)) = B_loc_row_idx(:);
                B_penalty_cols(penalty_triplet_iterator:(penalty_triplet_iterator + 4*dof^2 - 1)) = B_loc_col_idx(:);
                penalty_triplet_iterator = penalty_triplet_iterator + 4*dof^2;
            end

            % update boundary matrices
            [lower_boundary_element_idx, upper_boundary_element_idx] = obj.initial_mesh.getBoundaryElementIdx();
            [~, ~, elements] = obj.initial_mesh.getPet();
            c_vals_loc = [obj.wave_speed_cell{2}(lower_boundary_element_idx,1), obj.wave_speed_cell{1}(lower_boundary_element_idx,1);
                              obj.wave_speed_cell{2}(upper_boundary_element_idx,1), obj.wave_speed_cell{1}(upper_boundary_element_idx,1)];
            system_struct.B_bound(elements(lower_boundary_element_idx,:), elements(lower_boundary_element_idx,:)) = system_struct.B_bound(elements(lower_boundary_element_idx,:), elements(lower_boundary_element_idx,:))/c_vals_loc(1,2);
            system_struct.B_bound(elements(upper_boundary_element_idx,:), elements(upper_boundary_element_idx,:)) = system_struct.B_bound(elements(upper_boundary_element_idx,:), elements(upper_boundary_element_idx,:))/c_vals_loc(2,2);
            
            % rebuild sparse matrices
            n = size(obj.initial_mesh.nodes,1);
            system_struct.A = sparse(A_rows, A_cols, A_vals, n, n);
            system_struct.B_flux_int = sparse(B_flux_rows, B_flux_cols, B_flux_vals, n, n);
            system_struct.B_penalty_int = sparse(B_penalty_rows, B_penalty_cols, B_penalty_vals, n, n);
            
            system_struct.B = system_struct.A - system_struct.B_flux_int + system_struct.B_penalty_int + system_struct.B_bound;

            obj.precompute_matrix_resonator_masks(system_struct);

            % dissect sparse matrices
            [A_rows, A_cols, A_vals] = find(system_struct.A);
            [B_flux_rows, B_flux_cols, B_flux_vals] = find(system_struct.B_flux_int);
            [B_penalty_rows, B_penalty_cols, B_penalty_vals] = find(system_struct.B_penalty_int);
            
            % iterate over resonators
            for i = 0:size(obj.initial_mesh.resonators_matrix, 1)

                % current and old c_value in resonator
                res_el = find(obj.initial_mesh.element_idx_to_resonator_idx_map == i);
                c_res_initial = obj.wave_speed_cell{1}(res_el(1,1));
                c_multiplier = 1/c_res_initial;

                % find triplet indices, such that both rows and cols are in the resonator
                A_triplets_in_resonator_idx = obj.matrix_resonator_masks(i+1).A_masks;
                B_flux_triplets_in_resonator_idx = obj.matrix_resonator_masks(i+1).B_flux_masks;
                B_penalty_triplets_in_resonator_idx = obj.matrix_resonator_masks(i+1).B_penalty_masks;

                % update all entries for which both dofs (row and col) are in the resonator 
                A_vals(A_triplets_in_resonator_idx) = A_vals(A_triplets_in_resonator_idx)*c_multiplier;
                B_flux_vals(B_flux_triplets_in_resonator_idx) = B_flux_vals(B_flux_triplets_in_resonator_idx)*c_multiplier;
                B_penalty_vals(B_penalty_triplets_in_resonator_idx) = B_penalty_vals(B_penalty_triplets_in_resonator_idx)*c_multiplier;
            end

            % rebuild sparse matrices
            n = size(obj.initial_mesh.nodes,1);
            system_struct.A = sparse(A_rows, A_cols, A_vals, n, n);
            system_struct.B_flux_int = sparse(B_flux_rows, B_flux_cols, B_flux_vals, n, n);
            system_struct.B_penalty_int = sparse(B_penalty_rows, B_penalty_cols, B_penalty_vals, n, n);
            
            system_struct.B = system_struct.A - system_struct.B_flux_int + system_struct.B_penalty_int + system_struct.B_bound;

            obj.neutered_matrix_struct = system_struct;

        end

        function obj = assemble_initial_matrices(obj)
            % collects system matrices into a struct for later time adaptation
            % with the coefficients at time t0

            t0 = obj.pde_data.initial_time;

            % initialize c_vals
            [nodes_quad, ~, elements_quad] = obj.quadrature_mesh.getPet();
            c_evaluation_nodes = nodes_quad(elements_quad);
            if ~obj.pde_data.wave_speed_is_continuous
                c_evaluation_nodes(:, end) = c_evaluation_nodes(:, end) - 1e-14; 
                c_evaluation_nodes(:, 1) = c_evaluation_nodes(:, 1) + 1e-14; 
            end
            c_vals_temp = obj.pde_data.wave_speed_coeff_fun(c_evaluation_nodes, t0);
            obj.wave_speed_cell = {c_vals_temp, c_vals_temp};

            system_struct = obj.assemble_matrices_at_time(t0);

            f_values = obj.pde_data.rhs_fun(nodes_quad(elements_quad), t0);
            system_struct.load_vector = fem1d.loadVector1D(obj.initial_mesh.nodes, obj.initial_mesh.elements, f_values);
            
            obj.initial_matrix_struct = system_struct;

            if obj.matrix_update_type == "piecewise-const-coefficient-in-space"
                obj.neuter_initial_matrix_struct();
            end
        end

        function system_struct = assemble_matrices_at_time(obj, t)
            % similarly to assemble_initial_matrices assembles all required matrices into a struct
            % but with the coefficients a t time t
            % changes matrices at the boundary depending on the boundary condition

            % initialize c_vals
            [nodes_quad, ~, elements_quad] = obj.quadrature_mesh.getPet();
            c_evaluation_nodes = nodes_quad(elements_quad);
            if ~obj.pde_data.wave_speed_is_continuous
                c_evaluation_nodes(:, end) = c_evaluation_nodes(:, end) - 1e-14; 
                c_evaluation_nodes(:, 1) = c_evaluation_nodes(:, 1) + 1e-14; 
            end
            c_vals = obj.pde_data.wave_speed_coeff_fun(c_evaluation_nodes, t);

            % get matrices
            [nodes, ~, elements] = obj.initial_mesh.getPet();
            A = fem1d.stiffnessMatrix1D(nodes, elements, c_vals);
            if isempty(obj.initial_matrix_struct)
                M = fem1d.massMatrix1D(nodes, elements, ones(size(c_vals)));
            else
                M = obj.initial_matrix_struct.M;
            end
            B_flux_int = dg1d.interiorFluxMatrix1D(nodes, elements, c_vals);
            B_flux_bound = dg1d.boundaryFluxMatrix1D(nodes, elements, c_vals);
            B_penalty_int = dg1d.interiorPenaltyMatrix1D(nodes, elements, c_vals, obj.sigma);
            B_penalty_bound = dg1d.boundaryPenaltyMatrix1D(nodes, elements, c_vals, obj.sigma);
            B_bound = B_penalty_bound - B_flux_bound;

            system_struct = struct("A", A, "M", M, "B_flux_int", B_flux_int, "B_bound", B_bound, ...
                                         "B_penalty_int", B_penalty_int);
                                         
            % change matrices based on b.c.
            [~, ~, elements] = obj.initial_mesh.getPet();
            [lower_boundary_element_idx, upper_boundary_element_idx] = obj.initial_mesh.getBoundaryElementIdx();
            boundary_element_idx = [lower_boundary_element_idx, upper_boundary_element_idx];
            for i = 1:2
                boundary_condition = obj.pde_data.boundary_conditions{i};
                switch boundary_condition.bc_type
                    case "dirichlet"
                    case "neumann"
                        system_struct.B_bound(elements(boundary_element_idx(i),:), elements(boundary_element_idx(i),:)) = 0;
                    case "transparent"
                        system_struct.B_bound(elements(boundary_element_idx(i),:), elements(boundary_element_idx(i),:)) = 0;
                    otherwise
                        error("the boundary condition type << " + bc_type + " >> has not been implemented")
                end
            end
            system_struct.B = system_struct.A - system_struct.B_flux_int + system_struct.B_penalty_int + system_struct.B_bound;
        end

        function obj = calculate_stable_stepsize(obj)
            % calculates a stable stepsize dt if one has not already been chosen using the initial mass and stiffness matrices
            if obj.dt > 0
                return
            end
            B_loc = obj.initial_matrix_struct.A - obj.initial_matrix_struct.B_flux_int + obj.initial_matrix_struct.B_penalty_int + obj.initial_matrix_struct.B_bound;
            M_loc = obj.initial_matrix_struct.M;
            lMax = eigs(B_loc,1);
            lMin = eigs(M_loc,1,0);
            obj.dt = sqrt(lMin/lMax)/obj.initial_mesh.dof;

            % correct for case where eigs doesn't converge
            if isnan(obj.dt)
                obj.dt = obj.initial_mesh.h_min/obj.initial_mesh.dof;
            end
        end

        function system_struct = setup_system(obj, current_time)
            % sets up system by introducing time dependent components and boundary conditions

            % initialize c_vals
            [nodes_quad, ~, elements_quad] = obj.quadrature_mesh.getPet();
            c_evaluation_nodes = nodes_quad(elements_quad);
            if ~obj.pde_data.wave_speed_is_continuous
                c_evaluation_nodes(:, end) = c_evaluation_nodes(:, end) - 1e-14; 
                c_evaluation_nodes(:, 1) = c_evaluation_nodes(:, 1) + 1e-14; 
            end
            c_vals_temp = obj.pde_data.wave_speed_coeff_fun(c_evaluation_nodes, current_time);
            obj.wave_speed_cell{2} = c_vals_temp;

            switch obj.matrix_update_type
                case "brute-force"
                    system_struct = obj.assemble_matrices_at_time(current_time);
                case "piecewise-const-coefficient-in-space"
                    system_struct = obj.update_matrices_for_piecewise_const_coefficient(current_time);
                otherwise
                    system_struct = obj.initial_matrix_struct;
            end

            system_struct.time = current_time;
            
            % collect load vector
            if obj.pde_data.rhs_is_time_independent
                system_struct.load_vector = obj.initial_matrix_struct.load_vector;
            else
                f_values = obj.pde_data.rhs_fun(obj.quadrature_mesh.nodes(obj.quadrature_mesh.elements), current_time);
                system_struct.load_vector = fem1d.loadVector1D(obj.initial_mesh.nodes, obj.initial_mesh.elements, f_values);
            end
            
            % boundary conditions
            system_struct = obj.impose_boundary_conditions(system_struct);
        end

        function system_struct = update_matrices_for_piecewise_const_coefficient(obj, current_time)
            % updates the matrices in system_struct based on the resonator distribution in the mesh 
            % by dividing through the initial coefficient and multiplying the current 
            % for connecting matrices needs to subtract old values at intersection nodes

            system_struct = obj.neutered_matrix_struct;

            % dissect sparse matrices
            [A_rows, A_cols, A_vals] = find(system_struct.A);
            [B_flux_rows, B_flux_cols, B_flux_vals] = find(system_struct.B_flux_int);
            [B_penalty_rows, B_penalty_cols, B_penalty_vals] = find(system_struct.B_penalty_int);
            
            % iterate over resonators
            for i = 0:size(obj.initial_mesh.resonators_matrix, 1)

                % current and old c_value in resonator
                res_el = find(obj.initial_mesh.element_idx_to_resonator_idx_map == i);
                c_res_current = obj.wave_speed_cell{2}(res_el(1,1));
                c_multiplier = c_res_current;

                % find triplet indices, such that both rows and cols are in the resonator
                A_triplets_in_resonator_idx = obj.matrix_resonator_masks(i+1).A_masks;
                B_flux_triplets_in_resonator_idx = obj.matrix_resonator_masks(i+1).B_flux_masks;
                B_penalty_triplets_in_resonator_idx = obj.matrix_resonator_masks(i+1).B_penalty_masks;

                % update all entries for which both dofs (row and col) are in the resonator 
                A_vals(A_triplets_in_resonator_idx) = A_vals(A_triplets_in_resonator_idx)*c_multiplier;
                B_flux_vals(B_flux_triplets_in_resonator_idx) = B_flux_vals(B_flux_triplets_in_resonator_idx)*c_multiplier;
                B_penalty_vals(B_penalty_triplets_in_resonator_idx) = B_penalty_vals(B_penalty_triplets_in_resonator_idx)*c_multiplier;
            end

            % recalculate interior resonator interface values for flux and penalty matrix (updating connecting information between resonators)
            res_boundary_upper_element_idx = obj.initial_mesh.getInteriorResonatorBoundaryFace();

            % extend triplets
            dof = obj.initial_mesh.dof;
            empty_triplet_slots = zeros(dof^2*4*length(res_boundary_upper_element_idx), 1);
            flux_triplet_iterator = length(B_flux_rows) + 1;
            penalty_triplet_iterator = length(B_penalty_rows) + 1;
            B_flux_rows = [B_flux_rows; empty_triplet_slots];
            B_flux_cols = [B_flux_cols; empty_triplet_slots];
            B_flux_vals = [B_flux_vals; empty_triplet_slots];
            B_penalty_rows = [B_penalty_rows; empty_triplet_slots];
            B_penalty_cols = [B_penalty_cols; empty_triplet_slots];
            B_penalty_vals = [B_penalty_vals; empty_triplet_slots];

            % build local matrices for B_flux (still need to multiply c, (1/h) 
            [phi_val, dphi_val, ~] = common.getShapeFunctionValueMatrix(dof);
            B_flux_ref_1 = 1*(phi_val(:,end)*dphi_val(:,end).' + dphi_val(:, end)*phi_val(:,end).');
            B_flux_ref_21 = 1*phi_val(:,end)*dphi_val(:,1).';
            B_flux_ref_22 = (-1)*dphi_val(:,end)*phi_val(:,1).';
            B_flux_ref_3 = (-1)*(phi_val(:,1)*dphi_val(:,1).' + dphi_val(:, 1)*phi_val(:,1).');

            B_penalty_ref_1 = (1)*phi_val(:, end)*phi_val(:,end).';
            B_penalty_ref_2 = (-1)*phi_val(:, end)*phi_val(:,1).';
            B_penalty_ref_3 = (1)*phi_val(:,1)*phi_val(:,1).';
            for k = res_boundary_upper_element_idx.'
            
                % initializations for both 
                elements_loc = obj.initial_mesh.elements(k-1:k,:);
                nodes_loc = obj.initial_mesh.nodes(elements_loc);
                h_loc = [abs(nodes_loc(1,1) - nodes_loc(1,end)); abs(nodes_loc(2,1) - nodes_loc(2,end))];
                c_vals_loc = [obj.wave_speed_cell{2}(k-1,end), obj.wave_speed_cell{1}(k-1,end);
                              obj.wave_speed_cell{2}(k,1), obj.wave_speed_cell{1}(k,1)];
                bordering_elements = [elements_loc(1,:), elements_loc(2,:)];
                B_loc_row_idx = repmat(bordering_elements.', 1, 2*dof);
                B_loc_col_idx = repmat(bordering_elements, 2*dof, 1);
                
                % flux update
                B_flux_loc_plus = [c_vals_loc(1,1)*B_flux_ref_1/h_loc(1), (c_vals_loc(2,1)*B_flux_ref_21/h_loc(2) + c_vals_loc(1,1)*B_flux_ref_22/h_loc(1));
                                    (c_vals_loc(2, 1)*B_flux_ref_21/h_loc(2) + c_vals_loc(1, 1)*B_flux_ref_22/h_loc(1)).', c_vals_loc(2, 1)*B_flux_ref_3/h_loc(2)];
                B_flux_loc = B_flux_loc_plus;
                B_flux_vals(flux_triplet_iterator:(flux_triplet_iterator + 4*dof^2 - 1)) = B_flux_loc(:);
                B_flux_rows(flux_triplet_iterator:(flux_triplet_iterator + 4*dof^2 - 1)) = B_loc_row_idx(:);
                B_flux_cols(flux_triplet_iterator:(flux_triplet_iterator + 4*dof^2 - 1)) = B_loc_col_idx(:);
                flux_triplet_iterator = flux_triplet_iterator + 4*dof^2;

                % penalty update
                h_min_loc = min(h_loc);
                c_max_loc = max(c_vals_loc, [], 1);
                B_penalty_loc_plus = obj.sigma*c_max_loc(1)/h_min_loc*[B_penalty_ref_1, B_penalty_ref_2;
                                                                B_penalty_ref_2.', B_penalty_ref_3];
                B_penalty_loc = B_penalty_loc_plus;
                B_penalty_vals(penalty_triplet_iterator:(penalty_triplet_iterator + 4*dof^2 - 1)) = B_penalty_loc(:);
                B_penalty_rows(penalty_triplet_iterator:(penalty_triplet_iterator + 4*dof^2 - 1)) = B_loc_row_idx(:);
                B_penalty_cols(penalty_triplet_iterator:(penalty_triplet_iterator + 4*dof^2 - 1)) = B_loc_col_idx(:);
                penalty_triplet_iterator = penalty_triplet_iterator + 4*dof^2;
            end

            % update boundary matrices
            [lower_boundary_element_idx, upper_boundary_element_idx] = obj.initial_mesh.getBoundaryElementIdx();
            [~, ~, elements] = obj.initial_mesh.getPet();
            c_vals_loc = [obj.wave_speed_cell{2}(lower_boundary_element_idx,1), obj.wave_speed_cell{1}(lower_boundary_element_idx,1);
                              obj.wave_speed_cell{2}(upper_boundary_element_idx,1), obj.wave_speed_cell{1}(upper_boundary_element_idx,1)];
            system_struct.B_bound(elements(lower_boundary_element_idx,:), elements(lower_boundary_element_idx,:)) = system_struct.B_bound(elements(lower_boundary_element_idx,:), elements(lower_boundary_element_idx,:))*c_vals_loc(1,1);
            system_struct.B_bound(elements(upper_boundary_element_idx,:), elements(upper_boundary_element_idx,:)) = system_struct.B_bound(elements(upper_boundary_element_idx,:), elements(upper_boundary_element_idx,:))*c_vals_loc(2,1);
            
            % rebuild sparse matrices
            n = size(obj.initial_mesh.nodes,1);
            system_struct.A = sparse(A_rows, A_cols, A_vals, n, n);
            system_struct.B_flux_int = sparse(B_flux_rows, B_flux_cols, B_flux_vals, n, n);
            system_struct.B_penalty_int = sparse(B_penalty_rows, B_penalty_cols, B_penalty_vals, n, n);
            
            system_struct.B = system_struct.A - system_struct.B_flux_int + system_struct.B_penalty_int + system_struct.B_bound;

            % DEBUG
            % system_struct_debug = obj.assemble_matrices_at_time(current_time);
            % errors = [norm(full(system_struct_debug.A - system_struct.A));
            %           norm(full(system_struct_debug.B_bound - system_struct.B_bound));
            %           norm(full(system_struct_debug.B_flux_int - system_struct.B_flux_int));
            %           norm(full(system_struct_debug.B_penalty_int - system_struct.B_penalty_int));];
            % s = 1;
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
            
            % collect b.c vectors and matrix
            dirichlet_vector = dg1d.dirichletbcVector1D(nodes, elements, c_vals, g_vals,obj. sigma);
            neumann_vector = dg1d.neumannbcVector1D(nodes, elements, c_vals, g_vals);
            transparent_bc_matrix = sparse(size(nodes, 1), size(nodes, 1)); 

            for i = 1:2
                boundary_condition = obj.pde_data.boundary_conditions{i};
                switch boundary_condition.bc_type
                    case "dirichlet"
                        neumann_vector(elements(boundary_element_idx(i),:)) = 0;
                    case "neumann"
                        dirichlet_vector(elements(boundary_element_idx(i),:)) = 0;
                    case "transparent"
                        dirichlet_vector(elements(boundary_element_idx(i),:)) = 0;
                        neumann_vector(elements(boundary_element_idx(i),:)) = 0;
                        if boundary_condition.bc_location == obj.initial_mesh.lower_interval_bound
                            transparent_bc_matrix(1, 1) = sqrt(c_vals(1, 1));
                        elseif boundary_condition.bc_location == obj.initial_mesh.upper_interval_bound
                            transparent_bc_matrix(end, end) = sqrt(c_vals(end, end));
                        else
                            error(sprintf('transparent b.c. neither at the upper nor lower interval bound, but at: %d', boundary_condition.bc_location));
                        end
                    otherwise
                        error("the boundary condition type << " + bc_type + " >> has not been implemented")
                end
            end
            system_struct.R = transparent_bc_matrix; 
            system_struct.load_vector = system_struct.load_vector + dirichlet_vector + neumann_vector;
        end

        function obj = implement_initial_conditions(obj)
            % introduces the already known initial displacement into solution and calculates the solution at the time t0 + dt
            % using a taylor-expansion
            u0 = obj.pde_data.initial_displacement(obj.initial_mesh.nodes);
            v0 = obj.pde_data.initial_velocity(obj.initial_mesh.nodes); 
            system_struct = obj.setup_system(obj.time_vector(1));
            obj.solution(:, 1) = u0;
            system_matrix = system_struct.M;
            system_rhs = obj.dt^2/2*(system_struct.load_vector - system_struct.B*u0 - system_struct.R*v0);
            obj.solution(:, 2) = u0 + obj.dt*v0 + obj.solve_system(system_matrix, system_rhs);
        end

        function obj = setup_iteration(obj)
            % main preparation function for the iteration, initializes and preallocates everything needed for the actual 
            % iteration
            obj.initialize_quadrature();
            obj.setup_settings();
            obj.assemble_initial_matrices();
            obj.calculate_stable_stepsize();

            obj.time_vector = obj.pde_data.initial_time:obj.dt:obj.pde_data.final_time;
            obj.solution = zeros(size(obj.initial_mesh.nodes,1), size(obj.time_vector, 2));
            obj.implement_initial_conditions();
        end

        function obj = leap_frog_leap(obj)
            % main iteration, applies leapfrog time integration
            M = obj.initial_matrix_struct.M;
            numiter = length(obj.time_vector)-1;
            
            for i = 2:numiter
                system_struct = obj.setup_system(obj.time_vector(i));
                system_matrix = (M + obj.dt/2*system_struct.R);
                system_rhs = obj.dt^2*system_struct.load_vector + ...
                            (2*system_struct.M - obj.dt^2*system_struct.B)*obj.solution(:, i) - ...
                            (system_struct.M - obj.dt/2*system_struct.R)*obj.solution(:, i-1);
                % obj.solution(:, i+1) = obj.solve_system(system_matrix, system_rhs);
                obj.solution(:, i+1) = system_matrix\system_rhs;
            end
        end

        function obj = run(obj)
            obj.setup_iteration();
            obj.leap_frog_leap();
        end
    end
end