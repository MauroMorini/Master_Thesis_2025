classdef MeshIntervalDG1d < handle & matlab.mixin.Copyable
    % MeshIntervalDG1d is a class for creating and managing a 1D mesh
    % for discontinuous Galerkin methods.
    %
    % Properties:
    %   nodes - Array of node coordinates in the mesh.
    %   elements - Connectivity matrix defining the mesh elements.
    %   boundary_interface_node_idx - Indices of boundary nodes in the mesh.
    %   dof - Number of degrees of freedom per element.
    %   h_max - Maximum mesh size.
    %   h_min - Minimum mesh size.
    %   lower_interval_bound - Lower bound of the mesh interval.
    %   upper_interval_bound - Upper bound of the mesh interval.
    %   element_interface_nodes - Coordinates of the element interface nodes.
    %   face_node_to_element_map - Mapping of internal face nodes to elements.
    %
    % Methods:
    %   MeshIntervalDG1d - Constructor for the MeshIntervalDG1d class.
    %   updatePet - Updates the mesh properties and connectivity based on the
    %               current element interface nodes.
    %   getPet - Returns the current nodes, elements, and boundary node indices.
    %   plotMesh - Plots the mesh and highlights the boundary nodes.
    properties
        nodes                       % (num_nodes, 1) node matrix (with duplicates)
        elements                    % (num_el, dof) connectivity matrix
        boundary_node_idx           % (1,2) index vector with 
        dof
        h_max
        h_min
        lower_interval_bound
        upper_interval_bound
        element_interface_nodes
        face_node_to_element_map
        lower_boundary_element_idx
        upper_boundary_element_idx
        is_resonator_mesh = false;
        resonators_matrix                   % (num_res, 2) matrix which contains the intervals which correspond with a resonator 
                                            % each row has to be sorted from small to big
        element_idx_to_resonator_idx_map    % (num_el, 1) int vector containing number in {0,...,num_res} corresponding to resonator indices
                                            % and 0 corresponding to the background
    end

    methods
        function obj = MeshIntervalDG1d(varargin)
            switch (nargin)
                case 0
                    % default: equidistant mesh on [0,1] with mesh size 0.1
                    obj.lower_interval_bound = 0;
                    obj.upper_interval_bound = 1;
                    obj.h_max = 0.1;
                    obj.h_min = obj.h_max/100;    
                    obj.dof = 2;                
                case 2
                    % inputs: ([lower_interval_bound,upper_interval_bound],[h_max, h_min])
                    intBounds = varargin{1};
                    Hbounds = varargin{2};
                    obj.lower_interval_bound = intBounds(1);
                    obj.upper_interval_bound = intBounds(2);
                    obj.h_max = Hbounds(1);
                    obj.h_min = Hbounds(2);
                    obj.dof = 2;
                otherwise
                    error("Constructor has wrong number of inputs")
            end
            assert(obj.h_max >= 5*obj.h_min, "h_min should be small enough, maximally h_max/5")
            obj.element_interface_nodes = (obj.lower_interval_bound:obj.h_max:(obj.upper_interval_bound-obj.h_min))'; 
            obj.element_interface_nodes(end+1) = obj.upper_interval_bound;
            while abs(obj.element_interface_nodes(end) - obj.element_interface_nodes(end-1)) > obj.h_max
                obj.element_interface_nodes = [obj.element_interface_nodes(1:end-1); (obj.element_interface_nodes(end-1)+obj.upper_interval_bound)/2; obj.upper_interval_bound];
            end
            obj = updatePet(obj);
        end
% UPDATE, GET, SET -------------------------------------------------------------------------------------------------------------------------
        function obj = updatePet(obj)
            % updates connectivity matrix and edge matrix and sets nodes
            obj.element_interface_nodes = sort(obj.element_interface_nodes);        % possibly remove sorting
            if obj.element_interface_nodes(1) ~= obj.lower_interval_bound || obj.element_interface_nodes(end) ~= obj.upper_interval_bound
                error("boundary points of element_interface_nodes don't coincide with properties of class")
            end
            % check that h_min and h_max are maintained
            D = abs(diff(obj.element_interface_nodes));
            corrSize = all(obj.h_min-10*eps <= D & D <= obj.h_max+10*eps);
            if ~corrSize
                warning("There are elements which are too small or too big")
            end
            obj.setNodes();
            
            % connectivity matrix
            num_nodes = length(obj.nodes);
            N = length(obj.element_interface_nodes);
            elements_temp = 1:num_nodes;
            elements_temp = reshape(elements_temp, obj.dof, []);
            obj.elements = elements_temp.';

            % face_node_to_element_map
            obj.face_node_to_element_map = [NaN, 1; (1:N-2).', (2:N-1).'; N-1, NaN];
            
            % find boundary element indices
            obj.lower_boundary_element_idx = 1;
            obj.upper_boundary_element_idx = size(obj.elements, 1);

            % update resonator map
            obj.fill_element_idx_to_resonator_idx_map();
        end

        function obj = setNodes(obj)
            % creates interior nodes based on Gauss-Lobatto
            % and updates object
            N = length(obj.element_interface_nodes);
            [quad_nodes, ~] = common.getLobattoQuadrature(obj.dof);
            element_meshsizes = diff(obj.element_interface_nodes);
            element_midpoints = obj.element_interface_nodes(1:end-1) + element_meshsizes/2;
            nodes_matrix_extended = zeros(N-1, obj.dof);
            nodes_matrix_extended(:, 1) = obj.element_interface_nodes(1:end-1);
            nodes_matrix_extended(:, 2:end) = element_midpoints + quad_nodes(2:end).*element_meshsizes/2;
            nodes_matrix_extended = nodes_matrix_extended.';
            obj.nodes = nodes_matrix_extended(:);
            num_nodes = length(obj.nodes);
            obj.boundary_node_idx = [1, num_nodes];
        end

        function [nodes, boundary_interface_node_idx, elements] = getPet(obj)
            nodes = obj.nodes;
            elements = obj.elements;
            boundary_interface_node_idx = obj.boundary_node_idx;
            return 
        end

        function [lower_boundary_element_idx, upper_boundary_element_idx] = getBoundaryElementIdx(obj)
            lower_boundary_element_idx = obj.lower_boundary_element_idx;
            upper_boundary_element_idx = obj.upper_boundary_element_idx;
        end
% UTILS---------------------------------------------------------------------------------------------------

% REFINE ------------------------------------------------------------------------------------------------
        function obj = refineElementsByFact(obj, elements_to_be_refined_idx, refine_factor)
            % refines a list of given elements by a factor
            arguments (Input)
                obj
                elements_to_be_refined_idx      % (N, 1) index vector (could also boolean vector be) N<=num_el
                refine_factor double = 2
            end
            arguments (Output)
                obj
            end

            % calculate element sizes of to be refined elements and check h_min condition
            h_loc_temp = diff(obj.element_interface_nodes);
            h_loc_temp = h_loc_temp(elements_to_be_refined_idx);
            h_loc = h_loc_temp/refine_factor;
            non_refinable_element_idx = find(h_loc < obj.h_min - eps);
            refIdx = elements_to_be_refined_idx;
            if ~isempty(non_refinable_element_idx)
                warning("some elements were to small to be refined: " + non_refinable_element_idx)
                refIdx(non_refinable_element_idx) = [];
                h_loc(non_refinable_element_idx) = [];
            end

            % add aditional faces
            pLoc = zeros(length(refIdx), refine_factor);
            pLoc(:,1) = obj.element_interface_nodes(refIdx);
            for i = 1:(refine_factor-1)
                pLoc(:,i+1) = pLoc(:,i) + h_loc;
            end
            untouched_elements = obj.element_interface_nodes;
            untouched_elements(refIdx) = [];
            total_face_nodes = [untouched_elements;pLoc(:)];
            obj.element_interface_nodes = total_face_nodes;
        end

        function refinableTIdx = getRefinableElements(obj, refFactor)
            % returns index of elements which can be refined by the factor
            % given
            arguments (Input)
                obj
                refFactor double = 2 
            end
            arguments (Output)
                refinableTIdx
            end
            D = abs(diff(obj.element_interface_nodes));
            refinableTIdx = find(D/refFactor >= obj.h_min);
        end

        function obj = refineAll(obj, refine_factor)
            % refines all refinable elements and sets h_min and h_max smaller depending on the refinement factor
            obj.h_min = obj.h_min/refine_factor;
            ref_idx = obj.getRefinableElements(refine_factor);
            obj.refineElementsByFact(ref_idx, refine_factor);
            obj.h_max = obj.h_max/refine_factor;
        end

        function KIdx = findElementAt(obj, x)
            % finds element containing point x, if x is in two elements it
            % returns an array
            % returns index of element
            assert(obj.lower_interval_bound <= x & x <= obj.upper_interval_bound, "Point x is outside of the domain")
            K = [obj.element_interface_nodes(obj.elements(:,1)), obj.element_interface_nodes(obj.elements(:,2))];
            containsX = (K(:,1) <= x & x <= K(:,2)) | (K(:,2) <= x & x <= K(:,1));
            KIdx = find(containsX);
        end

% CREATE MESH ROUTINES --------------------------------------------------------------------------------------------------
        function obj = createUniformMesh(obj, h)
            % creates a uniform mesh for a given meshsize h in (h_min, h_max)
            arguments (Input)
                obj
                h
            end
            assert(h < obj.h_max && h > obj.h_min, "h needs to be in (h_min, h_max)")
            obj.element_interface_nodes = (obj.lower_interval_bound:h:obj.upper_interval_bound).';
        end
        function obj = createRngMesh(obj, seed)
            % Function creates a pseudorandom mesh with random inner points for
            % a given seed
            arguments (Input)
                obj
                seed double = 1 
            end
            arguments (Output)
                obj
            end
            r = rng(seed, "twister");
            h = rand(1);
            h = obj.h_min*h + obj.h_max*(1-h);
            pLoc = [obj.lower_interval_bound];
            i = 1;
            while pLoc(i)+h < obj.upper_interval_bound - obj.h_max
                i = i + 1;
                h = rand(1);
                h = obj.h_min*h + obj.h_max*(1-h);
                pLoc(i) = pLoc(i-1)+h;
            end
            while pLoc(i)+obj.h_max/2 < obj.upper_interval_bound-obj.h_min
                i = i+1;
                pLoc(i) = pLoc(i-1)+obj.h_max/2;
            end
            obj.element_interface_nodes = [pLoc'; obj.upper_interval_bound];
        end

        function obj = buildResonatorMesh(obj, resonators_mat, stepsizes_vector)
            % build a mesh with uniform stepsizes except some intervalls coresponding to resonators
            % which are more refined
            arguments (Input)
                obj
                resonators_mat (:,2)        % (num_res, 2) matrix which contains the intervals which correspond with a resonator 
                                            % each row has to be sorted from small to big
                stepsizes_vector (1,2)      % stepsize vector which has the bigger stepsize (background) in the first entry and 
                                            % the resonator stepsize in the second entry (will be sorted though)
            end
            arguments (Output)
                obj
            end

            if isempty(resonators_mat)
                return
            end

            % check that inputs are admissible
            if any(resonators_mat(:,1) >=  resonators_mat(:,2))
                error("the resonator pairs are not arranged properly, they should be of the form: resonators_mat = [r1, R1; r2, R2;...] where ri < Ri")
            end
            if max(resonators_mat, [], 'all') > obj.upper_interval_bound || min(resonators_mat, [], 'all') < obj.lower_interval_bound
                error("a resonator is place outside of the domain")
            end
            stepsizes_vector = sort(stepsizes_vector,'descend');
            h_background = stepsizes_vector(1);
            h_resonator = stepsizes_vector(2);
            if  h_resonator <= obj.h_min || h_background > obj.h_max
                error("the stepsizes are either too big or too small in comparison to h_max, h_min")
            end
            if h_resonator >= min(abs(resonators_mat(:,1)- resonators_mat(:,2)))
                error("the resonator stepsize is bigger or equal the length of one of the resonators")
            end

            % sort resonators 
            resonators_mat_sorted = resonators_mat;
            [resonators_mat_sorted(:, 1), sort_idx] = sort(resonators_mat(:,1));
            resonators_mat_sorted(:,2) = resonators_mat(sort_idx, 2);

            % build resonator mesh by setting face nodes
            last_node_in_chain = obj.lower_interval_bound;
            local_faces = [];
            num_res = size(resonators_mat_sorted,1);
            for resonator_idx = 1:num_res

                % build up to lower resonator boundary node
                resonators = resonators_mat_sorted(resonator_idx, :);
                local_faces = [local_faces, last_node_in_chain:h_background:resonators(1)-obj.h_min];
                while abs(resonators(1) - local_faces(end-1)) > h_background
                    local_faces = [local_faces, (local_faces(end)+resonators(1))/2];
                end
                last_node_in_chain = resonators(1);

                % build up to upper resonator boundary node
                local_faces = [local_faces, last_node_in_chain:h_resonator:resonators(2)-obj.h_min];
                while abs(resonators(2) - local_faces(end)) > h_resonator
                    local_faces = [local_faces, (local_faces(end)+resonators(2))/2];
                end          
                last_node_in_chain = resonators(2);    
            end

            % build up to upper boundary node
            local_faces = [local_faces, last_node_in_chain:h_background:obj.upper_interval_bound-obj.h_min];
            while abs(obj.upper_interval_bound - local_faces(end)) > h_background
                local_faces = [local_faces, (local_faces(end)+obj.upper_interval_bound)/2];
            end

            obj.element_interface_nodes = [local_faces';obj.upper_interval_bound];
            obj.resonators_matrix = resonators_mat_sorted;
            obj.is_resonator_mesh = true;
        end

        function obj = fill_element_idx_to_resonator_idx_map(obj)
            obj.element_idx_to_resonator_idx_map = zeros(size(obj.element_interface_nodes));
            if ~obj.is_resonator_mesh
                obj.element_idx_to_resonator_idx_map(end) = [];
                return
            end
            num_res = size(obj.resonators_matrix, 1);
            for i = 1:num_res
                resonator_loc = obj.resonators_matrix(i, :);
                res_idx = find(resonator_loc(1) <= obj.element_interface_nodes & obj.element_interface_nodes < resonator_loc(2));
                obj.element_idx_to_resonator_idx_map(res_idx) = i;
            end
            obj.element_idx_to_resonator_idx_map(end) = [];
        end

% VISUALIZE ------------------------------------------------------------------------------------------------------
        function f = plotMesh(obj, f)
            if nargin == 1
                f = figure;
            end
            figure(f);            
            plot(obj.element_interface_nodes,0,'b.','MarkerSize',10)
            hold on
            plot(obj.nodes(obj.boundary_node_idx), [0,0], 'rx','MarkerSize',10)
            hold off
        end

        function f = plotDGsol(obj, uh, f)
        % plots dg solution visualizing the discontinuity
            if nargin == 2
                f = figure;
            end
            % initialization
            figure(f); 
            element_faces = obj.element_interface_nodes;
            num_element_plot_points = obj.dof*10;
            bary_weights = common.calculateBarycentricWeights(obj.nodes, obj.elements);
            num_el = size(obj.elements, 1);

            hold on
            for k = 1:num_el
                element_nodes = obj.nodes(obj.elements(k,:));
                plot_nodes_loc = linspace(element_nodes(1), element_nodes(end), num_element_plot_points).';
                bary_weights_loc = bary_weights(obj.elements(k,:));
                Phi_loc = common.evaluateLagrangeBarycentric(plot_nodes_loc, bary_weights_loc, element_nodes);
                uh_vals_loc = sum(Phi_loc .* uh(obj.elements(k,:)), 1);
                plot(plot_nodes_loc, uh_vals_loc, 'Color', 'k', 'LineWidth', 2)
                plot([plot_nodes_loc(1), plot_nodes_loc(end)], [uh_vals_loc(1), uh_vals_loc(end)], 'o', 'Color', 'k', 'LineWidth', 2)
                xlabel("x")
                ylabel("uh")
            end
            hold off
        end
    end

    methods (Access = protected)
        function cp = copyElement(obj)
            cp = copyElement@matlab.mixin.Copyable(obj);
        end
    end
end