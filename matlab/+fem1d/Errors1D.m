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

        function obj = run(obj)
            obj.generate_quadrature_mesh();
            obj.interpolate_sol_at_quad_mesh();
            obj.calculate_errors();
        end
    end
end