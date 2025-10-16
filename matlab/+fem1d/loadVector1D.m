function load_vector = loadVector1D(nodes, elements, f_values)
    % calculates the load vector for finite elements (without b.c)
    arguments (Input)
        nodes           % (num_nodes, 1) node vector
        elements        % (num_el, dof) connectivity matrix
        f_values        % (num_el, num_quad) values of forcing term at nodes
    end
    arguments (Output)
        load_vector     % (num_nodes, 1) load vector 
    end

    % initializations
    nEl = size(elements, 1); 
    num_nodes = length(nodes);
    dof = size(elements, 2);
    num_quad = size(f_values, 2);

    % preallocation
    load_vector = zeros(num_nodes, 1);
    
    % collect quadrature information
    [phi_val, ~, quad_weights] = common.getShapeFunctionValueMatrix(dof, num_quad);

    % iterate over elements
    for k = 1:nEl
        % elementwise stepsize 
        h = abs(nodes(elements(k,end)) - nodes(elements(k,1)));
        
        loc_load_vector = zeros(dof,1);
        for p = 1:dof
            loc_load_vector(p) = (h/2)*phi_val(p,:).*f_values(k,:)*quad_weights.';
        end
        load_vector(elements(k,:)) = load_vector(elements(k,:)) + loc_load_vector;
    end
end