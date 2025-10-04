function u_projected = l2projection1D(nodes_origin, elements_origin, values_origin, nodes_target, elements_target, isNested)
    % calculates L^2-projection of a FEM-solution values_origin onto a target space 
    
    arguments (Input)
        nodes_origin (:,1)                 % (num_nodes_proj, 1) node vector for origin space
        elements_origin (:,:)
        values_origin (:,1)
        nodes_target (:,1)
        elements_target (:,:)
        isNested logical = false    
    end

    arguments (Output)
        u_projected (:,1)
    end

    % initialization
    transpose_mass_mat = false;
    values_origin_temp = values_origin;
    nodes_target_temp = nodes_target;
    elements_target_temp = elements_target;

    % find coarse and fine mesh
    if size(elements_origin,1) >= size(elements_target,1)
        nodes_coarse = nodes_target;
        elements_coarse = elements_target;
        nodes_fine = nodes_origin;
        elements_fine = elements_origin;       
    else
        nodes_coarse = nodes_origin;
        elements_coarse = elements_origin;
        nodes_fine = nodes_target;
        elements_fine = elements_target;
        transpose_mass_mat = true; 
    end
    
    if ~isNested
        % TODO implement non-nested case first create union space and project solution onto that one
    end

    % assemble matrices
    % M_mass_target = fem1d.massMatrix1D(nodes_target_temp, elements_target_temp, ones(size(nodes_target_temp)));
    M_mass_target = fem1d.nestedL2projectionMassMatrix1D(nodes_target_temp, elements_target_temp, nodes_target_temp, elements_target_temp);
    M_project = fem1d.nestedL2projectionMassMatrix1D(nodes_coarse, elements_coarse, nodes_fine, elements_fine);

    if transpose_mass_mat
        M_project = M_project.';
    end

    u_projected = M_mass_target\(M_project * values_origin_temp);