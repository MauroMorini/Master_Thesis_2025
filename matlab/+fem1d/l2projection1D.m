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
        u_projected (:,:)
    end
    
    if ~isNested
        % TODO implement non-nested case first create union space and project solution onto that one
    end

    

end