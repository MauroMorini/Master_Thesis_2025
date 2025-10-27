classdef SIPGWaveSolver
    properties
        initial_mesh        % MeshIntervalDG1d object
    end

    methods
        function obj = SIPGWaveSolver(inputArg1)
            obj.Property1 = inputArg1^2;
            disp(obj.Property1);
            
        end
    end
end