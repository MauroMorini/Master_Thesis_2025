% author: Mauro Morini
% last modified: 25.11.24
classdef WaveSolution
    % small object containing information of a certain numerical solution at a
    % time T and it's current mesh 

    properties
        t               % Time at which numerical sol is recorded
        uh               % (nP, 1) point array
        mesh            % Mesh1D object 
    end

    methods
        function obj = WaveSolution(mesh, t, uh)
            obj.mesh = mesh;
            obj.t = t;
            obj.uh = uh;
        end

        function [mesh, t, uh] = getSolution(obj)
            mesh = obj.mesh;
            t = obj.t;
            uh = obj.uh;
        end
    end
end