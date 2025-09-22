classdef LobattoQuadTests < matlab.unittest.TestCase
    addpath()
    methods (TestClassSetup)
        % Shared setup for the entire test class
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        % Test methods

        function testLobattoNodes(testCase)
            tol = 1e-9;
            realNodes = {[-1, 1],[-1, 0, 1],[-1, -1/sqrt(5), 0, 1/sqrt(5), 1]};
            for i = length(realNodes)
                quad_nodes = common.getLobattoQuadratureNodes(i+1);
                testCase.assertLessThan(norm(realNodes{i}-quad_nodes,2), tol)
            end
        end
    end

end