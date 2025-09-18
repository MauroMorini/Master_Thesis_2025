function b = loadVectorQuadratic1D_v0(x, T, f)
        % calculate the nx1 load vector for finite elements solution of the
        % poisson equation in 1D given a function handle f, a grid x and a
        % connectivity matrix T for quadratic elements
        
        % number of elements
        nEl = size(T, 1);       
        
        n = length(x);
        b = zeros(n, 1);
        
        % iterate over elements
        for i = 1:nEl
        
            % element 
            K = x(T(i,:));
            
            % elementwise stepsize 
            h = abs(K(end) - K(1));
        
            % shape function
            N = {@(xi) 1/2*(xi.^2 - xi); @(xi) 1-xi.^2; @(xi) 1/2*(xi.^2 + xi)};
        
            % calculate elementwise load vector using simpson rule
            w = [1, 4, 1]/3;
            bK = zeros(3,1);
            y = [-1, 0, 1];
            
            for p = 1:3
                bK(p) = (N{p}(y).*f(K))*w.';
            end
            bK = bK*h/2;
        
            % assemble global load vector
            b(T(i,:)) = b(T(i,:)) + bK;
        end
        end