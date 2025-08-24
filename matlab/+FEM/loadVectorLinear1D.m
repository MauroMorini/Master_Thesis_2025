function b = loadVectorLinear1D(x, T, f)
    % calculate the nx1 load vector for finite elements solution of the
    % poisson equation in 1D given a function handle f, a grid x and a
    % connectivity matrix T
    
    % number of elements
    nEl = size(T, 1);       
    
    n = length(x);
    b = zeros(n, 1);
    
    % iterate over elements
    for i = 1:nEl
        % elementwise stepsize 
        h = abs(x(T(i,end)) - x(T(i,1)));
        
        % middle of element
        m = (x(T(i,end)) + x(T(i,1)))/2;
        
        % calculate elementwise load vector using simpson rule
        bK = zeros(2,1);
        bK(1) = (1/3)*(h/2)*(f(x(T(i,1))) + 4*(f(m)*1/2));
        bK(2) = (1/3)*(h/2)*(f(x(T(i,end))) + 4*(f(m)*1/2));
        
        % calculate elementwise load vector using trapezoidal rule
        % bK(1) = (h/2)*f(x(i));
        % bK(2) = (h/2)*f(x(i+1));
        
        % g1 = @(x) f(m + x*h/2).*(1-x)/2;
        % g2 = @(x) f(m + x*h/2).*(1+x)/2;
        % bK(1) = h/2*integral(g1, -1, 1);
        % bK(2) = h/2*integral(g2, -1, 1);
        
        % assemble global load vector
        b(T(i,:)) = b(T(i,:)) + bK;
    end
end