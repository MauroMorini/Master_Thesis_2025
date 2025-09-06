function [L2err, H1err] = errorsLinear1D(T, x, uh, dudx, u)
    % calculates the error in L^2 and H^1 norm for linear FE with
    % uh: numerical solution
    % u: exact solution
    % dudx: exact 1-st derivative of u
    
    L2err = 0;
    H1err = 0;
    
    % iterate over elements
    for i = 1:size(T, 1)
    
        % element
        K = x(T(i, :));
    
        % element midpoint
        m = (K(1) + K(end))/2;
    
        % element length
        h = abs(K(1) - K(end));
        
        % border values of uh on element
        uh1 = uh(T(i, 1));
        uh2 = uh(T(i, end));
        
        % calculate elementwise error of ||u-uh||_L^2
        % with the trapezoid rule and add it to the
        % sum
        L2err = L2err + (h/2)*abs(real(uh1 - u(K(1)))^2 + ...
            real(uh2 - u(K(end)))^2);
    
        % calculate elementwise error of ||dxdu - dxduh||_L^2 with trapezoid
        H1err = H1err + (h/2)*( (real(uh1 - uh2)/h + real(dudx(K(1))))^2 + ...
            (real(uh1 - uh2)/h + real(dudx(K(end))))^2);
    
        % ===================================================================
        % calculate elementwise error with simpson with the given integral
        % formula from the lecture (it doesn't work)
        % L2err = L2err + (h/6)*((u(K(1) - uh1)^2 + 4*(u(m) - (uh1 + uh2)/2)^2 ...
        %     + (u(K(end)) - uh2)^2 ));
    
        % testing with integral command
        % N1 = @(x) (1-x)/2;
        % N2 = @(x) (1+x)/2;
        % F1 = @(x) (u(m + x*h/2) - (uh1*N1(x) + uh2*N2(x))).^2;
        % F2 = @(x) (dudx(m + x*h/2)*h/2-(-uh1 + uh2)/2).^2;
        % L2err = L2err + h/2*integral(F1, -1, 1);
        % H1err = H1err + h/2*integral(F2, -1, 1);
        %====================================================================
       
    end
    H1err = sqrt(H1err + L2err);
    L2err = sqrt(L2err);
        
end