% author: Mauro Morini  
% last modified: 08.05.24
classdef FEM1D
    % collection of FE methods in 1D

    properties
        Property1
    end

    methods(Static)
        % author: Mauro Morini
        % last modified: 22.03.2024
        function A = fem1d.stiffnessMatrix1D(x, T, c)
        % calculate the nxn stiffness matrix using a connectivity matrix T and a
        % possibly non equidistant grid x of size n for 2 or 3 degrees of freedom
        % A = [(c*phi'_i, phi'_j)_L^2]i,j
        %
        % Inputs :
        % T : (nEl, DoF) connectivity matrix
        % x : (1, n) point vector
        % c : function handle 
        %
        % Outputs:
        % A : (n,n) sparse stiffness matrix
                
        % number of elements 
        nEl = size(T, 1); 
        
        n = length(x);
        A = sparse(n,n);
        
        % iterate over all elements
        for i = 1:nEl
            
            % element
            K = x(T(i,:));
        
            % get element matrix 
            AK = fem1d.stiffnessElementMatrix1D(K, c);
        
            % assembling of stiffness matrix
            A(T(i,:), T(i,:)) = A(T(i,:), T(i,:)) + AK;
        
        end
        end

        % author: Mauro Morini
        % last modified: 24.03.24
        function AK = stiffnessElementMatrix1D(K, c)
        % calculates element matrix for a given element K on the grid, a function
        % handle c for 3 degrees of freedom (i.e. r=2). Code can be adapted for
        % higher grade polynomials
        %
        % Inputs:
        % K : (1,Dof) element point vector
        % 
        % Output:
        % AK : (Dof,Dof) local element matrix
        
        % Initializations
        Dof = size(K, 2);
        if Dof == 1
            K = K';
            Dof = size(K, 2);
        end
        AK = zeros(Dof);
        
        % initialize variables depending on Dof 
        switch Dof
            case 2
                % first derivative shape functions 
                N = {@(xi) -(1/2)*ones(size(xi)), @(xi) (1/2)*ones(size(xi))};
        
                % quadrature 
                y = [-1, 1];   
                w = [1, 1];     
                B = [N{1}(y); N{2}(y)];
            case 3
                % first derivative shape functions 
                N = {@(xi)1/2*(2*xi-1), @(xi) -2*xi, @(xi) 1/2*(1+2*xi)};
        
                % quadrature 
                y = [-1, 0, 1];   
                w = [1, 4, 1]/3;
                B = [N{1}(y); N{2}(y); N{3}(y)];
            otherwise
                error("local stiffness matrix has not been implemented for " + Dof + " degrees of freedom")
        end
        
        % element length
        h = abs(K(end) - K(1));
        
        % calculate functional values for K of c    
        cVal = c(K);
        
        % calculate integral for (A_K)_pq 
        for p = 1:Dof
            for q = 1:p
                AK(p,q) = (B(p,:).*B(q,:).*cVal)*w.';
                AK(q,p) = AK(p,q);
            end
        end
        AK = AK*(2/h);
        end
        
        % author: Mauro Morini  
        % last modified: 17.03.24
        function M = massMatrix1D(x, t, c)
        % calculates global mass matrix for linear or quadratic 
        % FE in 1D M(i,j) = int_Omega phi_i*phi_j*c
        % Inputs : 
        % x :  (1, nP) pointvector 
        % t : nEx3 connectivity matrix with elements in rows
        % c : function handle 
        %
        % Outputs : 
        % M : global stiffness matrix nPxnP
        
        % decide if quadratic or linear FE
        DoF = size(t, 2);
        
        nP = length(x);
        nE = size(t, 1);
        M = sparse(nP, nP);
        
        % iterate over elements
        for i = 1:nE
            
            % element
            K = x(t(i,:));
        
            % get element matrix 
            Mloc = fem1d.massElementMatrix1D(K, c);
        
            % Assembling
            M(t(i, :), t(i, :)) = M(t(i, :), t(i, :)) + Mloc;
        end
        end
    
        % author: Mauro Morini
        % last modified: 24.03.24
        function AK = massElementMatrix1D(K, c)
        % calculates element matrix for a given element K on the grid, a function
        % handle c for 3 degrees of freedom (i.e. r=3). Code can be adapted for
        % higher grade polynomials
        %
        % Inputs:
        % K : (1,Dof) element point vector
        % 
        % Output:
        % AK : (Dof,Dof) local element matrix
        
        % Initializations
        Dof = size(K, 2);
        if Dof == 1
            K = K';
            Dof = size(K, 2);
        end
        AK = zeros(Dof);
        
        % initialize variables depending on Dof 
        switch Dof
            case 2
                % shape functions 
                N = {@(xi) (1-xi)/2, @(xi) (1+xi)/2};
        
                % quadrature 
                quadK = [K(1), (K(1)+K(2))/2, K(2)];
                y = [-1, 0, 1];   
                w = [1, 4, 1]/3;     
                B = [N{1}(y); N{2}(y)];
            case 3
                % shape functions 
                N = {@(xi)1/2*(xi.^2 - xi), @(xi) 1 - xi.^2, @(xi) 1/2*(xi + xi.^2)};
        
                % quadrature 
                quadK = [K(1), (K(1)+K(2))/2, K(2), (K(2)+K(3))/2, K(3)];
                y = [-1, -0.5, 0, 0.5, 1];   
                w = [7, 32, 12, 32, 7]*(2/90);
                B = [N{1}(y); N{2}(y); N{3}(y)];
            otherwise
                error("local mass matrix has not been implemented for " + Dof + " degrees of freedom")
        end
        
        % element length
        h = abs(K(end) - K(1));
        
        % calculate functional values for K of c    
        cVal = c(quadK);
        
        % calculate integral for (A_K)_pq 
        for p = 1:Dof
            for q = 1:p
                AK(p,q) = (B(p,:).*B(q,:).*cVal)*w.';
                AK(q,p) = AK(p,q);
            end
        end
        AK = AK*(h/2);
        end
        
        % author: Mauro Morini
        % last modified: 08.12.24
        function b = loadVector1D(x, T, f)
            Dof = size(T,2);

            switch Dof
                case 1
                    error("connectivity matrix needs to have more than 1 node in element")
                case 2
                    b = fem1d.loadVectorLinear1D(x,T,f);
                case 3
                    b = fem1d.loadVectorQuadratic1D(x,T,f);
                otherwise
                    error("load vector for " + Dof + " : degrees of freedom per element has not been implemented")
            end
        end

        % author: Mauro Morini
        % last modified: 17.10.23
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
        
        % author: Mauro Morini
        % last modified: 28.10.23
        function b = loadVectorQuadratic1D(x, T, f)
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
            w = [1, 4, 1];
            bK = zeros(3,1);
            y = [-1, 0, 1];
            
            for p = 1:3
                bK(p) = (N{p}(y).*f(K))*w.';
            end
            bK = bK*h/6;
        
            % assemble global load vector
            b(T(i,:)) = b(T(i,:)) + bK;
        end
        end

        % author: Mauro Morini
        % last modified: 03.11.24
        function M = projMassMatrix1D(p1,t1,p2,t2)
        % V1 = span{phi_j^(1)}, V2 = span{phi_j^(2)} are the span of hat basis functions over nodes p1, p2
        % calculates mass matrix Mij = (phi_i^(2), phi_j^(1))_L2
        % note for projection here V2 is the projected space

        % Initializations
        N1 = length(p1);
        N2 = length(p2);
        M = sparse(N2, N1);

        % % functions
        % N0 = @(x) (1-x)/2;
        % N1 = @(x) (1+x)/2;
        % FKInv = @(x,h,m) 2*(x-m)/h;
        
        % iterate over elements of V1
        for i = 1:size(t1, 1)
            K1 = p1(t1(i,:));

            % h1 = abs(K1(1)-K1(end));
            % m1 = (K1(1)+K1(end))/2;

            % find elements in V2 which have a non empty intersection with
            % K
            p2El = [p2(t2(:,1)), p2(t2(:,end))];
            idxEmpty = p2El(:,2) <= K1(1) | K1(end) <= p2El(:,1);        % elements with empty intersection
            idxEl = find(~idxEmpty);
            for j = 1:length(idxEl)
                K2 = p2(t2(idxEl(j),:));

                % h2 = abs(K2(1)-K2(end));
                % m2 = (K2(1) + K2(end))/2;
                % 
                % % Intersection element
                % KInt = [max(K1(1),K2(1)), min(K1(2),K2(end))];
                % hInt = abs(KInt(1)-KInt(2));
                % mInt = (KInt(1) + KInt(2))/2;
                % f = @(x) [N0(FKInv(x,h2,m2))*N0(FKInv(x,h1,m1)), N0(FKInv(x,h2,m2))*N1(FKInv(x,h1,m1));
                %             N1(FKInv(x,h2,m2))*N0(FKInv(x,h1,m1)), N1(FKInv(x,h2,m2))*N1(FKInv(x,h1,m1))];
                % Mloc = hInt/6*(f(KInt(1)) + 4*f(mInt) + f(KInt(2)));

                Mloc = fem1d.projMassElementMatrix1D(K1, K2);
                M(t2(idxEl(j),:), t1(i,:)) = M(t2(idxEl(j),:), t1(i,:)) + Mloc;
            end          
        end
        end
        
        % author: Mauro Morini
        % last modified: 08.12.24
        function Mloc = projMassElementMatrix1D(K1, K2)
            Dof = length(K1);
            Mloc = zeros(Dof);
            FKInv = @(x,h,m) 2*(x-m)/h;

            h1 = abs(K1(1)-K1(end));
            m1 = (K1(1)+K1(end))/2;
            h2 = abs(K2(1)-K2(end));
            m2 = (K2(1)+K2(end))/2;

            % find intersection Element
            K = [max(K1(1),K2(1)), min(K1(end),K2(end))];
            h = abs(K(1)-K(2));
            m = (K(1)+K(2))/2;

            % initialize variables depending on Dof 
            switch Dof
                case 2
                    % shape functions 
                    N = {@(xi) (1-xi)/2, @(xi) (1+xi)/2};
            
                    % quadrature 
                    quadK = [K(1), m, K(2)];
                    y = [-1, 0, 1];   
                    w = [1, 4, 1]/3;     
                case 3
                    % shape functions 
                    N = {@(xi)1/2*(xi.^2 - xi), @(xi) 1 - xi.^2, @(xi) 1/2*(xi + xi.^2)};
            
                    % quadrature 
                    quadK = [K(1), (K(1)+m)/2, m, (K(2)+m)/2, K(2)];
                    y = [-1, -0.5, 0, 0.5, 1];   
                    w = [7, 32, 12, 32, 7]*(2/90);
                    B1 = [N{1}(y); N{2}(y); N{3}(y)];
                otherwise
                    error("local mass matrix has not been implemented for " + Dof + " degrees of freedom")
            end
            B1 = zeros(Dof, length(quadK));
            B2 = zeros(Dof, length(quadK));
            for i = 1:Dof
                B1(i,:) = N{i}(FKInv(quadK,h1,m1)); 
                B2(i,:) = N{i}(FKInv(quadK,h2,m2)); 
            end
            % calculate integral for (A_K)_pq 
            for p = 1:Dof
                for q = 1:Dof
                    Mloc(p,q) = (B1(q,:).*B2(p,:))*w.';
                end
            end
            Mloc = Mloc*(h/2);
            end
        
        % author: Mauro Morini
        % last modified: 08.12.24
        function [L2err, H1err] = errors1D(T, x, uh, dudx, u)
            Dof = size(T,2);

            switch Dof
                case 1
                    error("connectivity matrix needs to have more than 1 node in element")
                case 2
                    [L2err, H1err] = fem1d.errorsLinear1D(T, x, uh, dudx, u);
                case 3
                    [L2err, H1err] = fem1d.errorsQuad1D(T, x, uh, dudx, u);
                otherwise
                    error("load vector for " + Dof + " : degrees of freedom per element has not been implemented")
            end
        end

        % author: Mauro Morini
        % last modified: 24.03.24
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
        
        % author: Mauro Morini
        % last modified: 28.10.23
        function [L2err, H1err] = errorsQuad1D(T, x, uh, dudx, u)
        % calculates the error in L^2 and H^1 norm for quadratic FE with
        % uh: numerical solution
        % u: exact solution
        % dudx: exact 1-st derivative of u
        
        r = size(T, 2);
        
        L2err = 0;
        H1err = 0;
        
        % weights QF
        w = [1/3 4/3 1/3];
        
        % nodes of QF
        y = [-1, 0, 1];
        
        % shape functions
        N = {@(xi) 1/2*(xi.^2 - xi); @(xi) 1-xi.^2; @(xi) 1/2*(xi.^2 + xi)};
        dN = {@(xi)1/2*(2*xi-1), @(xi) -2*xi, @(xi) 1/2*(1+2*xi)};
        
        % iterate over elements
        for i = 1:size(T, 1)
        
            % element
            K = x(T(i, :));
        
            % element length
            h = abs(K(1) - K(end));
        
            % middle of element
            m = (K(end) + K(1))/2;
        
            % element map
            F = @(xi) m + xi*h/2;
            
            uhVal = 0;
            duhVal = 0;
            for p = 1:r
                uhVal = uhVal + uh(T(i, p))*N{p}(y);
                duhVal = duhVal + uh(T(i, p))*dN{p}(y);
            end
        
            L2err = L2err + h/2*real(u(F(y)) - uhVal).^2*w.';
            H1err = H1err + h/2*real(dudx(F(y)) - 2/h*duhVal).^2*w.';
          
        end
        H1err = sqrt(H1err + L2err);
        L2err = sqrt(L2err);
        
        end
        
        % author: Mauro Morini
        % last modified: 03.11.24
        function L2err = L2ProjectionErrorLinear(u1, p1, u2, p2)
            % calculate L2 norm of u1 - u2 where Pi is a projection onto
            % the space generated by linear basis functions on p2
            % p1 and p2 need to be on the same domain with the same
            % boundary points
            %
            % Inputs: 
            % u1: (N1,1) solution vector
            % u2: (N2,1) solution vector
            % p1: (N1,1) point vector of V1
            % p2: (N2,1) point vector of V2
            %
            % Outputs:
            % L2err: scalar L2 error

            % create span{V1,V2}
            p = sort([p1; p2]);
            p = unique(p);
            N = length(p);
            t = [(1:N-1)', (2:N)'];
            u1Interp = interp1(p1,u1,p);
            u2Interp = interp1(p2,u2,p);

            L2err = 0;
            for i = 1:size(t,1)
                h = abs(p(t(i,1))-p(t(i,2)));
                L2err = L2err + (h/2)*sum((u1Interp(t(i,:)) - u2Interp(t(i,:))).^2);
            end
            L2err = sqrt(L2err);
        end
        
    end
end