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