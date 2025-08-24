% author: Mauro Morini
% last modified: 26.11.24
classdef Mesh1d < handle
    % 

    properties
        a       % left interval border
        b       % right interval border
        Hmax    % maximal meshsize 
        Hmin    % minimal meshsize
        p       % point vector (N,1)
        e       % border points (indices)
        t       % connectivity matrix (N,r)
        r = 1   % degree of elements
    end

    methods
        function obj = Mesh1d(varargin) 
            % Constructor 
            % with exception of the default constructer, the inputs have to
            % be given in the following order
            % 
            % Inputs:
            % [a,b]: Intervall bounds
            % [Hmax, Hmin]: bounds for mesh size Hmin <= Hmax/5
            switch (nargin)
                case 0
                    % default: equidistant mesh on [0,1] with mesh size 0.1
                    obj.a = 0;
                    obj.b = 1;
                    obj.Hmax = 0.1;
                    obj.Hmin = obj.Hmax/100;                    
                case 2
                    % inputs: ([a,b],[Hmax, Hmin])
                    intBounds = varargin{1};
                    Hbounds = varargin{2};
                    obj.a = intBounds(1);
                    obj.b = intBounds(2);
                    obj.Hmax = Hbounds(1);
                    obj.Hmin = Hbounds(2);
                otherwise
                    error("Constructor has wrong number of inputs")
            end
            assert(obj.Hmax >= 5*obj.Hmin, "Hmin should be small enough, maximally Hmax/5")
            obj.p = (obj.a:obj.Hmax:(obj.b-obj.Hmin))'; 
            obj.p(end+1) = obj.b;
            while abs(obj.p(end) - obj.p(end-1)) > obj.Hmax
                obj.p = [obj.p(1:end-1); (obj.p(end-1)+obj.b)/2; obj.b];
            end
            obj = updatePet(obj);
        end

        function [p, e, t] = getPet(obj)
            % returns p,e,t values depending on obj.r the chosen element
            % degree
            p = obj.p;
            e = obj.e;
            t = obj.t;
            if obj.r == 1
                return
            end

            H = diff(p)/obj.r;
            for i = 1:(obj.r-1)
                pLoc = obj.p(1:end-1) + H*i;
                p = [p;pLoc];
            end
            p = sort(p);
            
            t = zeros(size(t,1), obj.r+1);
            % connectivity matrix
            for i = 1:obj.r+1
                t(:,i) = (i:obj.r:(length(p)-(obj.r+1-i)))';
            end
        end

        function obj = refine(obj, Href, refIdx)
            % refines called elements in refIdx to subelements of size
            % Href: Hmin <= Href < Hmax
            % refIdx: index array (not logical)
            if Href < obj.Hmin
                error("Refinement size cannot be smaller than Hmin")
            elseif Href >= obj.Hmax
                error("Refinement size cannot be bigger or equal to Hmax")
            end
            pLoc = obj.p;
            % iterate over to be refined elements add new points compatible
            % with Hmax and Hmin and add them to p then update
            for i = 1:length(refIdx)
                K = obj.p(obj.t(refIdx(i),:));
                xLoc = (K(1):Href:K(2))';           % new points in between element              
                if abs(xLoc(end) - K(2)) < obj.Hmin
                    xLoc = xLoc(1:end-1);
                end
                xLoc = xLoc(2:end);
                if isempty(xLoc)
                    warning("Href was too big, the element: " + refIdx(i) + " could not be refined")
                end
                pLoc(end+1:end+length(xLoc)) = xLoc;
            end
            obj.p = pLoc;
            obj = updatePet(obj);
        end

        function obj = refineByFact(obj, refIdx, refFact)
            % refines given elments by a factor refFact, default is 2, 
            %
            % Inputs: 
            % refIdx: index array (not logical)
            % refFact: natural number if not given default is 2
            if ~exist("refFact", 'var')
                refFact = 2;
            end

            pLoc = obj.p;
            % iterate over to be refined elements add new points compatible
            % with Hmax and Hmin and add them to p then update
            for i = 1:length(refIdx)
                K = obj.p(obj.t(refIdx(i),:));
                hK = abs(K(1) - K(2))/refFact;
                if hK < obj.Hmin - eps
                    warning("Href was too big, the element: " + refIdx(i) + " could not be refined")
                    continue
                end
                xLoc = linspace(K(1),K(2),refFact+1)';           % new points in between element              
                
                pLoc(end+1:end+(length(xLoc)-2)) = xLoc(2:end-1);
            end
            obj.p = pLoc;
            obj = updatePet(obj);

        end

        function [obj, isRemoved] = removePoints(obj, delIdx)
            % removes points in p with idx in delIdx from mesh if the
            % resulting bigger element does not suprass Hmax size, else
            % return the same object with a warning
            %
            % Inputs: 
            % delIdx: (m,1) idx vector with m < length(p) 
            
            isRemoved = true;
            N = length(obj.p);
            delete = false(N,1);
            delete(delIdx) = true;
            
            % Check that boundary nodes are not in delIdx
            assert(~delete(1) & ~delete(end), "The first and last point of p can't be eliminated");
            
            delete = delete(2:end-1);
            % check how big each element would be if points are removed
            pLeft = obj.p(1:end-2);
            pRight = obj.p(3:end);
            pDiff = abs(pLeft-pRight);

            % restrict delete vector to all points which when removed will
            % not leave a gap bigger than Hmax
            delete = (pDiff <= obj.Hmax + eps) & delete;
            if sum(delete) == 0
                warning("No point has been eliminated")
                isRemoved = false;
            end
            obj.p([false;delete;false]) = [];
            obj = updatePet(obj);
        end

        function [obj, isRemoved] = removeRand(obj, seed)
            % removes a random removable point in Mesh
            rng(seed)
            removableIdx = obj.getRemovablePoints();
            if isempty(removableIdx)
                isRemoved = false;
                return
            end
            rdmIdx = randi([1,length(removableIdx)],1,1);
            deleteIdx = removableIdx(rdmIdx);
            obj.p(deleteIdx) = [];
            obj = obj.updatePet();
        end

        function obj = shiftMesh(obj, Hshift)
            % shifts inner points by the param Hshift
            %
            % Inputs:
            % Hshift: scalar parameter in [Hmin, Hmax) \cup (-Hmax, -Hmin]
            if size(obj.p,1) == 2
                warning("There are only boundary points, no shift can happen")
                return
            end

            % force Hshift into interval [Hmin, Hmax]
            hMod = obj.Hmax - obj.Hmin;
            Hshift = (mod(abs(Hshift) - obj.Hmin, hMod) + obj.Hmin)*(-1)^(Hshift < 0);
            assert(obj.Hmax > abs(Hshift) && obj.Hmin <= abs(Hshift), "abs(Hshift) must be in [Hmin, Hmax)")

            pLoc = obj.p;

            % shift inner points to zero
            hInner = abs(obj.b - obj.a - 2*obj.Hmin) + realmin;
            displace = obj.a + obj.Hmin;
            pLoc(2:end-1) = mod(pLoc(2:end-1) - displace + Hshift, hInner) + displace;
            pLoc = sort(pLoc);

            % % tranform p onto unit interval
            % h = obj.b-obj.a;
            % pLoc = (pLoc - obj.a)/h;
            % 
            % % shift points
            % pLoc(2:end-1) = mod(pLoc(2:end-1) + Hshift, 1+realmin);
            % 
            % % transform back
            % pLoc = pLoc*h + obj.a;
            
            if false
                % removes points which got too close to the boundary points
                if abs(pLoc(1)-pLoc(2)) < obj.Hmin
                    pLoc = [pLoc(1); pLoc(3:end)];
                    %pLoc(2) = pLoc(2) + obj.Hmin;
                end
                if abs(pLoc(end)-pLoc(end-1)) < obj.Hmin
                    pLoc = [pLoc(1:end-2); pLoc(end)];
                    %pLoc(end-1) = pLoc(end-1) - obj.Hmin;
                end
                % add in points where the distance got too big
                if abs(pLoc(1)-pLoc(2)) > obj.Hmax
                    pLoc = [pLoc(1);(pLoc(2)+pLoc(1))/2 ;pLoc(2:end)];
                elseif abs(pLoc(3)-pLoc(2)) > obj.Hmax
                    pLoc = [pLoc(1:2);(pLoc(2)+pLoc(3))/2 ;pLoc(3:end)];
                end
                if abs(pLoc(end)-pLoc(end-1)) > obj.Hmax
                    pLoc = [pLoc(1:end-1);(pLoc(end)+pLoc(end-1))/2;pLoc(end)];
                elseif abs(pLoc(end-2)-pLoc(end-1)) > obj.Hmax
                    pLoc = [pLoc(1:end-2);(pLoc(end-2)+pLoc(end-1))/2;pLoc(end-1:end)];
                end
            end

            % check that Hmin and Hmax are maintained
            % D = abs(diff(pLoc));
            % corrSize = all(obj.Hmin-10*eps <= D & D <= obj.Hmax+10*eps);
            % if ~corrSize
            %     warning("There are elements which are too small or too big")
            % end
            
            obj.p = pLoc;
            obj = updatePet(obj);
        end

        function obj = createRngMesh(obj, seed)
            % Function creates a pseudorandom mesh with random inner points for
            % a given seed
            r = rng(seed, "twister");
            h = rand(1);
            h = obj.Hmin*h + obj.Hmax*(1-h);
            pLoc = [obj.a];
            i = 1;
            while pLoc(i)+h < obj.b - obj.Hmax
                i = i + 1;
                h = rand(1);
                h = obj.Hmin*h + obj.Hmax*(1-h);
                pLoc(i) = pLoc(i-1)+h;
            end
            while pLoc(i)+obj.Hmax/2 < obj.b-obj.Hmin
                i = i+1;
                pLoc(i) = pLoc(i-1)+obj.Hmax/2;
            end
            obj.p = [pLoc'; obj.b];
            obj = obj.updatePet();
        end

        function removablePIdx = getRemovablePoints(obj)
            % returns index vector of all points which can be removed
            % without breaking the Hmax condition (1,end) are never in it
            D = abs(diff(obj.p));
            D = D(2:end) + D(1:end-1);
            isRemovable = [false;D <= obj.Hmax + eps;false];
            removablePIdx = find(isRemovable);
        end

        function refinableTIdx = getRefinableElements(obj, refFactor)
            % returns index of elements which can be refined by the factor
            % given
            D = abs(obj.p(obj.t(:,1)) - obj.p(obj.t(:,2)));
            refinableTIdx = find(D/refFactor >= obj.Hmin);
        end

        function obj = updatePet(obj)
            % updates connectivity matrix and edge matrix 
            obj.p = sort(obj.p);
            if obj.p(1) ~= obj.a || obj.p(end) ~= obj.b
                error("boundary points of p don't coincide with properties of class")
            end
            % check that Hmin and Hmax are maintained
            D = abs(diff(obj.p));
            corrSize = all(obj.Hmin-10*eps <= D & D <= obj.Hmax+10*eps);
            if ~corrSize
                warning("There are elements which are too small or too big")
            end
            N = length(obj.p);
            obj.e = [obj.a, obj.b];
            obj.t = [(1:(N-1))', (2:N)'];
        end

        function KIdx = findElementAt(obj, x)
            % finds element containing point x, if x is in two elements it
            % returns an array
            % returns index of element
            assert(obj.a <= x & x <= obj.b, "Point x is outside of the domain")
            K = [obj.p(obj.t(:,1)), obj.p(obj.t(:,2))];
            containsX = (K(:,1) <= x & x <= K(:,2)) | (K(:,2) <= x & x <= K(:,1));
            KIdx = find(containsX);
        end

        function pIdx = findPointClosestTo(obj, x)
            % finds index of point in p closest to x
            assert(obj.a <= x & x <= obj.b, "Point x is outside of the domain")
            [~, pIdx] = min(abs(obj.p-x));
        end

        function f = plotMesh(obj, f)
            if nargin == 1
                f = figure;
            end
            figure(f);            
            plot(obj.p,0,'b.','MarkerSize',10)
            hold on
            plot(obj.e, [0,0], 'rx','MarkerSize',10)
            hold off
        end
    end
end