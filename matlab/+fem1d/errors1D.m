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