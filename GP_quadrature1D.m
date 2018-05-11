
function [gp,gpw]=GP_quadrature1D(ngp);  % gauss sampling points & weights for given the number of points (ngp)

if ngp == 1
    gp = [0];
    gpw = [2];

elseif ngp == 2
    gp = [-1/sqrt(3) 1/sqrt(3)];
    gpw = [1 1];

elseif ngp == 3
    gp = [-sqrt(0.6) 0 sqrt(0.6)];
    gpw = [5/9 8/9 5/9];
end
end