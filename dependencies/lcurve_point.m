function [x_coord, y_coord] = lcurve_point(x, b, A)

    [m, n] = size(A);

    %  3D diff operator
    N  = round(n^(1/3));
    e  = ones(n,1);
    F1 = (1/6)*spdiags([-e e], 0:1, n,n);
    F2 = (1/6)*spdiags([-e e], [0 N], n,n);
    F3 = (1/6)*spdiags([-e e], [0 round(N*N)], n,n);

    F1b = -F1';
    F2b = -F2';
    F3b = -F3';

    x_coord = sum((A*x - b).^2);
    x_coord = log(x_coord);
    y_coord = norm(F1*x,1) + norm(F2*x,1) ...
            + norm(F3*x,1) + norm(F1b*x,1) ...
            + norm(F2b*x,1)+ norm(F3b*x,1);
    y_coord = log(y_coord);
   
end