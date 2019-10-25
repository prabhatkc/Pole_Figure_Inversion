function [x, history] = admm_tv_3D(dataset, admm_opts, convergence_para)

% this fuction solve the following problem total variation problem via ADMM:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    minimize 0.5*|| Ax - b ||_2^2 + \lambda || Fx ||_1          (a)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Input:
% ------ 
%  dataset
%       * A : Forward projection matrix when applied to ODF yields
%       Polefigure. It has the dimension [m n]
%       * b : Polefigure of size [m 1]
%       * gt: denotes ground truth ODF of size [n 1]
%  admm_opts
%       * rho : is the augmented Lagrangian parameter. One must not assign a negative
%       value for rho. 
%       While determining value for rho we would suggest you to perform multiple
%       experiments starting with rho as 1. Subsequently maximize or minimize rho 
%       by a factor of 10 as you deem to see a fitting solution. 
%       * alpha: is the over-relaxation parameter (typical values for alpha are
%       between 1.0 and 1.8).
%
%
% Output:
% ------
%   * x       : denotes inverted solution and is returned in a vector form of size [n 1].
%   * history : is a structure that contains the cost function value, the primal and
%   dual residual norms, and the tolerances for the primal and dual residual
%   norms at each iteration. history also include, PSNR (Peak Signal to Noise Ratio)
%   & RMSE (Root Mean Squared Error) value at each iteration. 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t_start = tic;
QUIET    = 0;

%convergence parameters
MAX_ITER = convergence_para.MAX_ITER; 
ABSTOL   = convergence_para.ABSTOL;
RELTOL   = convergence_para.RELTOL;

%data encompassed by eq (a) above 
A  = dataset.A;
b  = dataset.b;
gt = dataset.gt;

%penalty/regularization parameters
lambda = admm_opts.lambda;
rho    = admm_opts.rho;
alpha  = admm_opts.alpha;

[m, n] = size(A);

%3D diff operator
N = round(n^(1/3));
e = ones(n,1);
F1 = spdiags([e -e], 0:1, n,n); %F1=0.5.*F1;
F2 = spdiags([e -e], [0 N], n,n); %F2=0.5.*F2;
F3 = spdiags([e -e], [0 round(N*N)], n,n); %F3=0.5.*F3;
FtF1 = F1'*F1;
FtF2 = F2'*F2;
FtF3 = F3'*F3;
FtF = FtF1+FtF2+FtF3;
                                                  
% save the backprojected result
Atb = A'*b;

%ADMM solver
x  = zeros(n,1);
z1 = zeros(n,1);
z2 = zeros(n,1);
z3 = zeros(n,1);
u1 = zeros(n,1);
u2 = zeros(n,1);
u3 = zeros(n,1);

H=A'*A+rho*FtF;

if strcmpi(admm_opts.submethod, 'inv')
    inv_H=inv(H);
else
    H = A'*A+rho*FtF;
    L = chol(H, 'lower');
    L = sparse(L);
    U = sparse(L');
end

ncub=(N-1)/2;
mask=cubicBoundryPrior(N, (ncub/2)+4);
mask=mask(:);

if ~QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'costFunc','pole err');
end


for k = 1:MAX_ITER

    % x - suboptimal 
    q = Atb + rho*(F1'*(z1 - u1)+F2'*(z2-u2)+F3'*(z3-u3)); 
    
    if strcmpi(admm_opts.submethod, 'inv')
        x = inv_H*q;
    elseif strcmpi(admm_opts.submethod, 'lu')
        x = U \ (L \ q);
    else
        [x, ~]=cgs(H, q, 1e-5, 200, L, U, x);
    end
    
    ind=find(x<0);
    x(ind)=0;
    x=x.*mask;
     
    % z - suboptimals
    z1old = z1;
    F1x_hat = alpha*F1*x +(1-alpha)*z1old;
    z1 = shrinkage(F1x_hat + u1, lambda/rho);
    
    z2old = z2;
    F2x_hat = alpha*F2*x +(1-alpha)*z2old;
    z2 = shrinkage(F2x_hat + u2, lambda/rho);
    
    z3old = z3;
    F3x_hat = alpha*F3*x +(1-alpha)*z3old;
    z3 = shrinkage(F3x_hat + u3, lambda/rho);
    
    % u - updates
    u1 = u1 + F1x_hat - z1;
    u2 = u2 + F2x_hat - z2;
    u3 = u3 + F3x_hat - z3;

    %primal & dual residuals
    r1=F1*x-z1;
    r2=F2*x-z2;
    r3=F3*x-z3;
    s1=-rho*F1'*(z1 - z1old);
    s2=-rho*F2'*(z2 - z2old);
    s3=-rho*F3'*(z3 - z3old);
    
    % diagnostics
    history.cost(k)    = costFunc(A, b, lambda, x, z1, z2, z3);
    history.pole_err(k)= odf_err(gt, x);
     
    history.r_norm(k)  = norm([r1' r2' r3']);
    history.s_norm(k)  = norm([s1' s2' s3']);
    history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm([z1' z2' z3']));
    history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*[u1' u2' u3']);
   
   % reporting 
    if ~QUIET
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n', k, ...
            history.r_norm(k), history.eps_pri(k), ...
            history.s_norm(k), history.eps_dual(k), history.cost(k), ...
            history.odf_err(k));
    end
    
    % termination checks
    if (history.r_norm(k) < history.eps_pri(k) && ...
       history.s_norm(k) < history.eps_dual(k))
         break;
    end

end

if ~QUIET
    toc(t_start);
end
end

function p = costFunc(A, b, lambda, x, z1, z2, z3)
    p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z1,1) + lambda*norm(z2,1) ...
        + lambda*norm(z3,1));
end
