function [x, history] = admm_tv_fb_3D(dataset, admm_opts, convergence_para)

% this fuction solves the following problem total variation problem via ADMM:
%========================================================================
%    minimize 0.5*|| Ax - b ||_2^2 + \lambda || Fx ||_1          (a)
%========================================================================
% 
% Input:
% ------ 
%  dataset
%       * A : Forward projection matrix when applied to ODF yields
%       Polefigure. It has the dimension [m n]
%       * b : Experimental Polefigure of size [m 1]
%
%  admm_opts
%       * rho : is the augmented Lagrangian parameter. assign rho > 0.
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
%   norms at each iteration. it also consists of pfsnr i.e. Polefigure signal to noise
%   ratio (SNR) gain after each iteration in [db].
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_start = tic;

%convergence parameters
MAX_ITER = convergence_para.MAX_ITER; 
ABSTOL   = convergence_para.ABSTOL;
RELTOL   = convergence_para.RELTOL;

%data encompassed by eq (a) above 
A  = dataset.A;
b  = dataset.b;
%gt = dataset.gt;
%b_tilde = dataset.b_sim;

%penalty/regularization parameters
lambda = admm_opts.lambda;
rho    = admm_opts.rho;
alpha  = admm_opts.alpha;

[m, n] = size(A);

%3D diff operator
N=round(n^(1/3));
e = ones(n,1);
F1 = (1/6)*spdiags([-e e], 0:1, n,n);
F2 = (1/6)*spdiags([-e e], [0 N], n,n);
F3 = (1/6)*spdiags([-e e], [0 round(N*N)], n,n);

F1b = -F1';
F2b = -F2';
F3b = -F3';

FtF1 = F1'*F1;
FtF2 = F2'*F2;
FtF3 = F3'*F3;

FtF1b = F1b'*F1b;
FtF2b = F2b'*F2b;
FtF3b = F3b'*F3b;

FtF = FtF1 + FtF2 + FtF3 + FtF1b + FtF2b + FtF3b;
                                                  
%W = spdiags(wi, [0], m, m);
% save the backprojected result
Atb = A'*b;

%ADMM solver
x   = zeros(n,1);
z1  = zeros(n,1);
z2  = z1;
z3  = z1;

z1b = z1;
z2b = z1;
z3b = z1;

u1  = zeros(n,1);
u2  = zeros(n,1);
u3  = zeros(n,1);

u1b = zeros(n,1);
u2b = zeros(n,1);
u3b = zeros(n,1);

H   = A'*A+rho*FtF;

if strcmpi(admm_opts.sub_optimal_method, 'inv')
    inv_H = inv(H);
elseif strcmpi(admm_opts.sub_optimal_method, 'lu')
    H = A'*A+rho*FtF;
    L = chol(H, 'lower');
    L = sparse(L);
    U = sparse(L');
else
    H = H;
end

ncub=(N-1)/2;
mask=cubicBoundryPrior(N, (ncub/2)+4);
mask=mask(:);

if strcmpi (convergence_para.disp, 'True')
    fprintf('\n------------------------------------ Begin iterative cost minimization of current PF ------------------------------ \n')
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'costFunc','SNR');
    fprintf('------------------------------------------------------------------------------------------------------------------- \n')
end

for k = 1:MAX_ITER

    % x - suboptimal 
    q = Atb + rho*(F1'*(z1 - u1)     +F2'*(z2 - u2)    +F3'*(z3-u3) ...
                   + F1b'*(z1b - u1b)+ F2b'*(z2b - u2b) +F3b'*(z3b-u3b));
                   
    if strcmpi(admm_opts.sub_optimal_method, 'inv')
        x = inv_H*q;
    elseif strcmpi(admm_opts.sub_optimal_method, 'lu')
        x = U \ (L \ q);
    else
        [x, ~]=cgs(H, q, 1e-5, 200, [], [], x);
    end
    
    ind    = find(x<0);
    x(ind) = 0;
    
    if strcmpi(convergence_para.mask, 'T')
        x = x.*mask;
    end
    
     
    z1old = z1;
    F1x_hat = alpha*F1*x +(1-alpha)*z1old;
    z1 = shrinkage(F1x_hat + u1, lambda/rho);
    
    z2old = z2;
    F2x_hat = alpha*F2*x +(1-alpha)*z2old;
    z2 = shrinkage(F2x_hat + u2, lambda/rho);
    
    z3old = z3;
    F3x_hat = alpha*F3*x +(1-alpha)*z3old;
    z3 = shrinkage(F3x_hat + u3, lambda/rho);
    
    z1bold = z1b;
    F1bx_hat = alpha*F1b*x +(1-alpha)*z1bold;
    z1b = shrinkage(F1bx_hat + u1b, lambda/rho);
    
    z2bold = z2b;
    F2bx_hat = alpha*F2b*x +(1-alpha)*z2bold;
    z2b = shrinkage(F2bx_hat + u2b, lambda/rho);
    
    z3bold = z3b;
    F3bx_hat = alpha*F3b*x +(1-alpha)*z3bold;
    z3b = shrinkage(F3bx_hat + u3b, lambda/rho);
    
    
    % u - updates
    u1 = u1 + F1x_hat - z1;
    u2 = u2 + F2x_hat - z2;
    u3 = u3 + F3x_hat - z3;

    u1b = u1b + F1bx_hat - z1b;
    u2b = u2b + F2bx_hat - z2b;
    u3b = u3b + F3bx_hat - z3b;
    
    %primal & dual residuals
    history.cost(k)  = costFunc(A, b, lambda, x, z1, z2, z3, z1b, z2b, z3b);

    r1=F1*x-z1;
    r2=F2*x-z2;
    r3=F3*x-z3;
    
    r1b=F1b*x-z1b;
    r2b=F2b*x-z2b;
    r3b=F3b*x-z3b;
    
    history.r_norm(k)  = norm([r1' r2' r3' r1b' r2b' r3b']);
    
    s1=-rho*F1'*(z1 - z1old);
    s2=-rho*F2'*(z2 - z2old);
    s3=-rho*F3'*(z3 - z3old);
    
    s1b=-rho*F1b'*(z1b - z1bold);
    s2b=-rho*F2b'*(z2b - z2bold);
    s3b=-rho*F3b'*(z3b - z3bold);
    history.s_norm(k)  = norm([s1' s2' s3' s1b' s2b' s3b']);

    history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm([z1' z2' z3' z1b' z2b' z3b']));
    history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*[u1' u2' u3' u1b' u2b' u3b']);
    history.snr(k)     = pf_snr(b, x, A);
    
  % reporting     
  if strcmpi (convergence_para.disp, 'True') 
      fprintf('%3d\t%10.6f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n', k, ...
            history.r_norm(k), history.eps_pri(k), ...
            history.s_norm(k), history.eps_dual(k), history.cost(k), ...
            history.snr(k));
  end

    
  % termination checks
    if (history.r_norm(k) < history.eps_pri(k) && ...
       history.s_norm(k) < history.eps_dual(k))
         break;
    end

end
toc(t_start);

end

function p = costFunc(A, b, lambda, x, z1, z2, z3, z1b, z2b, z3b)
    p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z1,1) + lambda*norm(z2,1) ...
        + lambda*norm(z3,1) + lambda*norm(z1b,1) + lambda*norm(z2b,1) ...
        + lambda*norm(z3b,1));
end
