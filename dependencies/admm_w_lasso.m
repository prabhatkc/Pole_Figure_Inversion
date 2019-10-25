function [x, history] = admm_w_lasso(dataset, admm_opts, convergence_para)

% admm_w_lasso solves the following problem via ADMM:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    minimize 0.5*Ax - b)'*W*(Ax-b)+ \lambda || Fx ||_1          (a)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Input:
% ------ 
%  dataset
%       * A : Forward projection matrix when applied to ODF yields
%       Polefigure. It has the dimension [m n]
%       * b : Polefigure of size [m 1]
%       * gt: if applicable (i.e. for simulated calculation) denotes ground
%       truth ODF of size [n 1]
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
%   norms at each iteration. In the cases where the Grouth truth is known, we
%   also include RMSE (Root Mean Squared Error) & PSNR (peak signal to noise
%   ratio) value at each iteration. 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_start = tic;

%convergence parameters
MAX_ITER = convergence_para.MAX_ITER; 
ABSTOL   = convergence_para.ABSTOL;
RELTOL   = convergence_para.RELTOL;

%data encompassed by eq (a) above 
A       = dataset.A;
b       = dataset.b;
if isfield(dataset, 'gt')
    gt      = dataset.gt;
    b_tilde = dataset.b_sim;
end

%penalty/regularization parameters
lambda = admm_opts.lambda;
rho    = admm_opts.rho;
alpha  = admm_opts.alpha;
[m, n] = size(A);


% save a matrix-vector multiply (Backprojected result)
wi  = zeros(m, 1);
ind = find(b~=0);
wi(ind) = 1./b(ind);
W = spdiags(wi, [0], m, m);
AWb = A'*W*b;

% ADMM solver
x = zeros(n,1);
z = zeros(n,1);
u = zeros(n,1);

% cache the factorization
[L U] = factor(A, W, rho);

if strcmp(convergence_para.mask(1), 'T')
    N    = round(n^(1/3));
    ncub = (N-1)/2;
    mask = cubicBoundryPrior(N, (ncub/2)+4);
    mask =  mask(:);
end

fprintf('\n------------------------------------ Begin iterative cost minimization of current PF ------------------------------ \n')
if isfield(dataset, 'gt')
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'costFunc','PSNR', 'ODF err');
else
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'costFunc','PSNR');
end 
fprintf('------------------------------------------------------------------------------------------------------------------- \n')



for k = 1:MAX_ITER
  
    % x - suboptimal
    q = AWb + rho*(z - u);    % temporary value
    
    % if fat
    x = q/rho - (A'*W*(U \ ( L \ (A*q) )))/rho^2;
   
    ind    = find(x<0);
    x(ind) = 0;
    
    if strcmp(convergence_para.mask(1), 'T')
        x = x.*mask;
    end
    % z - suboptimal
    zold = z;
    x_hat = alpha*x + (1 - alpha)*zold;
    z = shrinkage(x_hat + u, lambda/rho);

    % u - update
    u = u + (x_hat - z);

    % diagnostics, reporting, termination checks
    history.cost(k)  = costFunc(A, b, W, lambda, x, z);

    history.r_norm(k)  = norm(x - z);
    history.s_norm(k)  = norm(-rho*(z - zold));

    history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
    history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*rho*norm(u);
    
    if isfield(dataset, 'gt')
        history.odf_err(k) = odf_err(gt, x);
        %history.isnr(k)    = isnr(b, b_tilde, x, A);
        history.psnr(k) = pf_psnr(A*x, b_tilde);


        fprintf('%3d\t%10.6f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n', k, ...
            history.r_norm(k), history.eps_pri(k), ...
            history.s_norm(k), history.eps_dual(k), history.cost(k), ...
            history.psnr(k), history.odf_err(k));
        
        if (k < 50)
            if (history.r_norm(k) < history.eps_pri(k) && ...
                history.s_norm(k) < history.eps_dual(k))
                break;
            end  
        else
            if ((history.r_norm(k) < history.eps_pri(k) && ...
                history.s_norm(k) < history.eps_dual(k)) ...
                || (history.odf_err(k) > history.odf_err(k-1)))
                break;
            end
        end
    else 
        history.snr(k)     = pfsnr(b, x, A);
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n', k, ...
                history.r_norm(k), history.eps_pri(k), ...
                history.s_norm(k), history.eps_dual(k), history.cost(k), ...
                history.snr(k));
        if (history.r_norm(k) < history.eps_pri(k) && ...
            history.s_norm(k) < history.eps_dual(k))
            break;
        end
    end

end

end

function p = costFunc(A, b, W, lambda, x, z)
    p = ( 1/2*sum((A*x - b)'*W*(A*x-b)) + lambda*norm(z,1) );
end

function [L U] = factor(A, W, rho)
    [m, n] = size(A);
    % if fat
    L = chol( speye(m) + 1/rho*(A*A'*W), 'lower' );
    % force matlab to recognize the upper / lower triangular structure
    L = sparse(L);
    U = sparse(L');
end

