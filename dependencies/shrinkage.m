function z = shrinkage(x, kappa)

% This function solves the following optimization problem:
% ========================================================================
% minimize (\lambda |x_i| + (\rho/2)(x_i-v_i)
% ========================================================================
% x+ :=S_kappa (v_i);
% where, the softholding operator S is defined as
%	            { a ? kappa, a > kappa
% S_kappa (a) = { 0 	   , |a|<= kappa
%	            { a + kappa, a < -kappa
% where, kappa = lambda/rho

    z = max( 0, x - kappa ) - max( 0, -x - kappa );
end
