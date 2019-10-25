function [err]=polefig_err(gt1d, recon1d)
% 
% polefigure error for the reconstructed odf X_recon
% w.r.t the ground truth X_true is defined as:
%
% ========================================================================
% ||X_recon - X_true ||_2/ ||X_true||_2
% ========================================================================
%
%
	err = norm(gt1d-recon1d)/norm(gt1d);
end