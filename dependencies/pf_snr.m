function[ratio] = pf_snr(y, x_hat, A)
	ratio = -20*log10(norm(y - A*x_hat)/norm(y));
end
