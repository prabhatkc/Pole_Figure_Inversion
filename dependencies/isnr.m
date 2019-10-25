function [ratio]=isnr(y, y_tilde, x_hat, A)

    ratio=20*log10(norm(y-y_tilde)/norm(A*x_hat-y));
    
end