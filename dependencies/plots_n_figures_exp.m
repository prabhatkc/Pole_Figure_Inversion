function plots_n_figures_exp(history, im)

if history.plot
    K = length(history.cost);  
    
    f=figure;
    plot(1:K, history.snr)
    title('SNR');
    ylabel('-20log{||f-f_m||_2/||f||_2}'); xlabel('iter (k)');
    
    h = figure;
    plot(1:K, history.cost, 'k', 'MarkerSize', 10, 'LineWidth', 2); 
    ylabel('f(x^k) + g(z^k)'); xlabel('iter (k)');

    g = figure;
    subplot(2,1,1);                                                                                                                    
    semilogy(1:K, max(1e-8, history.r_norm), 'k', ...
    1:K, history.eps_pri, 'k--',  'LineWidth', 2); 
    ylabel('||r||_2'); 

    subplot(2,1,2);                                                                                                                    
    semilogy(1:K, max(1e-8, history.s_norm), 'k', ...
    1:K, history.eps_dual, 'k--', 'LineWidth', 2);   
    ylabel('||s||_2'); xlabel('iter (k)'); 

end

if im.plot
    figure, suptitle(im.title_r), plotLayers(im.recon);    
end
end