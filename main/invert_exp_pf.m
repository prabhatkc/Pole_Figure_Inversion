%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%   invert_exp_pf estimates ODF from the input polefigures in rastered form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Inputs:
%   --------
%   dataset
%       * Xdim: length of edges of 3D-ODF in x, y & z direction. The dimension 
%       in each direction is 2*ncub + 1
%       * Ydim: length of edges of 2D - PF in x & y direction.Concretely,the 
%       dimensions are 2*nlam +1
%       * A: Forward projection matrix when applied to ODF yields
%       Polefigure. It has the dimension [m n]
%       * b: Is the Polefigure obtained from experiments in rastered form
%
%   admm_opts  
%       * primary_method: string of either 'l1' or 'TV-3D'. l1 will
%       solve a lasso cost function to invert the PF where as TV-3D
%       will solve a total variation based cost function. l1 method is
%       used as a reference to get a sense of boundry parameters and
%       other admm parameters like rho, alpha & lambda. Thus,
%       determined parameters are used in TV-3D method to deduce a
%       more robust ODF solution.
%       * sub_optimal_method: sub_optimal_method is used to choose a 
%       specific route while solving the x-step of the ADMM method.
%       choose 'weighted' or 'non-weighted' if 'l1' is the primary method;
%       choose 'inv' or 'lu' or 'cgs' if 'TV-3D' is the primary method. 
%           * 'inv' : Proceeds by determining direct inverse in the x-step
%           * 'lu'  : Proceeds by determing Cholesky decomposition in the x-step (default)
%           * 'cgs' : Proceeds by matlab's inbuilt conjugate gradient descent.
%
% Optimization paramters:
% ----------------------- 
%       * follow guide on choosing paramters will be helpful 
%       * decrease lambda for slower convergence. It also aids in smooth minimization of cost
%       * lower rho or alpha for faster convergence 
%       * increase rho or alpha for slower but more accurate reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = invert_exp_pf(dataset, admm_opts)
    
    if nargin < 2
        error('not enough inputs, try again \n');
    end

    if ~isfield(admm_opts,'primary_method')
        primary_method = 'TV-3D';
    else
        primary_method = admm_opts.primary_method;
    end
    [m n] = size(dataset.A);
    
    switch primary_method
        
        case 'lsqr'
            damp = 0; %zero damp solves specifically argmin ||Ax -b||_2
            show = 1; %zero will switch of displays
            maxite = 5;
            
            t_start  = tic;
            fprintf("\n\t ------------------------------------------------------------------- \n");
            fprintf("\t  Working on the dataset with experimental Cu pole figure  \n");
            fprintf("\t ------------------------------------------------------------------- \n");
                
            [x, ~] = lsqrSOL(m, n, dataset.A, dataset.b, damp, [], [], [], maxite, show);
            ind    = find(x<0);
            x(ind) = 0;
                
            result.all_recon = x;
            fprintf("SNR = %f [dB]\n", pf_snr(dataset.b, x, dataset.A)); %11.3019
            toc(t_start)
            
            %figure, plotLayers(reshape(x, [21 21 21]))
            dataset = rmfield(dataset, {'A'});
            disp ('saving all the input datasets and the resulting inverse solutions from the LSQR method')
            mkdir('../results/Cu_incomplete/lsqr');
            save ('../results/Cu_incomplete/lsqr/lsqr_cu_inputs_n_results.mat', 'result', 'dataset');
        
        case 'l1'
            convergence_para.MAX_ITER = 100; %1200; %
            convergence_para.ABSTOL   = 1e-4;
            convergence_para.RELTOL   = 1e-2;
            convergence_para.mask     = 'F';
            convergence_para.disp     = 'True';

            admm_opts.rho    = 1e1; %
            admm_opts.alpha  = 1.8;
            admm_opts.lambda = 0.001312; %lambda_estimation(dataset, admm_opts, convergence_para); %1e-6*norm(dataset.A'*dataset.b, inf); 

            [x, history]            = admm_lasso(dataset, admm_opts, convergence_para);
            history.plot            = 1; % 1 -> yes, 0 -> no
            im.plot                 = 1; % 1 -> yes, 0 -> no
            im.recon                = reshape(x, dataset.Xdim);
            im.title_r              = 'Reconstructed ODF using L1 regularization';
            result.recon            = x;
            result.snr              = history.snr;
            plots_n_figures_exp(history, im);    
            %close all;
            dataset = rmfield(dataset, {'A'});
            disp ('saving all the resulting inverse solutions')
            mkdir('../results/Cu_incomplete/l1_rho1e1_ite100');
            save ('../results/Cu_incomplete/l1_rho1e1_ite100/l1_Cu_inputs_n_results.mat', 'result', 'dataset');
            
        case 'TV-3D'
            % convergence paramters
            convergence_para.MAX_ITER = 100; %1200; %
            convergence_para.ABSTOL   = 1e-4;
            convergence_para.RELTOL   = 1e-2;
            convergence_para.mask     = 'F';
            convergence_para.disp     = 'True';
                       
            admm_opts.rho       	= 1;
            admm_opts.alpha     	= 1.8;
            admm_opts.lambda    	= 0.000131; %lambda_estimation(dataset, admm_opts, convergence_para); %1e-5*norm(dataset.A'*dataset.b, inf);
            [x, history]            = admm_tv_fb_3D_exp(dataset, admm_opts, convergence_para);
            history.plot            = 1; % 1 -> yes, 0 -> no
            im.plot                 = 1; % 1 -> yes, 0 -> no
            im.recon                = reshape(x, dataset.Xdim);
            im.title_r              = 'Reconstructed ODF using TV-3D regularization';
            result.recon            = x;
            result.snr              = history.snr;
            plots_n_figures_exp(history, im);    
            %close all;
            dataset = rmfield(dataset, {'A'});
            disp ('saving all the resulting inverse solutions')
            mkdir('../results/Cu_incomplete/TV_rho1_ite100');
            save ('../results/Cu_incomplete/TV_rho1_ite100/tv_Cu_inputs_n_results.mat', 'result', 'dataset');
        end
end
                                           
function [lambda_old] = lambda_estimation (dataset, admm_opts, convergence_para)
    fprintf("\n .. Performing lambda estimation .. \n");
    convergence_para.disp = 'False';
    bp                    = dataset.A'*dataset.b;
    lambda_old            =  norm(bp, inf);
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  lambda bound estimation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda_list = [lambda_old];
    snr_list    = [-10]; %dummy placeholder 
    
    fprintf("\n .. (Part 1) 0.1 step size based bound estimation of lambda .. \n");
    while 1
        lambda_new = 0.1*lambda_old;
        admm_opts.lambda = lambda_new;
        if strcmpi(admm_opts.primary_method, 'l1')
           [x, history] = admm_lasso(dataset, admm_opts, convergence_para);
        else
           [x, history] = admm_tv_fb_3D_exp(dataset, admm_opts, convergence_para);
        end
        %[lx, ly] = lcurve_point(x, dataset.b, dataset.A); 
        snr      = pf_snr (dataset.b, x, dataset.A); 
        fprintf("(old, new) lambda = (%f,  %f)\n", lambda_old, lambda_new);
        fprintf("(old, new) SNR = (%f, %f)\n\n", snr_list(end), snr);
        
        lambda_list = [lambda_list lambda_new];
        snr_list  = [snr_list snr];
        
        %if (rmse_new > rmse_old && abs(snr_list(end) - snr_list(end-1))<=0.1)
        if (abs(snr_list(end) - snr_list(end-1))<=0.1 ||...
                (snr_list(end)-snr_list(end-1)<=-1) )
            break;
        end
        
        lambda_old   = lambda_new;
    end
    lambda_list = lambda_list(2:(end));   
    snr_list = snr_list(2:(end));
    
    fprintf("\n .. corresponding SNR values for lambda choices .. \n");
    fprintf('lambda:\t'); fprintf('%f\t', [lambda_list]); fprintf('\n');
    fprintf('SNR   :\t'); fprintf('%f\t', [snr_list]); fprintf('\n')
    para_est_plot(snr_list, lambda_list, 'bound');
    
  
end

function para_est_plot( snr_list, lambda_list, type)

    % RMSE plot
    figure;
    if strcmp(type, 'precise')
        set(gcf, 'Position', [110, 100, 1250, 500]); 
    else
        set(gcf, 'Position', [110, 100, 800, 500]); 
    end
   
    %SNR plot
    figure;
    if strcmp(type, 'precise')
        set(gcf, 'Position', [110, 100, 1250, 500]); 
    else 
        set(gcf, 'Position', [110, 100, 800, 500]); 
    end
    y = snr_list;
    x = lambda_list;
    h = plot(1:length(x), y, 'o');
    xlim([0 length(x)+1])
    set(gca,'xtick',[]);
    %set(h, 'marker', 'o');
    set(h, 'LineWidth', 2, 'markersize', 5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    ylabel("SNR [dB]");
    set(h, 'color', 'k')
    set(gca, 'Fontsize', 17);
    if strcmp(type, 'precise')
        set(gcf, 'Position', [110, 100, 1250, 500]); 
    end
    plot_prompt('F');

end

function [] = plot_prompt(disp)
    if strcmp(disp, 'F')
        for k =1:1
          close all;
          continue;
        end
    else 

        reply = input('Enter filename to save parameter est. plot ELSE enter [r] to continue: ', 's');
        for k = 1:1
            if strcmp(reply(1), 'r')
                close all;
                continue;
            else
                saveas(gca, fname);
                close;
            end 
        end
    end
end


