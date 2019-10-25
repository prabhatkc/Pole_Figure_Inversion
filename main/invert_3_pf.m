%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%   invert_3_pf estimates ODF from the 3 given polefigures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Inputs:
%   --------
%
%   dataset
%       * alpha, bkg: parameters to induce poissons. alpha = [0, 1] 
%       & bkg = 1 (default) 
%       * Xdim: length of edges of 3D-ODF in x, y & z direction. The dimension 
%       in each direction is 2*ncub + 1
%       * Ydim: length of edges of 2D - PF in x & y direction. Concretely, the 
%       dimensions are 2*nlam +1
%       * A: Forward projection matrix when applied to ODF yields
%       Polefigure. It has the dimension [m n]
%       * gt: is the simulated ODF that serves as the ground truth
%       * b_sim: is the pristine polefigure obtained as A*gt  without any
%       noise
%   admm_opts  
%       * primary_method: string of either 'l1' or 'TV-3D'. l1 will
%       solve a lasso cost function to invert the PF where as TV-3D
%       will solve a total variation based cost function. l1 method is
%       used as a reference to get a sense of boundary parameters and
%       other ADMM parameters like rho, alpha & lambda. Thus,
%       determined parameters are used in TV-3D method to deduce a
%       more robust ODF solution.
%       * sub_optimal_method: sub_optimal_method is used to choose a 
%       specific route while solving the x-step of the ADMM method.
%       choose 'weighted' or 'non-weighted' if 'l1' is the primary method;
%       choose 'inv' or 'lu' or 'cgs' if 'TV-3D' is the primary method. 
%           * 'inv' : proceeds by determining direct inverse in the x-step
%           * 'lu'  : proceeds by determining Cholesky decomposition in the x-step (default)
%           * 'cgs' : proceeds by matlab's inbuilt conjugate gradient descent.
%
% Optimization parameters:
% ----------------------- 
%        follow guide on choosing parameters will be helpful for a heuristic determination 
%       * decrease lambda for slower convergence. It also aids in smooth minimization of cost
%       * lower rho or alpha for faster convergence 
%       * increase rho or alpha for slower but more accurate reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = invert_3_pf(dataset, admm_opts)
    
    if nargin < 2
        error('not enough inputs, try again \n');
    end
    
    % default method
    if ~isfield(admm_opts,'primary_method')
        primary_method = 'TV-3D';
    else
        primary_method = admm_opts.primary_method;
    end
    
    [m n]            = size(dataset.A);
    alpha            = dataset.alpha;
    bkg              = dataset.bkg;
    
    result.all_recon = zeros(n, length(alpha));
    
    switch primary_method
        
        case 'lsqr'
            damp = 0; %zero damp solves specifically argmin ||Ax -b||_2
            show = 1; %zero will switch of displays
            maxite = 3;
            
            for i = 1:length(alpha) 
                t_start  = tic;
                fprintf("\n\t ------------------------------------------------------------------- \n");
                fprintf("\t  Working on the dataset with poisson noise (a = %.3f & b = %.2f)  \n", alpha(i), bkg);
                fprintf("\t ------------------------------------------------------------------- \n");
                
                [x, ~] = lsqrSOL(m, n, dataset.A, dataset.all_pf(:,i), damp, [], [], [], maxite, show);
                ind    = find(x<0);
                x(ind) = 0;
                
                result.all_recon(:, i) = x;
                fprintf("RMSE = %f\n", odf_err(dataset.gt, x));
                fprintf("PSNR = %f [dB]\n", pf_psnr(dataset.A*x, dataset.b_sim(:, 1)));
                fprintf("SNR = %f [dB]\n", pf_snr(dataset.all_pf(:, i), x, dataset.A));
                toc(t_start)
            end

            dataset = rmfield(dataset, {'A'});
            disp ('saving all the input datasets and the resulting inverse solutions from the LSQR method')
            mkdir('../results/Sf_complete/lsqr');
            save ('../results/Sf_complete/lsqr/lsqr_santaFe_inputs_n_results.mat', 'result', 'dataset');
        
        case 'l1'
            convergence_para.MAX_ITER = 3000;
            convergence_para.ABSTOL   = 1e-4;
            convergence_para.RELTOL   = 1e-2;
            convergence_para.mask     = "F";
            ncub                      = (dataset.Xdim(1)-1)/2;
            mask                      = cubicBoundryPrior( round(n^(1/3)), (ncub/2)+4);
            lambda_arr                = [28.078978, 100.798721, 107.142615, 173.452173];
           %time taken                   [42.611, 17.347, 15.27, 8.658 s]
            for i = 1:length(alpha)
                dataset.b               = dataset.all_pf(:, i);
                admm_opts.lambda        = lambda_arr(i); %1e-6*norm(dataset.A'*dataset.b, inf);
                admm_opts.rho           = 1e3; 
                admm_opts.alpha         = 1.8;
                [x, history]            = noise_level_optimization(dataset, alpha(i), bkg, admm_opts, convergence_para);
                result.all_recon(:, i)  = x; %.*mask(:);
                result.all_psnr(1:length(history.psnr), i) = history.psnr;
                result.all_snr(1:length(history.snr), i)   = history.snr;
            end                


            dataset = rmfield(dataset, {'b', 'A'});
            disp ('saving all the input datasets and the resulting inverse solutions')
            mkdir('../results/Sf_complete/L1');
            save ('../results/Sf_complete/L1/l1_santaFe_inputs_n_results.mat', 'result', 'dataset');
                                                
        case 'TV-3D'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  convergence paramters
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            convergence_para.MAX_ITER = 3000; 
            convergence_para.ABSTOL   = 1e-4;
            convergence_para.RELTOL   = 1e-2;
            convergence_para.mask     = "F";
           	ncub 			= (dataset.Xdim(1)-1)/2;
           	mask            = cubicBoundryPrior( round(n^(1/3)), (ncub/2)+4);
            
            %%%%%%%%%%%%%%%%%%%%
            %   a = 1.8
            %%%%%%%%%%%%%%%%%%%
            dataset.b               = dataset.all_pf(:, 1);
            admm_opts.lambda        = 10.747299; %0.9^50*0.0000001*norm(dataset.A'*dataset.b, inf); %; %45.02 %32.940500s %8.078978; %1e-8*norm(dataset.A'*dataset.b, inf); 8.078978; %79.872s
            admm_opts.rho           = 1e4; %1e4;
            admm_opts.alpha         = 1.8;
            [x, history]            = noise_level_optimization(dataset, alpha(1), bkg, admm_opts, convergence_para);
            result.all_recon(:, 1)  = x; %.*mask(:);
            result.all_psnr(1:length(history.psnr), 1) = history.psnr;
            result.all_snr(1:length(history.snr), 1)   = history.snr;
            result.all_rmse(1:length(history.odf_err),1)  = history.odf_err;
            %%%%%%%%%%%%%%%%%%%
            % a = 0.4
            %%%%%%%%%%%%%%%%%%%
%             dataset.b               = dataset.all_pf(:, 2);
%             admm_opts.lambda        = 718.849127; %17.26%15.093762%79.872125; %24.691s
%             admm_opts.rho           = 1e5; %1e6;
%             admm_opts.alpha         = 1.8;      
%             
%             [x history]             = noise_level_optimization(dataset, alpha(2), bkg, admm_opts, convergence_para);
%             result.all_recon(:, 2)  = x; %.*mask(:);
%             result.all_psnr(1:length(history.psnr), 2) = history.psnr; 
%             result.all_snr(1:length(history.snr), 2)   = history.snr;
%             %%%%%%%%%%%%%%%%%%%
%             %  a = 0.1
%             %%%%%%%%%%%%%%%%%%%
%             dataset.b           = dataset.all_pf(:, 3);
%             admm_opts.lambda    = 201.728533;%64.27 %53.69 %norm(dataset.A'*dataset.b, inf); %107.142615; %41.459s
%             admm_opts.rho       = 1e6; %1e6;
%             admm_opts.alpha     = 1.8;
%             
%             [x history]         = noise_level_optimization(dataset, alpha(3), bkg, admm_opts, convergence_para);
%             result.all_recon(:, 3)  = x; %.*mask(:);
%             result.all_psnr(1:length(history.psnr), 3)  = history.psnr; 
%             result.all_snr(1:length(history.snr), 3)   = history.snr;
%             
%             %%%%%%%%%%%%%%%%%%%
%             % a = 0.025
%             %%%%%%%%%%%%%%%%%%
%             dataset.b           = dataset.all_pf(:, 4);
%             
%             admm_opts.lambda    = 661.069557; %196.92s %173.452173; % 86.100s
%             admm_opts.rho       = 1e7; %1e6;
%             admm_opts.alpha     = 1.8;      
%            
%             [x history]             = noise_level_optimization(dataset, alpha(4), bkg, admm_opts, convergence_para);
%             result.all_recon(:, 4)	= x; %.*mask(:);
%             result.all_psnr(1:length(history.psnr), 4)	= history.psnr; 
%             result.all_snr(1:length(history.snr), 4)    = history.snr;
%             
%             dataset = rmfield(dataset, {'b', 'A'});
%             disp ('saving all the input datasets and the resulting inverse solutions');
%             mkdir('../results/Sf_complete/TV');
     		save ('../results/Sf_complete/TV/para_est/santaFe_inputs_n_results_lambda_best.mat', 'result', 'dataset', 'alpha', 'bkg');

        end
end
                                         

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ADMM Optimization calls
%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, history] = noise_level_optimization (dataset, noise_alpha, bkg, admm_opts, convergence_para)
    fprintf("\n\t ------------------------------------------------------------------- \n");
    fprintf("\t Working on the dataset with poisson noise (a = %.3f & b = %.2f)  \n", noise_alpha, bkg);
    fprintf("\t ------------------------------------------------------------------- \n");

    pf_string.title       = sprintf(['simulated PF corrupted with poisson noise \n \t a = ', num2str(noise_alpha), ' & b = ', num2str(bkg)]);
    pf_string.subtitle    = {'(2 0 0)', '(2 2 0)', '(1 1 1)'};
    
    if strcmpi(admm_opts.para_est, 'T')
        admm_opts.lambda = lambda_estimation(dataset, admm_opts, convergence_para);
    end
    
    fprintf("\nlambda determined: %f\n", admm_opts.lambda);
    convergence_para.disp = "True";
    if strcmpi(admm_opts.primary_method, 'l1')
        if strcmpi(admm_opts.sub_optimal_method(1), 'w')
            [x, history] = admm_w_lasso(dataset, admm_opts, convergence_para);
            im.title_r   = 'Reconstructed ODF using L1 regularization (lasso)';
        else
            [x, history] = admm_lasso(dataset, admm_opts, convergence_para);
            im.title_r   = 'Reconstructed ODF using L1 regularization (lasso)';
        end
    else
       [x, history] = admm_tv_fb_3D(dataset, admm_opts, convergence_para);
       im.title_r   = 'Reconstructed ODF using TV-3D regularization';
    end
    history.plot = 1; % 1 -> yes, 0 -> no
    im.plot      = 1; % 1 -> yes, 0 -> no
    im.recon     = reshape(x, dataset.Xdim);
    im.gt        = reshape(dataset.gt, dataset.Xdim);
    im.title_gt  = 'gt Layers';
    plots_n_figures(history, im);    
    prompt(dataset.prompt_disp);   
end     


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  lambda parameter search call
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
function [lambda_old] = lambda_estimation (dataset, admm_opts, convergence_para)
    fprintf("\n .. Performing lambda estimation .. \n");
    convergence_para.disp = 'False';
    bp                    = dataset.A'*dataset.b;
    rmse_old              =  odf_err( dataset.gt, bp);
    fprintf("rmse from comparision between GT and analytical method (i.e. A'*b): %f\n\n", rmse_old);
    lambda_old            =  0.01*norm(bp, inf);
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  lambda bound estimation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda_list = [lambda_old];
    rmse_list   = [rmse_old];
    snr_list    = [-10]; %dummy placeholder 
    lx_list     = [];  ly_list   = [];
    
    fprintf("\n .. (Part 1) 0.1 step size based bound estimation of lambda .. \n");
    while 1
        lambda_new = 0.9*lambda_old;
        admm_opts.lambda = lambda_new;
        if strcmpi(admm_opts.primary_method, 'l1')
            if strcmp(admm_opts.sub_optimal_method(1), 'w')
                [x, history] = admm_w_lasso(dataset, admm_opts, convergence_para);
            else
                [x, history] = admm_lasso(dataset, admm_opts, convergence_para);
            end
        else
           [x, history] = admm_tv_fb_3D(dataset, admm_opts, convergence_para);
        end
        %[lx, ly] = lcurve_point(x, dataset.b, dataset.A); 
        rmse_new = odf_err(dataset.gt, x);
        snr      = pf_snr (dataset.b, x, dataset.A); 
        fprintf("(old, new) lambda = (%f,  %f)\n", lambda_old, lambda_new);
        fprintf("(old, new) RMSE   = (%f,  %f), SNR = (%f, %f)\n\n", rmse_old, rmse_new, snr_list(end), snr);
        
        lambda_list = [lambda_list lambda_new];
        rmse_list = [rmse_list rmse_new];
        snr_list  = [snr_list snr];
        %lx_list = [lx_list lx]; ly_list = [ly_list ly];
        
        %if (rmse_new > rmse_old && abs(snr_list(end) - snr_list(end-1))<=0.1)
        if (abs(snr_list(end) - snr_list(end-1))<=0.1 ||...
                (snr_list(end)-snr_list(end-1)<=-1) )
            break;
        end
        
        lambda_old   = lambda_new;
        rmse_old     = rmse_new;
    end
    lambda_list = lambda_list(2:(end));   
    rmse_list = rmse_list(2:(end));
    snr_list = snr_list(2:(end));
    
    fprintf("\n .. corresponding RMSE and SNR values for lambda choices .. \n");
    fprintf('lambda:\t'); fprintf('%f\t', [lambda_list]); fprintf('\n');
    fprintf('RMSE  :\t'); fprintf('%f\t', [rmse_list]); fprintf('\n');
    fprintf('SNR   :\t'); fprintf('%f\t', [snr_list]); fprintf('\n')
    para_est_plot(rmse_list, snr_list, lambda_list, 'bound');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  lambda approx estimation
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     fprintf("\n .. (Part 2) precise lambda estimation .. \n");
%     lambda_old = lambda_list(end-1)+100;
%     admm_opts.lambda = lambda_old;
%     if strcmpi(admm_opts.primary_method, 'l1')
%         if strcmp(admm_opts.sub_optimal_method(1), 'w')
%             [x, history] = admm_w_lasso(dataset, admm_opts, convergence_para);
%         else
%             [x, history] = admm_lasso(dataset, admm_opts, convergence_para);
%         end
%     else
%        [x, history] = admm_tv_fb_3D(dataset, admm_opts, convergence_para);
%     end
%     
%     %[lx2, ly2] = lcurve_point(x, dataset.b, dataset.A);
%     rmse_old   = odf_err(dataset.gt, x);
%     snr_old    = pf_snr (dataset.b, x, dataset.A); 
% 
%     lambda_list2 = [lambda_old];
%     rmse_list2   = [rmse_old];
%     snr_list2    = [snr_old];
%     %lx_list2     = [lx2];  ly_list2   = [ly2];
% 
%     while 1
%         lambda_new = lambda_old-10;
%         admm_opts.lambda = lambda_new;
%         if strcmpi(admm_opts.primary_method, 'l1')
%             if strcmp(admm_opts.sub_optimal_method(1), 'w')
%                 [x, history] = admm_w_lasso(dataset, admm_opts, convergence_para);
%             else
%                 [x, history] = admm_lasso(dataset, admm_opts, convergence_para);
%             end
%         else
%            [x, history] = admm_tv_fb_3D(dataset, admm_opts, convergence_para);
%         end
%         %[lx2, ly2] = lcurve_point(x, dataset.b, dataset.A);
%         rmse_new = odf_err(dataset.gt, x);
%         snr      = pf_snr (dataset.b, x, dataset.A); 
% 
%         fprintf("(old, new) lambda = (%f,  %f)\n", lambda_old, lambda_new);
%         fprintf("(old, new) RMSE   = (%f,  %f), SNR = (%f, %f)\n\n", rmse_old, rmse_new, snr_list2(end), snr);
%         lambda_list2 = [lambda_list2 lambda_new];
%         rmse_list2   = [rmse_list2 rmse_new];
%         snr_list2    = [snr_list2 snr];
%         %lx_list2 = [lx_list2 lx2]; ly_list2 = [ly_list2 ly2];
%         
%         %if (rmse_new > rmse_old && abs(snr_list2(end) - snr_list2(end-1))<=0.1)
%         if ( abs(snr_list2(end) - snr_list2(end-1))<=0.1 || ...
%                 (snr_list2(end)-snr_list2(end-1)<=-0.5) || ...
%                 lambda_new<lambda_list(end))
%             break;
%         end
%         
%         lambda_old   = lambda_new;
%         rmse_old     = rmse_new;
% 
%     end
%     
%     fprintf("\n .. corresponding RMSE and SNR values for lambda choices .. \n");
%     fprintf('lambda:\t'); fprintf('%f\t', [lambda_list2]); fprintf('\n');
%     fprintf('RMSE  :\t'); fprintf('%f\t', [rmse_list2]); fprintf('\n');
%     fprintf('SNR   :\t'); fprintf('%f\t', [snr_list2]); fprintf('\n');
%     para_est_plot(rmse_list2, snr_list2, lambda_list2, 'precise');
end
   
function para_est_plot(rmse_list, snr_list, lambda_list, type)

    % RMSE plot
    figure;
    if strcmp(type, 'precise')
        set(gcf, 'Position', [110, 100, 1250, 500]); 
    else
        set(gcf, 'Position', [110, 100, 800, 500]); 
    end
    y = rmse_list;
    x = lambda_list;
    h = plot(1:length(x), y, 'o');
    xlim([0 length(x)+1])
    xticks(1:length(x));
    xticklabels({x'});

    if strcmp(type, 'bound')
        set(gca,'xticklabel',num2str(x','%.2e'));
    else 
        set(gca,'xticklabel',num2str(x','%.2f'));
    end

    %set(h, 'marker', 'o');
    set(h, 'LineWidth', 2, 'markersize', 5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    xlabel('\lambda')
    ylabel("RMSE");
    set(gca,'xaxisLocation','top')
    set(h, 'color', 'k')
    set(gca, 'Fontsize', 17);

    if strcmp(type, 'precise')
        set(gcf, 'Position', [110, 100, 1200, 500]); 
    end
    plot_prompt('F');

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

function [] = lcurve_plot(lx, ly)
figure, plot(lx, ly);
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