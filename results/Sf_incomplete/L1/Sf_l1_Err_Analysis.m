clear all; clc;
addpath('../../../dependencies/')
addpath('../../../dataset/')
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if (isOctave ==0)
    addpath('../../../dataset/matlab')
else
    addpath('../../../dataset/octave')
    pkg load image;
end 
load('l1_santaFe_inputs_n_results.mat');

m       = 2595;
n       = 9261;
BigA    = load('../../../dataset/matlab/BigA-SF-ncub10-unnormalized-Incomplete.txt');
A       = sparse(BigA(:,1), BigA(:,2), BigA(:,3), m, n); 
Xdim    = dataset.Xdim;
gt 	    = dataset.gt;
[m n]   = size(A);
alpha   = dataset.alpha;
bkg     = dataset.bkg;

nlevels = length(alpha);

E            = zeros(nlevels, 1); % ODF error for 4 different levels of Poisï¿½son's
recon_odf    = zeros(n, nlevels); % reconstucted ODF from PF for 5 error levels
psnr         = zeros(nlevels, 1) ; % maximum attained psnr of each noise levels;
snr          = zeros(nlevels, 1);

for i = 1:nlevels
    E(i)            = odf_err(gt, result.all_recon(:,i));
    psnr(i)         = pf_psnr(A*result.all_recon(:,i), dataset.b_sim);
    snr(i)          = pf_snr(dataset.all_pf(:, i), result.all_recon(:, i), A);
    recon_odf(:, i) = result.all_recon(:, i);
end

plotLayers(reshape(gt, Xdim));
print('gt_Sf_odf.eps', '-depsc2');  set(gcf, 'Position', [100, 100, 800, 750]);
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RMSE & ODF PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:nlevels
    pf_string.title     = sprintf(['L1 resolved ODF from incomplete PF \n corrupted with poisson noise a = ', num2str(alpha(i)), ' & b = ', num2str(bkg)]);
    figure, suptitle(pf_string.title);
    plotLayers(reshape(recon_odf(:, i), Xdim)); set(gcf, 'Position', [100, 100, 800, 750]);
    out_filename = sprintf(['l1_recovered_odf_a_',num2str(alpha(i)),'_b_',num2str(bkg), '.eps']);
    print(out_filename, '-depsc2');
    close;
end


pf_tick1 = sprintf(['(',num2str(bkg),', ' num2str(alpha(1)),')']);
pf_tick2 = sprintf(['(',num2str(bkg),', ' num2str(alpha(2)),')']);
pf_tick3 = sprintf(['(',num2str(bkg),', ' num2str(alpha(3)),')']);
pf_tick4 = sprintf(['(',num2str(bkg),', ' num2str(alpha(4)),')']);
pf_tick={pf_tick1, pf_tick2, pf_tick3, pf_tick4};

figure, h = plot(1:nlevels, E) ;
set(h, 'marker', 'o', 'Color', 'k');
set(h, 'LineWidth', 2, 'markersize', 5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
axis([ 1 (nlevels+1) min(E-0.03) max(E+0.03)]);
title("RMSE plot analysis")
ylabel('ODF Error values');
xlabel('Reconstructed ODF with different Poisson noise');

for i = 1:nlevels
    txt = sprintf(['\n', num2str(E(i))]);
    text(i+.05, E(i), txt);
end
xticks(1:length(alpha));
xticklabels(pf_tick);
grid on;
set(gca,'GridLineStyle','--');
print -depsc2 odf_err_l1_sf_incomplete.eps
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNR PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure, h = plot(1:nlevels, snr);
set(h, 'marker', 'o');
set(h, 'LineWidth', 2, 'markersize', 5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');

axis([ 0 (nlevels+1) min(snr-0.5) max(snr+0.5)]);
title("SNR plot analysis")
ylabel('SNR [dB]');
xlabel('Reconstructed PF corrupted with different Poisson noise');

xticks(1:length(alpha));
xticklabels(pf_tick);

for i = 1:nlevels
    txt = sprintf(['\n', num2str(snr(i))]);
    text(i+.05, snr(i), txt);
end
grid on;
set(gca,'GridLineStyle','--');

print -depsc2 pf_snr_l1_sf_incomplete.eps
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSNR PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure, h = plot(1:nlevels, psnr);
set(h, 'marker', 'o', 'Color', 'k');
set(h, 'LineWidth', 2, 'markersize', 5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');

axis([ 0 (nlevels+1) min(psnr-0.5) max(psnr+0.5)]);
title("PSNR plot analysis")
ylabel('PSNR [dB]');
xlabel('Reconstructed PF corrupted with different Poisson noise');
xticks(1:length(alpha));
xticklabels(pf_tick);

for i = 1:nlevels
    txt = sprintf(['\n', num2str(psnr(i))]);
    text(i+.1, psnr(i), txt);
end
grid on;
set(gca,'GridLineStyle','--');

print -depsc2 pf_psnr_l1_sf_incomplete.eps
close;

disp("all the Polefigure inputs and outputs are saved in the folder....");
disp(pwd);
