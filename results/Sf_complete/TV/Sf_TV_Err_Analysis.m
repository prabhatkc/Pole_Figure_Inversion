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
load('santaFe_inputs_n_results.mat');

m    = 3891;
n    = 9261;
BigA  = load('../../../dataset/matlab/BigA-SF-ncub10-unnormalized.txt');
A     = sparse(BigA(:,1), BigA(:,2), BigA(:,3), m, n); 
Xdim    = dataset.Xdim;
gt 	    = dataset.gt;
[m n]   = size(A);
alpha   = dataset.alpha;
bkg     = dataset.bkg;
nlevels = length(alpha);

E            = zeros(nlevels, 1); % ODF error for 4 different levels of Poisï¿½son's
recon_odf    = zeros(n, nlevels); % reconstucted ODF from PF for 5 error levels
psnr     = zeros(nlevels, 1) ; % maximum attained psnr of each noise levels;
snr      = zeros(nlevels, 1);
for i = 1:nlevels
    E(i)            = odf_err(gt, result.all_recon(:,i));
    string{i}       = ['alpha=', num2str(dataset.alpha(i)),', bkg=', num2str(dataset.bkg)];
    ind             = find(result.all_psnr(:,i)==0); result.all_psnr(ind,i)=NaN;
    ind             = find(result.all_snr(:,i)==0); result.all_snr(ind,i)=NaN;
    ind             = find(~isnan(result.all_psnr(:, i))); psnr(i) = (result.all_psnr(ind(end), i));
    ind             = find(~isnan(result.all_snr(:, i))); snr(i) = (result.all_snr(ind(end), i));
    recon_odf(:, i) = result.all_recon(:, i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psnr plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iteration plot
%-----------------
max_itr = length(result.all_psnr(:,1));
figure(1); clf; 
h=plot(1:max_itr, result.all_psnr(:, 1), 'k');
set(h,'LineWidth',2);
hold on; 
h=plot(1:max_itr, result.all_psnr(:, 2),'y');
set(h,'LineWidth',2);
h=plot(1:max_itr, result.all_psnr(:, 3),'r');
set(h,'LineWidth', 2);
h=plot(1:max_itr, result.all_psnr(:, 4),'g');
set(h,'LineWidth',2);

max_x = max_itr;
max_y = max(max(result.all_psnr))+2;
min_y = min(min(result.all_psnr))-2;
axis([ 0 max_x min_y max_y]); 
h=xlabel('Iteration (k)');
set(h,'Fontsize',14);
h=ylabel('PSNR (x^{(k)})');
set(h,'Fontsize',12);
legend({string{1}, string{2}, string{3}, string{4}},'Location', 'southeast');
set(gca,'Fontsize',13);
print -depsc2 ite_psnr_sf.eps
close;

% final psnr plot
%-----------------
figure, h = plot(1:nlevels, psnr);
set(h, 'marker', 'o', 'Color', 'k');
set(h, 'LineWidth', 2, 'markersize', 5, 'MarkerEdgeColor', 'k');
axis([ 0 (nlevels+1) min(psnr-0.5) max(psnr+0.5)]);
title("PSNR plot analysis")
ylabel('PSNR [dB]');
xlabel('Reconstructed PF corrupted with different Poisson noise');
pf_tick1 = sprintf(['(',num2str(bkg),', ' num2str(alpha(1)),')']);
pf_tick2 = sprintf(['(',num2str(bkg),', ' num2str(alpha(2)),')']);
pf_tick3 = sprintf(['(',num2str(bkg),', ' num2str(alpha(3)),')']);
pf_tick4 = sprintf(['(',num2str(bkg),', ' num2str(alpha(4)),')']);
pf_tick={pf_tick1, pf_tick2, pf_tick3, pf_tick4};

xticks(1:length(alpha));
xticklabels(pf_tick);

for i = 1:nlevels
    txt = sprintf(['\n', num2str(psnr(i))]);
    text(i+.1, psnr(i), txt);
end
grid on;
set(gca,'GridLineStyle','--');

print -depsc2 pf_psnr_TV_sf_complete.eps
close;

%%%%%%%%%%
% snr plot 
%%%%%%%%%%%
max_itr = length(result.all_snr(:,1));

figure(2); clf; 
h=plot(1:max_itr, result.all_snr(:, 1), 'k');
set(h,'LineWidth',2);
hold on; 
h=plot(1:max_itr, result.all_snr(:, 2),'y');
set(h,'LineWidth',2);
h=plot(1:max_itr, result.all_snr(:, 3),'r');
set(h,'LineWidth', 2);
h=plot(1:max_itr, result.all_snr(:, 4),'g');
set(h,'LineWidth',2);

max_x = max_itr;
max_y = max(max(result.all_snr))+2;
min_y = min(min(result.all_snr))-2;
axis([ 0 max_x min_y max_y]); 
h=xlabel('Iteration (k)');
set(h,'Fontsize',14);
h=ylabel('SNR (x^{(k)})');
set(h,'Fontsize',12);
legend({string{1}, string{2}, string{3}, string{4}},'Location', 'southeast');
set(gca,'Fontsize',13);
print -depsc2 ite_snr_sf.eps
close;
%----------------------------
% final overall snr plot
%----------------------------
figure, h = plot(1:nlevels, snr);
set(h, 'marker', 's', 'Color', 'k');
set(h, 'LineWidth', 2, 'markersize', 5, 'MarkerEdgeColor', 'k','MarkerFaceColor', 'k');
axis([ 0 (nlevels+1) min(snr-0.5) max(snr+0.5)]);
title("SNR plot analysis")
ylabel('SNR [dB]');
xlabel('Reconstructed PF corrupted with different Poisson noise');
pf_tick1 = sprintf(['(',num2str(bkg),', ' num2str(alpha(1)),')']);
pf_tick2 = sprintf(['(',num2str(bkg),', ' num2str(alpha(2)),')']);
pf_tick3 = sprintf(['(',num2str(bkg),', ' num2str(alpha(3)),')']);
pf_tick4 = sprintf(['(',num2str(bkg),', ' num2str(alpha(4)),')']);
pf_tick={pf_tick1, pf_tick2, pf_tick3, pf_tick4};

xticks(1:length(alpha));
xticklabels(pf_tick);

for i = 1:nlevels
    txt = sprintf([ num2str(snr(i))]);
    text(i+.1, snr(i), txt);
end
grid on;
set(gca,'GridLineStyle','--');
print -depsc2 pf_snr_tv_sf_complete.eps
close;


%%%%%%%%%%
% ODF plot 
%%%%%%%%%%%
plotLayers(reshape(gt, Xdim));
print('gt_Sf_odf.eps', '-depsc2');  set(gcf, 'Position', [100, 100, 800, 750]);
close;

for i = 1:nlevels
  
    pf_string.title     = sprintf(['TV resolved ODF, initially, corrupted with poisson noise a = ', num2str(alpha(i)), ' & b = ', num2str(bkg)]);
    figure, suptitle(pf_string.title);
    plotLayers(reshape(recon_odf(:, i), Xdim)); set(gcf, 'Position', [100, 100, 800, 750]);
    out_filename = sprintf(['tv_recovered_odf_a_',num2str(alpha(i)),'_b_',num2str(bkg), '.eps']);
    print(out_filename, '-depsc2');
    close;
end

figure, h2 = plot(1:nlevels, E, 'LineWidth', 1.5);
set(h2, 'marker', 's', 'Color', 'k');
set(h2, 'LineWidth', 2, 'markersize', 5, 'MarkerEdgeColor', 'k','MarkerFaceColor', 'k');
axis([ 1 (nlevels+1) min(E-0.03) max(E+0.03)]);
title("RMSE plot analysis")
ylabel('ODF Error values');
xlabel('Reconstructed ODF with different Poisson noise');
for i = 1:nlevels
    txt = sprintf([num2str(E(i))]);
    text(i+.05, E(i), txt);
end
xticks(1:length(alpha));
xticklabels(pf_tick);
print -depsc2 odf_err_sf.eps
close;

clf;
hold on;
x = 1:1:nlevels;
y1 = psnr;
y2 = E;
[hax, h1, h2] = plotyy(x, y1, x, y2);
set(hax,'xtick',[1 2 3 4])
pf_tick1 = sprintf(['(',num2str(bkg),', ' num2str(alpha(1)),')']);
pf_tick2 = sprintf(['(',num2str(bkg),', ' num2str(alpha(2)),')']);
pf_tick3 = sprintf(['(',num2str(bkg),', ' num2str(alpha(3)),')']);
pf_tick4 = sprintf(['(',num2str(bkg),', ' num2str(alpha(4)),')']);
xticklabels(hax, {pf_tick1, pf_tick2, pf_tick3, pf_tick4});
set([h1, h2], 'linestyle', '-');
set(h1, 'marker', 'o');
set (h2, 'marker', 's');
set([h1, h2], 'LineWidth', 2, 'markersize', 9, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
ylabel(hax(1), "PSNR (dB)");
ylabel(hax(2), "RMSE");
xlabel('Noise parameters (\beta, \gamma)');
grid on;
set(gca,'GridLineStyle','--')
box on;
set([h1, h2], 'color', 'k')
set([hax(1), hax(2)], 'ycolor', 'k')
set([hax(1), hax(2)], 'FontSize', 17);
legend({'PSNR', 'RMSE'},'Location', 'south');
print -depsc2 Sf_rmse_v_psnr.eps
close;

disp("all the Polefigure inputs and outputs are saved in the folder....");
disp(pwd);
