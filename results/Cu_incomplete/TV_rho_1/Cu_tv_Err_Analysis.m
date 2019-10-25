clear all; clc;
addpath('../../../dependencies/')
addpath('../../../dataset/Cu-incomplete')
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if (isOctave ==0)
    addpath('../../dataset/matlab')
else
    addpath('../../dataset/octave')
end 

load("tv_Cu_inputs_n_results.mat");
BigA = load("BigA-Cu-ncub10.txt");

m    = length(dataset.b);
n    = length(result.recon);
A    = sparse(BigA(:,1), BigA(:,2), BigA(:,3), m, n);
% 139.605s for 2000 ite
Xdim    = dataset.Xdim;
Ydim    = dataset.Ydim;
input_pf  = dataset.b;
recon_pf  = A*result.recon;
recon_odf = result.recon;
snr       = result.snr';

figure(1); clf; 
h=plot(1:length(snr), snr, 'k');
set(h,'LineWidth',2);
max_x = length(snr);
max_y = max(snr)+0.5;
min_y = min(snr);
axis([ 0 max_x min_y max_y]); 
h=xlabel('Iteration');
h=ylabel('SNR [dB]');
set(gca,'Fontsize',17);
grid on;
set(gca,'GridLineStyle','--')
print -depsc2 tv_snr_cu.eps
close;

figure, suptitle({'reconstructed ODF for cu'}),
plotLayers(reshape(recon_odf(:, 1), Xdim)); set(gcf, 'Position', [100, 100, 800, 800]);
print('reconstructed_cu_tv_odf.eps', '-depsc2');
close all;

disp("all the Polefigure inputs and outputs are saved in the folder....");
disp(pwd);
