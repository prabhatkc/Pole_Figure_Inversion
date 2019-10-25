clear all; clc;
addpath('../../../dependencies/')
addpath('../../../dataset/Cu-incomplete')
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if (isOctave ==0)
    addpath('../../dataset/matlab')
else
    addpath('../../dataset/octave')
end 

load("lsqr_cu_inputs_n_results.mat");
load("BigA-Cu-ncub10.txt");

m    = length(dataset.b);
n    = length(result.all_recon);
A    = sparse(BigA(:,1), BigA(:,2), BigA(:,3), m, n);
% 139.605s for 2000 ite
Xdim    = dataset.Xdim;
Ydim    = dataset.Ydim;
input_pf  = dataset.b;
recon_pf  = A*result.all_recon;
recon_odf = result.all_recon;

figure, suptitle({'reconstructed ODF for cu'}),
plotLayers(reshape(recon_odf(:, 1), Xdim)); set(gcf, 'Position', [100, 100, 800, 800]);
print('reconstructed_cu_odf.eps', '-depsc2');
close all;

disp("all the Polefigure inputs and outputs are saved in the folder....");
disp(pwd);
