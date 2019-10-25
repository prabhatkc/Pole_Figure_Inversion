%% Prepare pole figure data from the output of the optimization
clc
clear all
mtex_path = '../../../../../mtex-5.1.1/';
addpath(mtex_path);
startup;

addpath('../../../dataset/Cu-incomplete')
mkdir 'temp_txt';
mkdir 'pf_tv_figs'
temp_dir = 'temp_txt';
imgdir  = 'pf_tv_figs';

% read angles
angles = load('angles.txt');
Na = length(angles);
load('tv_Cu_inputs_n_results.mat');
BigA= load("BigA-Cu-ncub10.txt");

m    = length(dataset.b);
n    = length(result.recon);
A    = sparse(BigA(:,1), BigA(:,2), BigA(:,3), m, n);

pf_opt = A*result.recon;
hkl = {'111', '200', '220'};
% write out files to be read with MTEX
for i = 1:1
    for j = 1:3
        fname = strcat(temp_dir, '/test_(',hkl{j},')_nl',num2str(i),'.txt');
        fid = fopen(fname,'w');
        for k = 1:Na
            fprintf(fid,'%f\t%f\t%f\n',angles(k,1),angles(k,2),pf_opt((j-1)*Na+k,i));
        end
        fclose(fid);
    end
end

%% Import Script for PoleFigure Data
%
% This script was automatically created by the import wizard. You should
% run the whoole script or parts of it in order to import your data. There
% is no problem in making any changes to this script.

%% Specify Crystal and Specimen Symmetries

% crystal symmetry
CS = crystalSymmetry('m-3m', [1 1 1]);

% specimen symmetry
SS = specimenSymmetry('1');

% plotting convention
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','outOfPlane');

%% Specify File Names

% path to files

% which files to be imported
fnamenl1 = {...
  [temp_dir '/test_(111)_nl1.txt'],...
  [temp_dir '/test_(200)_nl1.txt'],...
  [temp_dir '/test_(220)_nl1.txt'],...
  };

%% Specify Miller Indice

h = { ...
  Miller(1,1,1,CS),...
  Miller(2,0,0,CS),...
  Miller(2,2,0,CS),...
  };

%% Import the Data

% create a Pole Figure variable containing the data
pfnl1 = loadPoleFigure(fnamenl1,h,CS,SS,'interface','generic',...
  'ColumnNames', { 'Polar Angle' 'Azimuth Angle' 'Intensity'});

% finally plot everything out

figure;
plot(pfnl1,'contourf')
CLim(gcm,'equal')
mtexColorbar('LineWidth',2.5,'FontWeight', 'bold');
saveFigure([imgdir, '/exp_Cu_tv_recon_pf.png']);

A = dir(temp_dir);
for k = 3:length(A)
    delete([ temp_dir '/' A(k).name])
end
close all;
rmdir 'temp_txt';