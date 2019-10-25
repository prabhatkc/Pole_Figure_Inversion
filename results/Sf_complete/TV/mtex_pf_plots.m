%% Prepare pole figure data from the output of the optimization
clc
clear all
mtex_path = '../../../../../mtex-5.1.1/';
addpath(mtex_path);
startup;

mkdir 'temp_txt';
mkdir 'pf_err_figs'
temp_dir = 'temp_txt';
imgdir  = 'pf_err_figs';

% read angles
angles = load('angles.txt');
%pf_opt = load('pfs_1.6_0.8_0.4_0.1.txt');
load('santaFe_inputs_n_results.mat')
pf_opt = dataset.all_pf;
hkl = {'111', '200', '220'};
% write out files to be read with MTEX
for i = 1:4
    for j = 1:3
        fname = strcat(temp_dir, '/test_(',hkl{j},')_nl',num2str(i),'.txt');
        fid = fopen(fname,'w');
        for k = 1:1297
            fprintf(fid,'%f\t%f\t%f\n',angles(k,1),angles(k,2),pf_opt((j-1)*1297+k,i));
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

% second set of pole figures for noise level 2
fnamenl2 = {...
  [temp_dir '/test_(111)_nl2.txt'],...
  [temp_dir '/test_(200)_nl2.txt'],...
  [temp_dir '/test_(220)_nl2.txt'],...
  };

% and so on
fnamenl3 = {...
  [temp_dir '/test_(111)_nl3.txt'],...
  [temp_dir '/test_(200)_nl3.txt'],...
  [temp_dir '/test_(220)_nl3.txt'],...
  };

% final set
fnamenl4 = {...
  [temp_dir '/test_(111)_nl4.txt'],...
  [temp_dir '/test_(200)_nl4.txt'],...
  [temp_dir '/test_(220)_nl4.txt'],...
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

pfnl2 = loadPoleFigure(fnamenl2,h,CS,SS,'interface','generic',...
  'ColumnNames', { 'Polar Angle' 'Azimuth Angle' 'Intensity'});

pfnl3 = loadPoleFigure(fnamenl3,h,CS,SS,'interface','generic',...
  'ColumnNames', { 'Polar Angle' 'Azimuth Angle' 'Intensity'});

pfnl4 = loadPoleFigure(fnamenl4,h,CS,SS,'interface','generic',...
  'ColumnNames', { 'Polar Angle' 'Azimuth Angle' 'Intensity'});

% finally plot everything out

figure;
plot(pfnl1,'contourf')
CLim(gcm,'equal')
mtexColorbar('LineWidth',2.5,'FontWeight', 'bold');
saveFigure([imgdir, '/noiselevel1.png']);

figure;
plot(pfnl2,'contourf')
CLim(gcm,'equal')
mtexColorbar('LineWidth',2.5,'FontWeight', 'bold');
saveFigure([imgdir, '/noiselevel2.png'])

figure;
plot(pfnl3,'contourf')
CLim(gcm,'equal')
mtexColorbar('LineWidth',2.5,'FontWeight', 'bold');
saveFigure([imgdir, '/noiselevel3.png'])

figure;
plot(pfnl4,'contourf')
CLim(gcm,'equal')
mtexColorbar('LineWidth',2.5,'FontWeight', 'bold');
saveFigure([imgdir, '/noiselevel4.png'])

A = dir(temp_dir);
for k = 3:length(A)
    delete([ temp_dir '/' A(k).name])
end
close all;
rmdir 'temp_txt';