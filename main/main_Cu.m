
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   main_Sf makes use ADMM method to invert simulated polefigures (PFs)
%   obtained from the popular SantaFe dataset. Furthermore, these PFs
%   are corrupted with different levels of poisson's noise 
%	
%	References:
% 	------------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
addpath('../dependencies/')
addpath('../dataset/Cu-incomplete/')
addpath('../dependencies/lsqrSOL/');

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if (isOctave ==0)
    addpath('../dataset/matlab')
else
    pkg load statistics;
    addpath('../dataset/octave')
end

%%%%%%%%%%%%%%%%%%%%%%
%  Upload datasets 
%%%%%%%%%%%%%%%%%%%%%%
BigA  = load('BigA-Cu-ncub10.txt');
b_exp = load('b.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Dimensions of the ODF and its PoleFigures 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncub = 10;
nlam = 17;
m    = length(b_exp);
n    = (2*ncub+1)^3;

Xdim = [(2*ncub+1) (2*ncub+1) (2*ncub+1)]; % dim of ODF
Ydim = round(m/3); % dim of Polefigure (PF)
A     = sparse(BigA(:,1), BigA(:,2), BigA(:,3), m, n);
clear BigA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input dataset/parameters to the solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataset.Xdim        = Xdim;
dataset.Ydim        = Ydim;
dataset.A           = A;
dataset.b           = b_exp;
dataset.prompt_disp	= 'T'; %T to display; F to hide (for TV options)

admm_opts.primary_method      ='l1';  % l1 or TV-3D
admm_opts.sub_optimal_method  ='lu'; % 'non-weighted' or weighted for l1
                                     %  inv or cgs or lu for TV-3D
invert_exp_pf(dataset, admm_opts);
