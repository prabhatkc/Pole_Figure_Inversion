%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   main_Sf makes use of ADMM method to invert simulated polefigures (PFs)
%   obtained from the popular SantaFe dataset. Furthermore, these PFs
%   are corrupted with different levels of poisson's noise 
%	
%	References:
% 	------------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
addpath('../dependencies/');
addpath('../dependencies/lsqrSOL/');
addpath('../dataset/');
addpath('../dataset/SF-complete/');
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if (isOctave ==0)
    addpath('../dataset/matlab')
else
    addpath('../dataset/octave')
	pkg load statistics;
	pkg load image;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Dimensions of the ODF and its PoleFigures 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m    = 3891;
n    = 9261;

Xdim = int8(n^(1/3)); % dim of ODF
Ydim = int8(m/3); % dim of Polefigure (PF)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Upload datasets 
%  * BigA : is determined from the software EMsoft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BigA  = load('BigA-SF-ncub10-unnormalized.txt');
A     = sparse(BigA(:,1), BigA(:,2), BigA(:,3), m, n);                                
gt    = load('SF-GT.txt'); % upload Groundtruth ODF
b_sim = A*gt; % generates simulated PF of size (m, 1)
clear BigA;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input dataset/parameters to the solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('SantaFe_input_noise_pf_complete.mat');
dataset.all_pf 	    = all_pf;
dataset.Xdim        = [Xdim Xdim Xdim];
dataset.Ydim        = Ydim;
dataset.A           = A;
dataset.b_sim       = b_sim;
dataset.gt          = gt;
dataset.prompt_disp	= 'T'; %T to display; F to hide (for TV options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  parameters to induce possion's noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataset.alpha  = alpha; 
dataset.bkg    = bkg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimization methods (admm_opts)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% primary_methods | sub_optimal_method 
% ----------------|-------------------------
%	'TV-3D'		  | 'lu' or 'cgs' or 'inv'*
%	'l1'		  | 'non-weighted' or weighted
%	'lsqr'        | -
% 
%  *('inv' may be infeasible for very large 
%	A matrix) 
%
% para_est: T or F
%         : if T then lambda search routine is
%           switched on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

admm_opts.primary_method      ='TV-3D'; 
admm_opts.sub_optimal_method  ='lu'; 
admm_opts.para_est            = 'F';
invert_3_pf(dataset, admm_opts);







