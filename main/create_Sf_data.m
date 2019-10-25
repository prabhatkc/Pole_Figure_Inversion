%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   create_Sf_data forward projects the SantaFe ODF to the PoleFigures
% 	domain and adds poisson's noise 
%	
%	References:
% 	------------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

addpath('../dependencies/')
addpath('../dataset/SF-complete/')
addpath('../dataset/')
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

%%%%%%%%%%%%%%%%%%%%%%
%  Upload datasets 
%%%%%%%%%%%%%%%%%%%%%%
BigA  = load('BigA-SF-ncub10-unnormalized.txt'); % upload forward projection matrix
                                  % determined using software EmSoft
A     = sparse(BigA(:,1), BigA(:,2), BigA(:,3), m, n);                                
gt    = load('SF-GT.txt'); % upload Groundtruth ODF
b_sim = A*gt; % generates simulated PF of size (m, 1)

clear BigA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  parameters to induce possion's noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = [1.6, 0.4, 0.1, 0.025]; 
bkg   = 1;

all_pf  = zeros(m, length(alpha));
for i = 1:length(alpha)
    b_sim_noisy     = addPoissonNoise(alpha(i), bkg, b_sim);
    b               = m_normalize(min(b_sim), max(b_sim), b_sim_noisy);
    all_pf(:, i)    = b;
end

disp ('saving the simulated Polefigures corrupted with noise')
if (isOctave ==0)
    save ('../dataset/matlab/SantaFe_input_noise_pf_complete.mat', 'all_pf', 'alpha', 'bkg');
else
    save ('../dataset/octave/SantaFe_input_noise_pf_complete.mat', 'all_pf','alpha', 'bkg' );
end

