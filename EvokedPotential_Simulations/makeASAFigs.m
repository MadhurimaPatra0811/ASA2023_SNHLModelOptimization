% Andrew Sivaprakasam
%Description: Plots Simulated and physiological data, in addition to other
% poster-related figures.

%TODO: 
% - Pass a labels vector from simulation code for easier plotting

%% Adding Paths
clear all, close all

addpath('Stimulus_Generation')
addpath('SimulatedData')
addpath('PhysiologicalData')
addpath('BEZ2018model/')
addpath('Functions')

fig_dir = 'Figures';
sim_name = 'sim_202305051452.mat';
phys_name = 'CA_Processed.mat';
phys_pitch_name = 'pitch_ca_tts_nh_efrs.mat';

if ~isfolder(fig_dir)
    mkdir(fig_dir);
end

cwd = pwd();
%% Load Model and Physiological Data:
%model
mod_data = load(sim_name);

%physiological
phys_data = load(phys_name);
phys_pitch_data = load(phys_pitch_name);

%% Extract Normal and Impaired EFRs
t_win = [0,.15]; 
fs_mod = mod_data.psth_fs;
fs = 8e3;

%model
n_mod_efr = mod_data.grand_envs_n;
n_mod_efr = resample(n_mod_efr,fs,fs_mod);
n_mod_efr = n_mod_efr((fs*t_win(1)+1):fs*t_win(2),:);
i_mod_efr = mod_data.grand_envs_i;
i_mod_efr = i_mod_efr(:,(fs*t_win(1)+1):fs*t_win(2),:);
model_normals = max(n_mod_efr,[],1);

%amplitude normalized to have maximum value of 1
n_mod_efr = n_mod_efr./model_normals;
i_mod_efr = bsxfun(@rdivide, i_mod_efr, permute(model_normals,[1,3,2]));

%physiology
% Concat data in a similar manner
% Col 1 = SAM, Col 2 = SQ25, Col 3 = Harmonic rank 5, Col 4 = Rank 13
% n_phys_efr = ;



%% Normalization by Norm Simulations and Baseline Values to Cross-Compare

%% Spectral Analysis

%% Plotting Figures

%% Transduction Figures
%% Save