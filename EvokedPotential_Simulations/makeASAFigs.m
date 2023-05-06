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

%% Extract Normal and Impaired EFRs, consolidate data consistently
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

%de-mean
n_mod_efr = n_mod_efr-mean(n_mod_efr,1);
i_mod_efr = i_mod_efr-mean(i_mod_efr,2);

%physiology
% Concat data in a similar manner
% Col 1 = SAM, Col 2 = SQ25, Col 3 = Harmonic rank 5, Col 4 = Rank 13
% data in slightly diff formats, but have fs = 8kHz.
n_sam_p = mean(phys_data.data_out.sam_all_n,1);
n_sam_p = n_sam_p((fs*t_win(1)+1):fs*t_win(2))';
n_sq25_p = mean(phys_data.data_out.sq25_all_n,1);
n_sq25_p = n_sq25_p((fs*t_win(1)+1):fs*t_win(2))';
i_sam_p = mean(phys_data.data_out.sam_all_i,1);
i_sam_p = i_sam_p((fs*t_win(1)+1):fs*t_win(2))';
i_sq25_p = mean(phys_data.data_out.sq25_all_i,1);
i_sq25_p = i_sq25_p((fs*t_win(1)+1):fs*t_win(2))';

n_pitch_1_p = phys_pitch_data.pool_efr_T_NH(:,6);
n_pitch_1_p = n_pitch_1_p((fs*t_win(1)+1):fs*t_win(2));
n_pitch_2_p = phys_pitch_data.pool_efr_T_NH(:,1);
n_pitch_2_p = n_pitch_2_p((fs*t_win(1)+1):fs*t_win(2));
i_pitch_1_p = phys_pitch_data.pool_efr_T_CA(:,6);
i_pitch_1_p = i_pitch_1_p((fs*t_win(1)+1):fs*t_win(2));
i_pitch_2_p = phys_pitch_data.pool_efr_T_CA(:,1);
i_pitch_2_p = i_pitch_2_p((fs*t_win(1)+1):fs*t_win(2));

%pool, normalize to 1, de-mean
n_phys_efr = [n_sam_p,n_sq25_p,n_pitch_1_p,n_pitch_2_p];
phys_norms = max(n_phys_efr,[],1);
n_phys_efr = n_phys_efr./phys_norms;
n_phys_efr = n_phys_efr-mean(n_phys_efr,1);

i_phys_efr = [i_sam_p,i_sq25_p,i_pitch_1_p,i_pitch_2_p];
i_phys_efr = i_phys_efr./phys_norms;
i_phys_efr = i_phys_efr-mean(i_phys_efr,1);

%% Spectral Analysis
%Consider using a slight delay or taper to avoid onset 

%Model FFT
nfft = 2^nextpow2(size(n_mod_efr,1));
f = linspace(0,fs/2,nfft/2);

n_mod_fft = abs(fft(n_mod_efr,nfft));
n_mod_fft = n_mod_fft(1:end/2,:)*2;
i_mod_fft = abs(fft(i_mod_efr,nfft));
i_mod_fft = i_mod_fft(1:end/2,:)*2;

%Physiology FFT
n_phys_fft = abs(fft(n_phys_efr,nfft));
n_phys_fft = n_phys_fft(1:end/2,:)*2;
i_phys_fft = abs(fft(i_phys_efr,nfft));
i_phys_fft = i_phys_fft(1:end/2,:)*2;

%% Plotting Figures

%Plot Parameters:

%colors
%font size
%linewidth
l_wdth = 1.5;

%fig position/size
fig_dims = [1540 397 1393 564];


%2 panel figure with model and physiological time waveforms side by side.
%Choose a single cihc

%reformatting model vectors to match physiology
%Col 1 = SAM, Col 2 = Sq25, Col 3 = Rank 5, Col 4 = Rank 13
sim_indexes = [10,12,1,6];

ihc = 4;
ihc_val = mod_data.ihc_grades(ihc);

n_to_plot_mod = n_mod_efr(:,sim_indexes);
i_to_plot_mod = squeeze(i_mod_efr(ihc,:,sim_indexes));

%get the cihc value specified, can loop thru this and generate fig if
%appropriate

buff = 0.2;
phys_scale_fact = 0.75; %arbitrary at this point, for vis purposes.

buff = buff*1:size(n_to_plot_mod,2);

t_waveform_model_phys = figure();
subplot(1,2,1);
title('Model')
hold on
plot(n_to_plot_mod+buff,'k','linewidth',l_wdth)
plot(i_to_plot_mod+buff,'r','linewidth',l_wdth)
hold off
ylim([-1,4.5]);
%FIX
legend('Normal',['C_{ihc} = ',num2str(ihc_val)]);

subplot(1,2,2)
title('Physiology')
hold on
plot(phys_scale_fact*n_phys_efr+buff,'k','linewidth',l_wdth)
plot(phys_scale_fact*i_phys_efr+buff,'r','linewidth',l_wdth)
hold off
legend('Normal',['Carboplatin']);
ylim([-1,4.5]);
set(gcf,'Position',fig_dims)
%% Transduction Figures
%% Save
