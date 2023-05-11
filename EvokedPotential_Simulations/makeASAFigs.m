% Andrew Sivaprakasam
%Description: Plots Simulated and physiological data, in addition to other
% poster-related figures. This could all certainly be more efficient. 

%TODO: 
% - Pass a labels vector from simulation code for easier plotting

%% Adding Paths
clear all, close all

addpath('Stimulus_Generation')
addpath('SimulatedData')
addpath('PhysiologicalData')
addpath('Functions')

fig_dir = 'Figures/';
sim_name = 'sim_allcihc.mat';
sim_dir = [sim_name(1:end-4),'/'];
phys_name = 'CA_Processed.mat';
phys_pitch_name = 'pitch_ca_tts_nh_efrs.mat';

if ~isfolder([fig_dir,sim_dir])
    mkdir([fig_dir,sim_dir]);
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
t = t_win(1):1/fs:t_win(2)-1/fs;

%model
n_mod_efr = mod_data.grand_envs_n;
n_mod_efr = resample(n_mod_efr,fs,fs_mod);
n_mod_efr = n_mod_efr((fs*t_win(1)+1):fs*t_win(2),:);
i_mod_efr = mod_data.grand_envs_i;
i_mod_efr = resample(i_mod_efr,fs,fs_mod,'Dimension',2);
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
L = fs*(t_win(2)-t_win(1));

n_mod_fft = abs(fft(n_mod_efr,nfft)/L);
n_mod_fft = n_mod_fft(1:end/2,:)*2;
i_mod_fft = abs(fft(i_mod_efr,nfft,2)/L);
i_mod_fft = i_mod_fft(:,1:end/2,:)*2;

%Physiology FFT
n_phys_fft = abs(fft(n_phys_efr,nfft)/L);
n_phys_fft = n_phys_fft(1:end/2,:)*2;
i_phys_fft = abs(fft(i_phys_efr,nfft)/L);
i_phys_fft = i_phys_fft(1:end/2,:)*2;

%% R_DFT computation (come back to this another time)

% F0 = 103;
% 
% %physiology
% for s = 1:size(n_phys_fft,2)
%     [PKS,LOCS,n_floor] = getPeaks(f,n_phys_fft(:,s),F0,16); 
%     PKS_NH_p(:,s) = PKS-n_floor(LOCS); 
%     LOCS_NH_p(:,s) = LOCS;
%     nfloor_n(:,s) = n_floor;
% end 
% 
% Rats_NH_p = sum(PKS_NH_p(3:16,:))./sum(PKS_NH_p(1:2,:));
% 
% for s = 1:size(i_phys_fft,2)
%     [PKS,LOCS,n_floor] = getPeaks(f,i_phys_fft(:,s),F0,16); 
%     PKS_CA_p(:,s) = PKS-n_floor(LOCS); 
%     LOCS_CA_p(:,s) = LOCS;
%     nfloor_n(:,s) = n_floor;
% end 
% 
% Rats_CA_p = sum(PKS_CA_p(3:16,:))./sum(PKS_CA_p(1:2,:));

% for s = 1:size(fft_i,2)
%     [PKS,LOCS,n_floor] = getPeaks(f,fft_i(:,s),F0,16); 
%     PKS_IHC(:,s) = PKS-n_floor(LOCS); 
%     LOCS_IHC(:,s) = LOCS;
%     nfloor_i(:,s) = n_floor;
% end 
% 
% Rats_NH = sum(PKS_NH(3:16,:))./sum(PKS_NH(1:2,:));
% Rats_IHC = sum(PKS_IHC(3:16,:))./sum(PKS_IHC(1:2,:));
% 
% all_nh = sum(PKS_NH(1:16,:));
% all_ihc = sum(PKS_IHC(1:16,:));

%% Estimating 'best fit' cihc
%Since all signals are de-meaned, normalized by max of the norm, and time-locked 
%estimating a best fit using the mean absolute difference
%Some slight latency shifts...but model should account for it

%Computing for physiological data:
% diffs_phys  = n_phys_efr-i_phys_efr;

%try minimizing the best fit mag squared coherence between normal and
%impaired

%account for model data using 103 f0 for everything, and phys data being
%100 for sq and SAM, but 103 for everything else. Make sure to use correct
%indices!!


%Compute NH vs CA coherence (physiology)
for j = 1:size(n_phys_efr,2)
   
   %if SAM or Sq25, use multiples of 100
   if j == 1 || j==2
      freqs = 100*(1:16);
   else
      freqs = 103*(1:16);
   end

   [cohere, F] = mscohere(n_phys_efr(:,j),i_phys_efr(:,j),[],[],freqs,fs);
   phys_cohere(j) = mean(cohere);
end

phys_cohere = phys_cohere';

sim_indexes = [10,12,1,6];
n_mod_reform = n_mod_efr(:,sim_indexes);

%should be done more efficiently!!
%columns correspond to each cihc, row is for a stimulus
for c = 1:length(mod_data.ihc_grades)
    
    i_mod_reform = squeeze(i_mod_efr(c,:,sim_indexes));
    for j = 1:size(n_mod_reform,2)
       [cohere, F] = mscohere(n_mod_reform(:,j),i_mod_reform(:,j),[],[],freqs,fs);
       mod_cohere(j,c) = mean(cohere);
    end

end

err = abs(phys_cohere-mod_cohere);
[min_err, locs] = min(err,[],2);

best_fits = mod_data.ihc_grades(locs); 
%% Plotting Figures

%Plot Parameters:
%colors
blck = [0.25, 0.25, 0.25];
colors_ca = [0.8500, 0.3250, 0.0980, .85]; %last number is the face alpha/transparency

f_size = 13; %font size
l_wdth = 1.5; %linewidth
fig_dims = [1540 397 1393 564]; %fig position/size
ylims_t = [1,7]; %waveform ylim
ylims_s = [0.05,.55]; %spectral ylim

%2 panel figure with model and physiological time waveforms side by side.
%Choose a single cihc

%reformatting model vectors to match physiology
%Col 1 = SAM, Col 2 = Sq25, Col 3 = Rank 5, Col 4 = Rank 13
labs = ["SAM","Sq25","F_0 103 | Rank 3","F_0 103 | Rank 13"];
sim_indexes = [10,12,1,6];

for i = 1:length(mod_data.ihc_grades)
    ihc = i;
    ihc_val = mod_data.ihc_grades(i);
    
    n_to_plot_mod = n_mod_efr(:,sim_indexes);
    i_to_plot_mod = squeeze(i_mod_efr(ihc,:,sim_indexes));
    
    %get the cihc value specified, can loop thru this and generate fig if
    %appropriate
    
    buff = 1.5;
    phys_scale_fact = 0.65; %arbitrary at this point, for vis purposes.
    
    buff = buff*(1:size(n_to_plot_mod,2));
    
    t_waveform_model_phys = tiledlayout(1,2,'TileSpacing','compact');
    nexttile;
    title('Model Simulations')
    hold on
    plot(t,n_to_plot_mod+buff,'color',blck,'linewidth',l_wdth)
    plot(t,i_to_plot_mod+buff,'color',colors_ca,'linewidth',l_wdth)
    hold off
    yticks(buff*1.1);
    yticklabels([labs,'FontWeight','Bold']);
    ytickangle(45)
    ylim(ylims_t);
    %better way to do this?
    legend('Normal','','','','',['C_{ihc} = ',num2str(ihc_val)]);
    
    nexttile;
    title('In-Vivo')
    hold on
    plot(t,phys_scale_fact*n_phys_efr+buff,'color',blck,'linewidth',l_wdth)
    plot(t,phys_scale_fact*i_phys_efr+buff,'color',colors_ca,'linewidth',l_wdth)
    hold off
    legend('Normal','','','','','Carboplatin - IHC Damage');
    ylim(ylims_t);
    yticks([]);
    
    %setting final attributes
    % han=axes(t_waveform_model_phys,'visible','off'); 
    % han.XLabel.Visible='on';
    xlabel(t_waveform_model_phys,'Time (s)','FontWeight','Bold','FontSize',f_size)
    set(gcf,'Position',fig_dims)
    set(findall(gcf,'-property','FontSize'),'FontSize',f_size);
    set(findall(gcf,'-property','FontWeight'),'FontWeight','Bold');
    
    %Spectral Figures
    
    n_to_plot_mod = n_mod_fft(:,sim_indexes);
    i_to_plot_mod = squeeze(i_mod_fft(ihc,:,sim_indexes));
    
    buff = .1;
    buff = buff.*(1:size(n_to_plot_mod,2));
    
    figure;
    spect_model_phys = tiledlayout(1,2,'TileSpacing','compact');
    nexttile;
    title('Model Simulations')
    hold on
    plot(f,n_to_plot_mod+buff,'color',blck,'linewidth',l_wdth)
    plot(f,i_to_plot_mod+buff,'color',colors_ca,'linewidth',l_wdth)
    hold off
    yticks(buff*1.1);
    yticklabels([labs,'FontWeight','Bold']);
    ytickangle(45)
    ylim(ylims_s);
    xlim([0,2000]);
    %better way to do this?
    legend('Normal','','','','',['C_{ihc} = ',num2str(ihc_val)]);
    
    nexttile;
    title('In-Vivo')
    hold on
    plot(f,phys_scale_fact*n_phys_fft+buff,'color',blck,'linewidth',l_wdth)
    plot(f,phys_scale_fact*i_phys_fft+buff,'color',colors_ca,'linewidth',l_wdth)
    hold off
    legend('Normal','','','','','Carboplatin - IHC Damage');
    ylim(ylims_s);
    xlim([0,2000]);
    yticks([]);
    
    %setting final attributes
    % han=axes(t_waveform_model_phys,'visible','off'); 
    % han.XLabel.Visible='on';
    xlabel(spect_model_phys,'Frequency (Hz)','FontWeight','Bold','FontSize',f_size)
    set(gcf,'Position',fig_dims)
    set(findall(gcf,'-property','FontSize'),'FontSize',f_size);
    set(findall(gcf,'-property','FontWeight'),'FontWeight','Bold');
    
    %Save figs
    exportgraphics(t_waveform_model_phys,[fig_dir,sim_dir,'simulatedWformFig_cihc_',num2str(i),'.png'],'Resolution',300) 
    exportgraphics(spect_model_phys,[fig_dir,sim_dir,'simulatedSpectFig_cihc_',num2str(i),'.png'],'Resolution',300) 
end

%% Figure showing how cihc affects EFR and Spectrum
fig_dims =  [1404 129 641 549];

stim = 12;
norm = repmat(n_mod_efr(:,stim),1,length(mod_data.ihc_grades));
impaired = squeeze(i_mod_efr(:,:,stim))';

buff = 1.5;
buff = buff.*(1:size(norm,2));

cihc_compare_t_fig = tiledlayout(1,1,"TileSpacing","tight");
nexttile;
hold on;
plot(t, norm+buff,'color',blck,'linewidth',l_wdth);
plot(t, flip(impaired,2)+buff,'color',colors_ca,'linewidth',l_wdth);
ylabel('C_{ihc}');
xlabel('Time (s)')
yticks(buff*1.1);
labs = flip(string(mod_data.ihc_grades));
title('SQ25')
set(gcf,'Position',fig_dims)
set(findall(gcf,'-property','FontSize'),'FontSize',f_size);
set(findall(gcf,'-property','FontWeight'),'FontWeight','Bold')
yticklabels([labs,'FontWeight','Bold']);
exportgraphics(cihc_compare_t_fig,[fig_dir,sim_dir,'wform_fxn_of_cihc.png'],'Resolution',300) 

%Spectral Version
norm = repmat(n_mod_fft(:,stim),1,length(mod_data.ihc_grades));
impaired = squeeze(i_mod_fft(:,:,stim))';

buff = .1;
buff = buff.*(1:size(norm,2));

cihc_compare_s_fig = tiledlayout(1,1,"TileSpacing","tight");
nexttile;
hold on;
plot(f, norm+buff,'color',blck,'linewidth',l_wdth);
plot(f, flip(impaired,2)+buff,'color',colors_ca,'linewidth',l_wdth);
ylabel('C_{ihc}');
xlabel('Frequency (Hz)')
yticks(buff*1.1);
labs = flip(string(mod_data.ihc_grades));
title('SQ25')
set(gcf,'Position',fig_dims)
set(findall(gcf,'-property','FontSize'),'FontSize',f_size);
set(findall(gcf,'-property','FontWeight'),'FontWeight','Bold')
yticklabels([labs,'FontWeight','Bold']);

exportgraphics(cihc_compare_s_fig,[fig_dir,sim_dir,'spect_fxn_of_cihc.png'],'Resolution',300) 

%% Making a single spectral figure for a single stim to inset on poster (clean up later)

colors_ca = [0.8500, 0.3250, 0.0980, .85]; %last number is the face alpha/transparency

figure;
sing_spect = tiledlayout(1,1,'TileSpacing','compact');
nexttile;
stim = 2; %sq25
hold on
plot(f,n_phys_fft(:,stim),'Color',blck,'LineWidth',l_wdth);
plot(f,i_phys_fft(:,stim),'Color',colors_ca,'LineWidth',l_wdth)
hold off
yticks([])
title('Physiology');
xlabel('Frequency (Hz)');
set(findall(gcf,'-property','FontSize'),'FontSize',f_size);
set(findall(gcf,'-property','FontWeight'),'FontWeight','Bold')
yticklabels([labs,'FontWeight','Bold']);
set(gcf,'Position',[675 647 488 240]);
exportgraphics(sing_spect,[fig_dir,sim_dir,'phys_spect_inset.png'],'Resolution',300) 


%% Transduction Figures

trans_fxn = tiledlayout(1,1,'TileSpacing','Tight');
nexttile;

cihc_vals = mod_data.ihc_grades;
x = -50:0.1:50;

hold on;
nrm = transduc_nL(x,.1,3);
plot(x,nrm,'Linewidth',4,'Color',blck);
text(10,max(nrm)*.98,'C_{ihc} = 1','Color',blck);
for i = 1:length(cihc_vals)
    to_plot = transduc_nL(x*cihc_vals(i),.1,3);
    colors_ca(4) = colors_ca(4)*.75; 
    plot(x,to_plot,'Linewidth',4,'Color',colors_ca);
    text(25,max(to_plot)*1.01,['C_{ihc} = ',num2str(cihc_vals(i))],'Color',colors_ca);
end
yax = linspace(min(nrm),max(nrm),5);
xax = linspace(min(x),max(x),5);
plot(zeros(length(yax),1),yax,'k--','linewidth',2);
plot(xax,zeros(length(xax),1),'k--','linewidth',2);
title('IHC Transduction Nonlinearity');
set(gca, 'XTick', []);
set(gca, 'YTick', []);%% stim plot
h(1) = xlabel('IHC Input Level');
h(2) = ylabel('IHC Potential (mV)');
set(findall(gcf,'-property','FontSize'),'FontSize',f_size);
set(findall(gcf,'-property','FontWeight'),'FontWeight','Bold');

exportgraphics(trans_fxn,[fig_dir,'transduc_fxn_cihc.png'],'Resolution',300);

