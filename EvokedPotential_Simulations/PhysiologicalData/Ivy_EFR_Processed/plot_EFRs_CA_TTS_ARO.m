clear
close all

load('CA_Processed.mat');
CA = data_out;
load('TTS_Processed.mat');
TTS = data_out;

blck = [0.25, 0.25, 0.25];
colors_ca = [0.8500, 0.3250, 0.0980];
colors_tts = [0, 0.4470, 0.7410];

clear data_out
%% Time and PLV 
%filter signal for ease of visualization
fs = 8e3;
[b,a] = butter(6,[50,800]/(fs/2));

tts_filt_pre(:,1) = filtfilt(b,a,mean(TTS.sam_all_n));
tts_filt_pre(:,2) = filtfilt(b,a,mean(TTS.sq50_all_n));
tts_filt_pre(:,3) = filtfilt(b,a,mean(TTS.sq25_all_n));

tts_plv_pre(:,1) = mean(TTS.plv_base_SAM,2);
tts_plv_pre(:,2) = mean(TTS.plv_base_sq50,2);
tts_plv_pre(:,3) = mean(TTS.plv_base_sq25,2);

tts_filt_post(:,1) = filtfilt(b,a,mean(TTS.sam_all_i));
tts_filt_post(:,2) = filtfilt(b,a,mean(TTS.sq50_all_i));
tts_filt_post(:,3) = filtfilt(b,a,mean(TTS.sq25_all_i));

tts_plv_post(:,1) = mean(TTS.plv_exp_SAM,2);
tts_plv_post(:,2) = mean(TTS.plv_exp_sq50,2);
tts_plv_post(:,3) = mean(TTS.plv_exp_sq25,2);

ca_filt_pre(:,1) = filtfilt(b,a,mean(CA.sam_all_n));
ca_filt_pre(:,2) = filtfilt(b,a,mean(CA.sq50_all_n));
ca_filt_pre(:,3) = filtfilt(b,a,mean(CA.sq25_all_n));

ca_plv_pre(:,1) = mean(CA.plv_base_SAM,2);
ca_plv_pre(:,2) = mean(CA.plv_base_sq50,2);
ca_plv_pre(:,3) = mean(CA.plv_base_sq25,2);

ca_filt_post(:,1) = filtfilt(b,a,mean(CA.sam_all_i));
ca_filt_post(:,2) = filtfilt(b,a,mean(CA.sq50_all_i));
ca_filt_post(:,3) = filtfilt(b,a,mean(CA.sq25_all_i));

ca_plv_post(:,1) = mean(CA.plv_exp_SAM,2);
ca_plv_post(:,2) = mean(CA.plv_exp_sq50,2);
ca_plv_post(:,3) = mean(CA.plv_exp_sq25,2);

tts_filt_pre = tts_filt_pre.*1e6;
tts_filt_post = tts_filt_post.*1e6;

ca_filt_pre = ca_filt_pre.*1e6;
ca_filt_post = ca_filt_post.*1e6;

t = CA.t;
f = CA.f;
%TTS
%Baseline Plots
base_exp_TTS = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');
sam_t = nexttile;
plot(t, tts_filt_pre(:,1),'Color',blck, 'LineWidth',2);
xlim([0.01,.26]);
ylim([-.5,.5]);
title('SAM','FontSize',15);
ylabel('Amplitude (\muV)','FontWeight','Bold','FontSize',13);

sq50_t = nexttile;
plot(t, tts_filt_pre(:,2),'Color',blck, 'LineWidth',2);
xlim([0.01,.26]);
ylim([-.5,.5]);
title('SQ50','FontSize',15);
xlabel('Time (s)','FontSize',13,'FontWeight','Bold');

sq25_t = nexttile;
plot(t, tts_filt_pre(:,3),'Color',blck, 'LineWidth',2);
xlim([0.01,.26]);
ylim([-.5,.5]);
title('SQ25','FontSize',15);

% PLV
sam_plv = nexttile;
plot(f, tts_plv_pre(:,1),'Color',blck, 'LineWidth',2);
xlim([0,4e3]);
ylim([0,1]);
ylabel('PLV','FontWeight','Bold','FontSize',13)

sq50_plv = nexttile;
plot(f, tts_plv_pre(:,2),'Color',blck, 'LineWidth',2);
xlim([0,4e3]);
ylim([0,1]);

sq25_plv = nexttile;
plot(f, tts_plv_pre(:,3),'Color',blck, 'LineWidth',2);
xlim([0,4e3]);
ylim([0,1]);

base_exp_TTS.XLabel.String = 'Frequency (Hz)';
base_exp_TTS.XLabel.FontWeight = 'Bold';
base_exp_TTS.XLabel.FontSize = 13;

set(gcf,'Position',[1228 32 1322 705])
exportgraphics(base_exp_TTS,'TTS_PLV_Base_Ivy.png','Resolution',300)

%Exposed Plots
axes(sam_t);
hold on
plot(t, tts_filt_post(:,1),'Color',colors_tts, 'LineWidth',2);
hold off
axes(sq50_t);
hold on
plot(t, tts_filt_post(:,2),'Color',colors_tts, 'LineWidth',2);
hold off
axes(sq25_t);
hold on
plot(t, tts_filt_post(:,3),'Color',colors_tts, 'LineWidth',2);
hold off

%Plv
axes(sam_plv);
hold on
plot(f, tts_plv_post(:,1),'Color',colors_tts, 'LineWidth',2);
hold off

axes(sq50_plv);
hold on
plot(f, tts_plv_post(:,2),'Color',colors_tts, 'LineWidth',2);
hold off

axes(sq25_plv);
hold on
plot(f, tts_plv_post(:,3),'Color',colors_tts, 'LineWidth',2);
hold off
legend('Pre','Post','Fontsize',11,'Location','NorthEast')

exportgraphics(base_exp_TTS,'TTS_PLV_Base_Exp_Ivy.png','Resolution',300)

%CA

base_exp_CA = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');
sam_t = nexttile;
plot(t, ca_filt_pre(:,1),'Color',blck, 'LineWidth',2);
xlim([0.01,.26]);
ylim([-.5,.5]);
title('SAM','FontSize',15);
ylabel('Amplitude (\muV)','FontWeight','Bold','FontSize',13);

sq50_t = nexttile;
plot(t, ca_filt_pre(:,2),'Color',blck, 'LineWidth',2);
xlim([0.01,.26]);
ylim([-.5,.5]);
title('SQ50','FontSize',15);
xlabel('Time (s)','FontSize',13,'FontWeight','bold');

sq25_t = nexttile;
plot(t, ca_filt_pre(:,3),'Color',blck, 'LineWidth',2);
xlim([0.01,.26]);
ylim([-.5,.5]);
title('SQ25','FontSize',15);

% PLV
sam_plv = nexttile;
plot(f, ca_plv_pre(:,1),'Color',blck, 'LineWidth',2);
xlim([0,4e3]);
ylim([0,1]);
ylabel('PLV','FontWeight','Bold','FontSize',13)

sq50_plv = nexttile;
plot(f, ca_plv_pre(:,2),'Color',blck, 'LineWidth',2);
xlim([0,4e3]);
ylim([0,1]);

sq25_plv = nexttile;
plot(f, ca_plv_pre(:,3),'Color',blck, 'LineWidth',2);
xlim([0,4e3]);
ylim([0,1]);

base_exp_CA.XLabel.String = 'Frequency (Hz)';
base_exp_CA.XLabel.FontWeight = 'Bold';
base_exp_CA.XLabel.FontSize = 13;

set(gcf,'Position',[1228 32 1322 705])
exportgraphics(base_exp_CA,'CA_PLV_Base_Ivy.png','Resolution',300)

%Exposed Plots
axes(sam_t);
hold on
plot(t, ca_filt_post(:,1),'Color',colors_ca, 'LineWidth',2);
hold off
axes(sq50_t);
hold on
plot(t, ca_filt_post(:,2),'Color',colors_ca, 'LineWidth',2);
hold off
axes(sq25_t);
hold on
plot(t, ca_filt_post(:,3),'Color',colors_ca, 'LineWidth',2);
hold off

%Plv
axes(sam_plv);
hold on
plot(f, ca_plv_post(:,1),'Color',colors_ca, 'LineWidth',2);
hold off

axes(sq50_plv);
hold on
plot(f, ca_plv_post(:,2),'Color',colors_ca, 'LineWidth',2);
hold off

axes(sq25_plv);
hold on
plot(f, ca_plv_post(:,3),'Color',colors_ca, 'LineWidth',2);
hold off
legend('Pre','Post','Fontsize',11,'Location','NorthEast')

exportgraphics(base_exp_CA,'CA_PLV_Base_Exp_Ivy.png','Resolution',300)

%% EFR Demo fig
purp = [0.4940 0.1840 0.5560];
gren = [0.4660 0.6740 0.1880];
[~,~,n_floor_pre] = getPeaks(f,ca_plv_pre(:,3),100,16);
[~,~,n_floor_post] = getPeaks(f,ca_plv_post(:,3),100,16);

n_floor_dem = mean(n_floor_pre+n_floor_post)/2;

figure;
hold on;
plot(f,ca_plv_pre(:,3),'linewidth',1.5,'Color',blck)
plot(f,ca_plv_post(:,3),'linewidth',1.5,'Color',colors_ca)
plot([5:240],.1*ones(1,236),'color',gren,'linewidth',5);
text(5,.05,'Low','color',gren,'fontweight','bold');
plot([280:1600],.1*ones(1,1321),'color',purp,'linewidth',5);
text(700,.05,'High','color',purp,'fontweight','bold');
rectangle('Position',[0,0,4000,n_floor_dem],'EdgeColor',[blck],'FaceColor',[blck,0.05],'LineWidth',1.5);
text(3700,.13,'floor','color',blck,'fontweight','bold');
ylim([0,1])
title('Carboplatin PLV - SQ25')
ylabel('PLV')
xlabel('Frequency(Hz)')
legend('Baseline','Exposed','FontSize',12)

set(gcf,'Position',[1540 430 727 531])
exportgraphics(gcf,'demo_plv_ratio.png','Resolution',300);

%% Plotting Ratios
R_plvs = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');
x_bar = [1,2];
ybound = [1,6.65];
ytick = 1:6;

nexttile;
hold on
errorbar(x_bar(1),mean(TTS.pre_post_SAM_rat(:,1)),std(TTS.pre_post_SAM_rat(:,1))./sqrt(4),'s','Color',blck,'MarkerEdgeColor',blck,'MarkerFaceColor',blck,'MarkerSize',11,'Linewidth',1.5)
errorbar(x_bar(2),mean(TTS.pre_post_SAM_rat(:,2)),std(TTS.pre_post_SAM_rat(:,2))./sqrt(4),'s','Color',colors_tts,'MarkerEdgeColor',colors_tts,'MarkerFaceColor',colors_tts,'MarkerSize',11,'Linewidth',1.5)
title('SAM','FontSize',15);
legend('Baseline','Synaptopathy','FontSize',11,'Location','NorthWest')
xlim([0,3]);
xticks([1,2]);
yticks(ytick);
xticklabels([]);

set(gca,'FontSize',12)
ylim(ybound);
grid on
hold off

nexttile;
hold on
errorbar(x_bar(1),mean(TTS.pre_post_sq50_rat(:,1)),std(TTS.pre_post_sq50_rat(:,1))./sqrt(4),'s','Color',blck,'MarkerEdgeColor',blck,'MarkerFaceColor',blck,'MarkerSize',11,'Linewidth',1.5)
errorbar(x_bar(2),mean(TTS.pre_post_sq50_rat(:,2)),std(TTS.pre_post_sq50_rat(:,2))./sqrt(4),'s','Color',colors_tts,'MarkerEdgeColor',colors_tts,'MarkerFaceColor',colors_tts,'MarkerSize',11,'Linewidth',1.5)
title('SQ50','FontSize',15);
xlim([0,3]);
ylim(ybound);
yticklabels([])
xticks([1,2]);
xticklabels([]);
yticks(ytick);

set(gca,'FontSize',12)
grid on
hold off

nexttile;
hold on
errorbar(x_bar(1),mean(TTS.pre_post_sq25_rat(:,1)),std(TTS.pre_post_sq25_rat(:,1))./sqrt(4),'s','Color',blck,'MarkerEdgeColor',blck,'MarkerFaceColor',blck,'MarkerSize',11,'Linewidth',1.5)
errorbar(x_bar(2),mean(TTS.pre_post_sq25_rat(:,2)),std(TTS.pre_post_sq25_rat(:,2))./sqrt(4),'s','Color',colors_tts,'MarkerEdgeColor',colors_tts,'MarkerFaceColor',colors_tts,'MarkerSize',11,'Linewidth',1.5)
title('SQ25','FontSize',15);
xlim([0,3]);
xticks([1,2]);
ylim(ybound);
yticklabels([])
set(gca,'FontSize',12)
xticklabels([]);
yticks(ytick);

grid on
hold off

R_plvs.YLabel.String = 'R_{PLV}';
R_plvs.YLabel.FontWeight = 'bold';
R_plvs.YLabel.FontSize = 13;

R_plvs.Title.String = 'Ratio of Upper to Lower Harmonic Sums';
R_plvs.Title.FontWeight = 'bold';
R_plvs.Title.FontSize = 14;

exportgraphics(R_plvs,'PLV_ratios_TTS.png','Resolution',600);

nexttile;
hold on
errorbar(x_bar(1),mean(CA.pre_post_SAM_rat(:,1)),std(CA.pre_post_SAM_rat(:,1))./sqrt(4),'s','Color',blck,'MarkerEdgeColor',blck,'MarkerFaceColor',blck,'MarkerSize',11,'Linewidth',1.5)
errorbar(x_bar(2),mean(CA.pre_post_SAM_rat(:,2)),std(CA.pre_post_SAM_rat(:,2))./sqrt(4),'s','Color',colors_ca,'MarkerEdgeColor',colors_ca,'MarkerFaceColor',colors_ca,'MarkerSize',11,'Linewidth',1.5)
xlim([0,3]);
xticks([1,2]);
yticks(ytick);
xticklabels(["Pre","Post"]);
set(gca,'FontSize',12)
legend('Baseline','IHC Damage','FontSize',11,'Location','NorthWest')
ylim(ybound);
grid on
hold off

nexttile;
hold on
errorbar(x_bar(1),mean(CA.pre_post_sq50_rat(:,1)),std(CA.pre_post_sq50_rat(:,1))./sqrt(4),'s','Color',blck,'MarkerEdgeColor',blck,'MarkerFaceColor',blck,'MarkerSize',11,'Linewidth',1.5)
errorbar(x_bar(2),mean(CA.pre_post_sq50_rat(:,2)),std(CA.pre_post_sq50_rat(:,2))./sqrt(4),'s','Color',colors_ca,'MarkerEdgeColor',colors_ca,'MarkerFaceColor',colors_ca,'MarkerSize',11,'Linewidth',1.5)
xlim([0,3]);
ylim(ybound);
yticklabels([])
xticks([1,2]);
yticks(ytick);
xticklabels(["Pre","Post"]);
set(gca,'FontSize',12)
grid on
hold off

nexttile;
hold on
errorbar(x_bar(1),mean(CA.pre_post_sq25_rat(:,1)),std(CA.pre_post_sq25_rat(:,1))./sqrt(4),'s','Color',blck,'MarkerEdgeColor',blck,'MarkerFaceColor',blck,'MarkerSize',11,'Linewidth',1.5)
errorbar(x_bar(2),mean(CA.pre_post_sq25_rat(:,2)),std(CA.pre_post_sq25_rat(:,2))./sqrt(4),'s','Color',colors_ca,'MarkerEdgeColor',colors_ca,'MarkerFaceColor',colors_ca,'MarkerSize',11,'Linewidth',1.5)
xlim([0,3]);
xticks([1,2]);
ylim(ybound);
yticklabels([])
yticks(ytick);
xticklabels(["Pre","Post"]);
set(gca,'FontSize',12)
grid on
hold off

exportgraphics(R_plvs,'PLV_ratios_TTS_CA.png','Resolution',300);

%% Only CA:
R_plvs_CA_only = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');

nexttile;
hold on
errorbar(x_bar(1),mean(CA.pre_post_SAM_rat(:,1)),std(CA.pre_post_SAM_rat(:,1))./sqrt(4),'s','Color',blck,'MarkerEdgeColor',blck,'MarkerFaceColor',blck,'MarkerSize',11,'Linewidth',1.5)
errorbar(x_bar(2),mean(CA.pre_post_SAM_rat(:,2)),std(CA.pre_post_SAM_rat(:,2))./sqrt(4),'s','Color',colors_ca,'MarkerEdgeColor',colors_ca,'MarkerFaceColor',colors_ca,'MarkerSize',11,'Linewidth',1.5)
xlim([0,3]);
xticks([1,2]);
yticks(ytick);
xticklabels(["Pre","Post"]);
set(gca,'FontSize',12)
legend('Baseline','IHC Damage','FontSize',11,'Location','NorthWest')
ylim(ybound);
grid on
hold off

nexttile;
hold on
errorbar(x_bar(1),mean(CA.pre_post_sq50_rat(:,1)),std(CA.pre_post_sq50_rat(:,1))./sqrt(4),'s','Color',blck,'MarkerEdgeColor',blck,'MarkerFaceColor',blck,'MarkerSize',11,'Linewidth',1.5)
errorbar(x_bar(2),mean(CA.pre_post_sq50_rat(:,2)),std(CA.pre_post_sq50_rat(:,2))./sqrt(4),'s','Color',colors_ca,'MarkerEdgeColor',colors_ca,'MarkerFaceColor',colors_ca,'MarkerSize',11,'Linewidth',1.5)
xlim([0,3]);
ylim(ybound);
yticklabels([])
xticks([1,2]);
yticks(ytick);
xticklabels(["Pre","Post"]);
set(gca,'FontSize',12)
grid on
hold off

nexttile;
hold on
errorbar(x_bar(1),mean(CA.pre_post_sq25_rat(:,1)),std(CA.pre_post_sq25_rat(:,1))./sqrt(4),'s','Color',blck,'MarkerEdgeColor',blck,'MarkerFaceColor',blck,'MarkerSize',11,'Linewidth',1.5)
errorbar(x_bar(2),mean(CA.pre_post_sq25_rat(:,2)),std(CA.pre_post_sq25_rat(:,2))./sqrt(4),'s','Color',colors_ca,'MarkerEdgeColor',colors_ca,'MarkerFaceColor',colors_ca,'MarkerSize',11,'Linewidth',1.5)
xlim([0,3]);
xticks([1,2]);
ylim(ybound);
yticklabels([])
yticks(ytick);
xticklabels(["Pre","Post"]);
set(gca,'FontSize',12)
grid on
hold off
exportgraphics(R_plvs_CA_only,'PLV_ratios_TTS_CA.png','Resolution',300);


%% Combined Time/Spectral figures:

%Carboplatin
purp = [0.4940 0.1840 0.5560];
gren = [0.4660 0.6740 0.1880];
[~,~,n_floor_pre] = getPeaks(f,ca_plv_pre(:,3),100,16);
[~,~,n_floor_post] = getPeaks(f,ca_plv_post(:,3),100,16);

n_floor_dem = mean(n_floor_pre+n_floor_post)/2;

figure;
% subplot(1,2,1)
hold on;
plot(f,ca_plv_pre(:,3),'linewidth',1.5,'Color',blck)
plot(f,ca_plv_post(:,3),'linewidth',1.5,'Color',colors_ca)
plot([5:240],.1*ones(1,236),'color',gren,'linewidth',5);
text(5,.05,'Low','color',gren,'fontweight','bold');
plot([280:1600],.1*ones(1,1321),'color',purp,'linewidth',5);
text(700,.05,'High','color',purp,'fontweight','bold');
rectangle('Position',[0,0,4000,n_floor_dem],'EdgeColor',[blck],'FaceColor',[blck,0.05],'LineWidth',1.5);
text(3200,.13,'floor','color',blck,'fontweight','bold');
hold off;
ylim([0,1])
title('IHC Damage - SQ25')
ylabel('PLV','FontWeight','bold')
xlabel('Frequency(Hz)','FontWeight','bold')
legend('Baseline','CA','FontSize',12,'Location','Northwest')

xstart = .6;
xend = .9;
ystart = 0.6;
yend = .9;

axes('Position',[xstart ystart xend-xstart yend-ystart])
box on
hold on
plot(t, ca_filt_pre(:,3),'Color',blck, 'LineWidth',2);
plot(t, ca_filt_post(:,3),'Color',colors_ca, 'LineWidth',2);
xlim([0.1,.26]);
ylim([-.5,.6]);
yticks([-.5,0,.5])
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude \muV','FontWeight','bold')
hold off
set(gcf,'Position',[1557 538 560 420])
print('CA_PLV_TIME','-dsvg');

figure;
hold on;
plot(f,tts_plv_pre(:,3),'linewidth',1.5,'Color',blck)
plot(f,tts_plv_post(:,3),'linewidth',1.5,'Color',colors_tts)
plot([5:240],.1*ones(1,236),'color',gren,'linewidth',5);
text(5,.05,'Low','color',gren,'fontweight','bold');
plot([280:1600],.1*ones(1,1321),'color',purp,'linewidth',5);
text(700,.05,'High','color',purp,'fontweight','bold');
rectangle('Position',[0,0,4000,n_floor_dem],'EdgeColor',[blck],'FaceColor',[blck,0.05],'LineWidth',1.5);
text(3200,.13,'floor','color',blck,'fontweight','bold');
hold off;
ylim([0,1])
title('Synaptopathy - SQ25')
ylabel('PLV','FontWeight','bold')
xlabel('Frequency(Hz)','FontWeight','bold')
legend('Baseline','TTS','FontSize',12,'Location','Northwest')
axes('Position',[xstart ystart xend-xstart yend-ystart])
box on
hold on
plot(t, tts_filt_pre(:,3),'Color',blck, 'LineWidth',2);
plot(t, tts_filt_post(:,3),'Color',colors_tts, 'LineWidth',2);
xlim([0.1,.26]);
ylim([-.5,.6]);
yticks([-.5,0,.5])
xlabel('Time(s)','FontWeight','bold');
ylabel('Amplitude \muV','FontWeight','bold')
hold off
set(gcf,'Position',[1557 538 560 420])
print('TTS_PLV_TIME','-dsvg');
% set(gcf,'Position',[1540 430 727 531])

%% Some stats:

% %pooling everything
% CA_pre_all = [CA.pre_post_SAM_rat(:,1);CA.pre_post_sq25_rat(:,1);CA.pre_post_sq50_rat(:,1)];
% CA_post_all = [CA.pre_post_SAM_rat(:,2);CA.pre_post_sq25_rat(:,2);CA.pre_post_sq50_rat(:,2)];
% 
% TTS_pre_all = [TTS.pre_post_SAM_rat(:,1);TTS.pre_post_sq25_rat(:,1);TTS.pre_post_sq50_rat(:,1)];
% TTS_post_all = [TTS.pre_post_SAM_rat(:,2);TTS.pre_post_sq25_rat(:,2);TTS.pre_post_sq50_rat(:,2)];
% 
% [h_CA,p_CA,~,~] = ttest(CA_pre_all,CA_post_all)
% [h_TTS,p_TTS,~,~] = ttest(TTS_pre_all,TTS_post_all)


