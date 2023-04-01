ag_fs = [125 250 500 1e3 2e3 4e3 8e3];
%ag_dbloss = [0 0 0 20 40 60 80];
 ag_dbloss = [0 0 0 0 0 0 0];
%ag_dbloss = [80 60 40 20 0 0 0];
species = 2; % Human cochlear tuning (Shera et al.,  2002)
Fs = 100e3;
fc=ag_fs;
fm=50;
fs=Fs;
%pure tone, 1000Hz, 70dB/30dB, histogram(normal), 0dB: uniform distribution(period is fc(carrier))
%sam, 1000Hz, fm= 100, uniform for 0(till fm), phase locking present at around 30dB, uniform line again at around 70dB 
%compute synchrony vs level fucntion(like rate level curve: saturate and drop a little bit), sam tone: like sine wave: increase then decrease (max at 25-30dB above threshold)
%jorris 1992; sam tone modulation data curves
spont = 100;   % spontaneous firing rate
tabs   = 0.6e-3; % Absolute refractory period
trel   = 0.6e-3;
modDepth=[0,0.0313,0.0625,0.125,0.25,0.5,1];
dur=0.6;
ondelay = 10e-3;
rt = 2.5e-3; % rise/fall time in seconds
t = 0:1/Fs:dur-1/Fs; % time vector
mxpts = length(t);
irpts = rt*Fs;
onbin = round(ondelay*Fs);
phi_c=0;
phi_m=0;
stimdb=10:5:40;
amp=sqrt(2)*20e-6*10.^(stimdb/20);
all_stim=cell(length(amp),length(ag_fs));  %2d cell array
rep_n=20;
%[cohcs,cihcs,OHC_Loss]=fitaudiogram2(ag_fs,ag_dbloss,2);
cohcs=[1,1,1,1,1,1,1];
cihcs=[1,1,1,1,1,1,1];


for ii=1:length(ag_fs)  %rows: cf
    for kk=1:length(amp) %cols: stimdb
        all_stim{ii,kk}=[];
        for jj=1:length(modDepth)  %for all mod depths corresponding to each cf and stimdb cell array
            pin = zeros(1,onbin+mxpts);
            %pin(onbin+1:onbin+mxpts) = amp(kk)*sin(2*pi*ag_fs(ii)*t); % unramped pure tone
            %pin(onbin+1:onbin+mxpts) = 0;  % zero stim

            pin(onbin+1:onbin+mxpts) = amp(kk)*sin(2*pi*ag_fs(ii)*t).*(1+modDepth(jj)*sin(2*pi*fm*t));  % unramped sam, first 1ms is silence
            pin(onbin+1:onbin+irpts)= pin(onbin+1:onbin+irpts).*(0:(irpts-1))/irpts; %ramping
            pin(onbin+(mxpts-irpts):onbin+mxpts)=pin(onbin+(mxpts-irpts):onbin+mxpts).*(irpts:-1:0)/irpts;
            %[sig, time]= create_SAM(fc(ii), fm, fs, modDepth(jj), dur, amp(kk), phi_c, phi_m);
            all_stim{ii,kk}=[all_stim{ii,kk};pin];  %stimulus
        end
    end
end
dt=1/Fs;
all_spks=cell(length(amp),length(ag_fs)); 

for ii=1:length(ag_fs)
    for kk=1:length(amp)
        cut=zeros(rep_n,size(all_stim{1,1},2),length(modDepth));
        for ll=1:length(modDepth)
            for jj=1:rep_n
                vihc = model_IHC_BEZ2018(all_stim{ii,kk}(ll,:),ag_fs(ii),1,dt,4*dur,cohcs(ii),cihcs(ii),2); 
                [psth,meanrate, varrate, synout, trd_vector,trel_vector] = model_Synapse_BEZ2018(vihc,ag_fs(ii),1,dt,1,0,spont,tabs,trel);
                cut(jj,:,ll)=psth(1:size(all_stim{1,1},2)); %spike train length taken till length of stimulus
            end
        end
        all_spks{ii,kk}=cut;
    end
end

all_vs=cell(length(ag_fs),length(amp));
for ii=1:length(ag_fs)
    for jj=1:length(amp)
        ap1=[];
        av=all_spks{ii,jj}; %spk trains for diff mod depth
        for kk=2:size(av,3) %for mod depths
            av1=av(:,onbin+1:end,kk); %spk train for time where stim present
            all_phs=cell(1,size(av1,1));
            for ll=1:size(av1,1) %for 20 reps
                spt=t(find(av1(ll,:)==1)); %find times at which spike occured
                ang=2*pi*fm*spt;  %to angles
                phs=ang-floor(ang/(2*pi))*(2*pi); %phase rounded to nearest 2*pi multiple
                all_phs{1,ll}=phs;
            end
            VS_t=cellfun(@(w) sqrt(sum(cos(w))^2+sum(sin(w))^2)/length(w),all_phs);%vs for each iter
            phi_t=cellfun(@(w) atan((sum(sin(w)))/sum(cos(w))),all_phs); %trial wise, for each iter(20)
            phi_c=atan(sum(sin(cell2mat(all_phs)))/sum(cos(cell2mat(all_phs)))); %mean phase for all trials 
            VS_pp=VS_t.*cos(phi_t-phi_c);
            ap1=[ap1;VS_pp];
        end
        all_vs{ii,jj}=ap1;
    end
end
mn=cellfun(@(w) [nanmean(w,2) nanstd(w,[],2)./sqrt(size(w,2))],all_vs,'UniformOutput',false); %std error value
figure  %for vs_pp, x axis: mod depths, row: cfs, col: intensities
kp=1;
for ii=1:7
    for jj=1:7
        subplot(7,7,kp)
        errorbar([1:6],mn{ii,jj}(:,1),mn{ii,jj}(:,2)); % sem(standard error of mean) across 20 repitions
        ylim([-0.02 0.5])
        kp=kp+1;
    end
end
all_dp=cell(length(ag_fs),length(amp));
for ii=1:length(ag_fs)  
    for jj=1:length(amp)
        ap1=[];
        av=all_spks{ii,jj};
        mu_n=nanmean(squeeze(sum(av(:,1:onbin,:),2))./ondelay);
        mu_s=nanmean(squeeze(sum(av(:,onbin+1:end,:),2))./dur);
        st_n=nanstd(squeeze(sum(av(:,1:onbin,:),2))./ondelay);
        st_s=nanstd(squeeze(sum(av(:,onbin+1:end,:),2))./dur);
        dp=(mu_s-mu_n)./sqrt(st_n.^2+st_s.^2);
        all_dp{ii,jj}=dp;
    end
end
figure  
kp=1;
for ii=1:7
    for jj=1:7
        subplot(7,7,kp)
        plot(all_dp{ii,jj});
        kp=kp+1;
    end
end




            
                
                
            





                
                
