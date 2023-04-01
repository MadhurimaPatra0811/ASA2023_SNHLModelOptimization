ag_fs = [125 250 500 1e3 2e3 4e3 8e3];
ag_dbloss = [0 0 0 20 40 60 80];
species = 2; % Human cochlear tuning (Shera et al.,  2002)
Fs = 100e3;
fc=ag_fs;
fm=50;
fs=Fs;
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
[cohcs,cihcs,OHC_Loss]=fitaudiogram2(ag_fs,ag_dbloss,2);
for ii=1:length(ag_fs)
    for kk=1:length(amp)
        all_stim{ii,kk}=[];
        for jj=1:length(modDepth)
            pin = zeros(1,onbin+mxpts);
            pin(onbin+1:onbin+mxpts) = amp(kk)*sin(2*pi*ag_fs(ii)*t).*(1+modDepth(jj)*sin(2*pi*fm*t)); % unramped stimulus
            pin(onbin+1:onbin+irpts)= pin(onbin+1:onbin+irpts).*(0:(irpts-1))/irpts;
            pin(onbin+(mxpts-irpts):onbin+mxpts)=pin(onbin+(mxpts-irpts):onbin+mxpts).*(irpts:-1:0)/irpts;
            %[sig, time]= create_SAM(fc(ii), fm, fs, modDepth(jj), dur, amp(kk), phi_c, phi_m);
            all_stim{ii,kk}=[all_stim{ii,kk};pin];
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
                cut(jj,:,ll)=psth(1:size(all_stim{1,1},2));
            end
        end
        all_spks{ii,kk}=cut;
    end
end

                
                
