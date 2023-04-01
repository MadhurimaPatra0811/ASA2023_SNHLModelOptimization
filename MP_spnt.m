ag_fs = [125 250 500 1e3 2e3 4e3 8e3];
ag_dbloss(1,:) = [0 0 0 20 40 60 80];
ag_dbloss(2,:) = [80 60 40 20 0 0 0];
ag_dbloss(3,:) = [0 0 0 0 0 0 0];
Fs = 100e3;dt=1/Fs;
spont = .4;   % spontaneous firing rate
tabs   = 0.6e-3; % Absolute refractory period
trel   = 0.6e-3;
dur=0.6;
%[cohcs,cihcs,OHC_Loss]=fitaudiogram2(ag_fs,ag_dbloss,2);
all_psth=cell(3,7);
for kk=1:3
    [cohcs,cihcs,OHC_Loss]=fitaudiogram2(ag_fs,ag_dbloss(kk,:),2);
    for ii=1:7
        for jj=1:20
            
            empty_stim=zeros(1,61000);
            vihc = model_IHC_BEZ2018(empty_stim,ag_fs(ii),1,dt,4*dur,cohcs(ii),cihcs(ii),2);
            [psth,meanrate, varrate, synout, trd_vector,trel_vector] = model_Synapse_BEZ2018(vihc,ag_fs(ii),1,dt,1,0,spont,tabs,trel);
            all_psth{kk,ii}=[all_psth{kk,ii};psth(1:61000)];
        end
    end
end

aa=cell2mat(cellfun(@(w) (sum(w(:,1001:end),2)./0.6)',all_psth,'UniformOutput',false));
aa=aa';
bb=[aa(:) [ones(140,1);2*ones(140,1);3*ones(140,1)]];

ecdf(aa')
aa=cellfun(@(w) sum(reshape(nanmean(w),200,61000/200))./0.002,all_psth,'UniformOutput',false);
figure
for ii=1:7
    subplot(2,4,ii)
    plot(aa{1,ii})
end