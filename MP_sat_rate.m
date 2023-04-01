load('nodblossdata.mat')
no_loss=cellfun(@(w) squeeze(sum(w(:,onbin:end,:),2))./0.6,all_spks,'UniformOutput',false);
clearvars -except no_loss
load('highfreqdblossdata.mat')
hi_loss=cellfun(@(w) squeeze(sum(w(:,onbin:end,:),2))./0.6,all_spks,'UniformOutput',false);
clearvars -except no_loss hi_loss
load('lowfreqdblossdata.mat')
lo_loss=cellfun(@(w) squeeze(sum(w(:,onbin:end,:),2))./0.6,all_spks,'UniformOutput',false);
%all_loss=cellfun(@(w,x,y) reshape([w x y],140,3),no_loss,hi_loss,lo_loss,'UniformOutput',false);
clearvars -except no_loss hi_loss lo_loss
low_loss=cellfun(@(w) squeeze(sum(w(:,onbin:end,:),2))./0.6,all_spks,'UniformOutput',false);
all_loss=cellfun(@(w,x,y,z) reshape([w x y z],140,4),no_loss,hi_loss,lo_loss,low_loss,'UniformOutput',false);

cls='rgbm';
figure
for ii=1:length(ag_fs)
    a1=all_loss(ii,:);
    subplot(2,4,ii)
    for jj=1:4
        a2=cell2mat(cellfun(@(w) w(:,jj),a1,'UniformOutput',false));
        errorbar(stimdb,nanmean(a2),nanstd(a2)./sqrt(140),cls(jj))
        hold on
    end
    xlabel('Intensity(dbSPL)')
    ylabel('driven rate')
    ylim([80 190])
    title(ag_fs(ii))
end
   







