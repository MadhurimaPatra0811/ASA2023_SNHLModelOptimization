load('all10dbloss.mat')
mn1=cellfun(@(w) [nanmean(w,2) nanstd(w,[],2)./sqrt(size(w,2))],all_vs,'UniformOutput',false); %std error value
clearvars -except mn1

load('all40dbloss.mat')
mn2=cellfun(@(w) [nanmean(w,2) nanstd(w,[],2)./sqrt(size(w,2))],all_vs,'UniformOutput',false); %std error value
clearvars -except mn2 mn1

load('nodblossdata.mat')
mn3=cellfun(@(w) [nanmean(w,2) nanstd(w,[],2)./sqrt(size(w,2))],all_vs,'UniformOutput',false); %std error value

figure  %for vs_pp, x axis: mod depths, row: cfs, col: intensities
kp=1;
for ii=1:7
    for jj=1:7
        subplot(7,7,kp)
        hold on
        errorbar([1:6],mn1{ii,jj}(:,1),mn1{ii,jj}(:,2));
        errorbar([1:6],mn2{ii,jj}(:,1),mn2{ii,jj}(:,2));
        errorbar([1:6],mn3{ii,jj}(:,1),mn3{ii,jj}(:,2));
        kp=kp+1;
        ylim([-.02 0.6])
    end
end