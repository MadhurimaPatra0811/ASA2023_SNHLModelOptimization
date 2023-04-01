load('cihc50.mat')  % spont rate,evoked rate comparison, 
clearvars -except all_spks
bs_dp=cell(7,7);
nbs=5000;
for ii=1:7 % 7 freqs
    for jj=1:7 %7 intensity
        [ii jj]
        for kk=1:7 %7 mod depth
            ap=all_spks{ii,jj}(:,:,kk); % 20x6100(data points)
            mns=fix(sum(nanmean(ap))); %bs will take no of spikes=mean spike rate 
            fnf=find(ap==1); % time index where spike occured, runs columnwise
            bs_rts=[];
            for nn=1:nbs
                aq=datasample(fnf,mns,'Replace',false); %mns no of data points taken 
                [row, col]=ind2sub(size(ap),aq);%indexes of fnf changed to get column(time posititon in matrix)
                zrs=zeros(1,size(ap,2)); %synthetic spike train(1x6100)
                zrs(1,col)=1;% replace 1 from col values
                bs_rts=[bs_rts;[sum(zrs(1:1000)) sum(zrs(1001:end))]]; %1:1000(1ms) is spont rate, remaining is evoked rate
            end
            dp=(nanmean(bs_rts(:,2))-nanmean(bs_rts(:,1)))/sqrt(nanmean([nanvar(bs_rts(:,2)) nanvar(bs_rts(:,1))]));
            bs_dp{ii,jj}=[bs_dp{ii,jj} dp];
        end
    end
end
bs1=cell2mat(bs_dp(:)); 
figure
errorbar([1:7],nanmean(bs1),nanstd(bs1)./sqrt(size(bs1,1)));
                
