clear
f1 = dir('*.mat');
nm={f1.name};
%cihc=zeros(1,19);
all_rt=zeros(21,21,19);  %freq,db,cihc vals

for ii=1:size(nm,2)
    ii
    
    spk=load(nm{1,ii});
    %     if ii==5
    %         aa=spk.all_spk;
    %     else
    %         aa=spk.kk;
    %     end
    aa=spk.x(ii,:);
    clearvars spk all_spk
    %cihc=[cihc str2double(nm{1,ii}(17:20))];
    ab=cellfun(@(w) nanmean(squeeze(sum(w(:,1001:end,:),2))./0.6),aa,'UniformOutput',false);
    all_rt(:,:,ii)=cell2mat(ab');
    clearvars ab aa
end


clrs='bycbycbycbycbycbyck';
figure
for ii=1:19
    ak=all_rt(:,:,ii);
    errorbar([1:21],nanmean(ak),nanstd(ak)./sqrt(size(ak,1)),clrs(ii));
    hold on
end












