function spkk = getspktime(spk)
% get spike time from all_apks

%     spk = all_spks{1,1};
    spks= spk(:,1001:61000,5);  % spk train for 1 mod depth
    spk1=zeros(size(spks,1),fix(size(spks,1)/100));
    
    for ii=1:size(spks,2)/100
        spk1(:,ii)= sum(spks(:,ii*100-99:ii*100),2);
    end
       
   f1=find(spk1==2);
   spk1(f1)=1;
        
        
        
      spkk=[];
    
    for ii=1:size(spk1,1)
        nsp=find(spk1(ii,:)==1);
        spkt=[];
        spkt(1:length(nsp),1)=ii;
        spkt(1:length(nsp),2)=nsp;
        spkk=[spkk;spkt];
    end
end
    
        

        
        