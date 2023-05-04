function [pin]=SAM_RAM(n_pulse,gap,fm,fc,md,amp)
% n_pulse:number of pulse of SAM
% gap:gap duration in seconds
% fm:modulation freq in Hz
% fc:carrier freq in Hz
% md:modulation depth (between 0 and 1)
% amp:amplitude (not stimdb!!)
% comment 2 lines after aa if you dint want ramp.
Fs=100e3;
dur=1/fm;
rt=2.5e-3; 
irpts = rt*Fs;
t=0:1/Fs:dur-1/Fs;
aa=amp*sin(2*pi*fc*t).*(1+md*sin(2*pi*fm*t-pi/2));
aa(1:irpts)=aa(1:irpts).*(0:(irpts-1))/irpts;
aa(length(aa)-irpts:length(aa))=aa(length(aa)-irpts:length(aa)).*(irpts:-1:0)/irpts;
gp_s=zeros(1,gap*Fs);
pin=[repmat([aa gp_s],1,n_pulse-1) aa];



