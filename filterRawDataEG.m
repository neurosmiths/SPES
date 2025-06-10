% [20250220] PLotting some raw data from the SPES session for 202407 in
% order to appease reviewers at Brain Stimulation. 

fName = 'D:\Data\UIC202407\SOZStim\20240412-152907\20240412-152907.ns6';

NSX = openNSx('read',fName,'c:18');

[b,a] = butter(4,300./(3e4/2),'high');

hpLFP = filtfilt(b,a,double(NSX.Data));



tSec = linspace(0,length(hpLFP)./3e4,length(hpLFP));

plot(tSec,hpLFP,'k')