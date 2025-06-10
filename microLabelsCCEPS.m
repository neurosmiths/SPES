function [microLabels,microPts,SOZ] = microLabelsCCEPS(ptID)

% MICROLABELSCCEPS outputs a list of patients with microwires and 
%   channel labels for microwires
%
% if you just want a list of micro patients, use 'NaP' as the input arg.

% author: EHS20200624
% added atlas labels: EHS20220217
%changed for CCEPS project 10192022

microPts = {'202001','202002','202006','202007','202009','202014',...
     '202016','202105','202107','202110','202114','202117','202202',...
     '202208','202210','202213','202215','202302','202306','202307',...
     '202308', '202314', '202314b', '202401', '202404', '202405',...
     '202406', '202407','202409', '202413'};

if any(contains(microPts,ptID))
    switch ptID
        case {'202001'}
            microLabels = {'left Amygdala','left Anterior Cingulate'};
            SOZ = [false, false]; % left frontal operculum LFOP 7-10
        case {'202002'}
            microLabels = {'left medial Orbital Gyrus','left Anterior Cingulate', 'empty_BankC', 'left Amygdala'};
             SOZ = [false, false, false, true]; % left anterior hippocampus LAHIP 1-4, left amygdala
        case {'202006'}
            microLabels = {'left Hippocampus','left dorsal Anterior Cingulate'};
            SOZ = [true, false]; % left anterior hippocampus LAHIP 
        case {'202007'}
            microLabels = {'left subcallosal area (vmPFC)','left Anterior Hippocampus'};
            SOZ = [false, true]; % left anterior hippocampus LAHIP 1-3
        case {'202009'}
            microLabels = {'right Gyrus Rectus','right dorsal Anterior Cingulate'};
            SOZ = [false, false]; % right hippocampus RHIP 1-3, left hippocampus.
        case {'202014'}
            microLabels = {'right OFC','right Hippocampus'};
            SOZ = [false, true]; % right amygdala & right hippocampus 1-4 (usually 1-2)
        case {'202016'}
            microLabels = {'right OFC','right Hippocampus'};
             SOZ = [false, true]; % right hippocampus 1-2, right amygdala (later seizures)
         case {'202105'}
             microLabels = {'left OFC','right Hippocampus'};
             SOZ = [false, false]; % left hippocampus 1-4
        case {'202107'}
            microLabels= {'left subgenual Cingulate','left Anterior Cingulate'};
            SOZ = [false, false]; % right hippocampus 1-4, right amygdala.
        case {'202110'}
            microLabels= {'left OFC','left subgenual Cingulate'};
            SOZ =[false, false]; % LAHIP left hippocampus 1-4
        case {'202114'}
            microLabels= {'right OFC','right Hippocampus'};
            SOZ = [false, true]; % right hippocampus 1-4 and left hippocampus, left ventral cingulate
        case {'202117'}
            microLabels= {'right OFC','right Hippocampus'};
            SOZ = [false, false]; % RPO and RCM
        case {'202202'}
            microLabels= {'left OFC','left ventral Cingulate','left dorsal Anterior Cingulate','right Anterior Hippocampus'};
            SOZ = [false, false, false, false]; % left insula and left hipp
        case {'202208'}
            microLabels= {'right OFC','right dorsal Anterior Cingulate','left Anterior Hippocampus'};
            SOZ = [true, false, false]; % right amygdala and right hipp, right ofc (subclinical).
        case {'202210'}
            microLabels= {'right OFC','right Hippocampus'};
            SOZ = [true, true]; % right hippocampus and right ofc (spread)
        case {'202213'}
            microLabels= {'right ventral Cingulate', 'left ventral Cingulate', 'left OFC'};
            SOZ = [false, false, false]; % right anterior hippocampus, right amygdala
        case {'202215'}
            microLabels= {'right OFC', 'right ventral Cingulate'};
            SOZ = [false, false]; % right amygdala 
        case {'202302'}
            microLabels= {'right OFC', 'right dorsal Anterior Cingulate', 'right Anterior Hippocampus'};
            SOZ = [false, false, true]; % right amygdala discharges. some in right hipp.
        case {'202306'}
            microLabels= {'right Anterior Cingulate', 'right Anterior Hippocampus', 'left Hippocampus'};
            SOZ = [false, false, false]; % left amygdala. 
        case {'202307'}
            microLabels= {'right Hippocampus', 'left Anterior Hippocampus', 'left Anterior Cingulate'};
            SOZ = [true, true, false]; % bilateral hippocampi and amygdala
        case {'202308'}
            microLabels= {'left OFC', 'left mid Cingulate', 'right Hippocampus'};
            SOZ = [false, false, false]; % left hippocampus.
        case {'202314'}
            microLabels = {'left OFC','left Anterior Cingulate', 'right Amygdala'};
            SOZ = [false, false, false]; % left anterior and posterior hippocampus 
        case {'202314b'}
            microLabels = {'left OFC','left Anterior Cingulate', 'right Amygdala'};
            SOZ = [false, false, false]; % No LTM report (get)
        case {'202401'}
            microLabels = {'left OFC','left ventral Cingulate', 'right Amygdala'};
            SOZ = [false, false, false]; % Left Hippocampus giving ictal discharges constantly.
        case {'202404'}
            microLabels = {'right Posterior Cingulate'};
            SOZ = false; % heschels gyrus
        case {'202405'}
            microLabels = {'right OFC', 'right  Anterior Cingulate', 'left Hippocampus'};
            SOZ = [false, false, false]; % right hippocampus and amygdala
        case {'202406'}
            microLabels = {'left OFC', 'left  Anterior Cingulate', 'right Hippocampus'};
            SOZ = [false, false, false]; % left hippocampus, left amygdala
        case {'202407'}
            microLabels = {'left OFC', 'left  Anterior Cingulate', 'left Amygdala'};
            SOZ = [false, false, false]; % right anterior hippocampus, left hippocampus
        case {'202409'}
            microLabels = {'left OFC', 'right Hippocampus','left  Anterior Cingulate'};
            SOZ = [false, true, false]; % RHIP, LHIP, RAmy, LAmy
        case {'202413'}
            microLabels = {'right OFC', 'right Anterior Cingulate','left Anterior Hippocampus'};
            SOZ = [false, true, false]; % right anterior cingulate
    end

    % actually finding the micro labels from the channel maps, and their
    % associated locations...
    try
        load(['\\155.100.91.44\D\Data\UIC' ptID '\Imaging\Registered\ChannelMap.mat'])
    catch
        load(['\\155.100.91.44\D\Data\CS' ptID '\Imaging\Registered\ChannelMap.mat'])
    end

    if ~exist('ChanMap','var')
        microChans = ChannelMap1(contains(LabelMap,'m'));
        % TODO:: allow user to pick which atlas to use. so far, just using NMM.
%         tmp = ElecAtlasProj(microChans,1);
%         % just picking the 4th electrode...
%         locIdcs = [0:length(microLabels)*2].*4;
%         microLocs = tmp(locIdcs(2:2:end-1))';
    else
        try
            microChans = ChanMap.ChannelMap1(contains(ChanMap.LabelMap,'m'));
            % TODO:: allow user to pick which atlas to use. so far, just using NMM.
            tmp = ChanMap.ElecNMMProj(microChans);
            % just picking the 4th electrode...
            locIdcs = [1:length(microLabels)]*4;
            microLocs = tmp(locIdcs)';
        catch
            microLocs = microLabels;
        end
    end
elseif strcmp(ptID,'NaP')
    fprintf('\njust returning list of patients.\n')
    microLabels = {};
    microLocs = {};
else
    fprintf('\nThis patient may not have had micros...\n')
    microLabels = {};
    microLocs = {};
end
