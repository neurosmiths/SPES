% There are some dependencies:
%   - chronux toolbox
%   - NPMK blackrock loading toolbox
%   - plot spike raster
%   - microLabelsCCEPS:
%       + custom function that provides micro labels
% Need to make it so 202314 and 202314b are treated as separate pts.
% c = uisetcolor % color picking tool!
set(0,'defaultfigurerenderer','painters') % set details for figures

%% ~~~~~~~~~~~~~~~~~~  FIRING RATE VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~ %%

% R = psthBins firing rate [trials x time]
% RtimeAll = mean of R over time. [1 x time]
% RtrialsAll = psth mean of R over trial [1 x X]
% RtrialsAllz = zscored RtrialsAll with poisson SD
% RprepostAll = R over pre-post tWin [trials x time(pre-post)]
% RprepostAllz = zscored with poisson SD RprepostAll [trials x time(pre-post)]
% RpreAllz = mean RprepostAllz over pre window [trials x 1]
% RpostAllz = mean RprepostAllz over post window [trials x 1]
% RdiffAllz = zscore mean difference btwn RpreAllz and RpostAllz
% RpreAllSD = std of mean R
% RpreAllS_3 = 3*RpreAllSD

%% ~~~~~~~~~~~~~~~~~~ define working data path ~~~~~~~~~~~~~~~~~~~~~~~~ %%

mainPath = '\\155.100.91.44\D\Data\Rhiannon\CCEPS\SPES_UNITS'; % load path for data

patientID = []; % initializing patientID

% search through main path to find data files.
nevList = subdir([mainPath '\*.nev']);
unitCount = 0;

SPES = false; % Change to False to update SPES structure (to add more units)
genPlots = false; % True to replot unit figures

%ptID = '202314b'; % uncomment for trouble shooting.

if SPES
    load("\\155.100.91.44\D\Data\Rhiannon\CCEPS\SPES_UNITS\SPESstruct.mat") % Loading SPES Structure
    % load WFprobs.mat; % loading WF metrics and neuron probabilities from CCEPunitWFs
    % logical_IN = isolatedCluster; % can change to over some probability score; FSprobs > 0.75
    %  isolatedCluster; % for small isolated cluster
    %  KMeanTwo; % load cluster information
    %  KMeanThree; % load cluster information

elseif ~SPES % Creating new SPES Structure

    % load WFprobs.mat; % loading WF metrics and neuron probabilities from CCEPunitWFs
    % logical_IN = isolatedCluster; % can change to over some probability score; FSprobs > 0.75
    % isolatedCluster; % for small isolated cluster
    % KMeanTwo; % load cluster information
    % KMeanThree; % load cluster information

    % looping over files.
    for fl = 1:length(nevList) %need to run over 22 and 23 (202314 and 202314b, to generate separate plots)
        %         % Check if the index is 17
        %         if fl == 17
        %             % Skip this iteration
        %             continue;
        %         end
        nevFile = (nevList(fl).name);
        nevFile(strfind(nevFile,'\'))='/';

        % get action potential times
        NEV = openNEV(nevFile,'nomat','nosave');
        TimeRes = NEV.MetaTags.TimeRes;
        ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode)' double(NEV.Data.Spikes.Unit)' (double(NEV.Data.Spikes.TimeStamp)./TimeRes)'];

        % get stim times from NEV dates
        [path,fname,ext] = fileparts(nevFile);

        %loading stim data
        matList =  dir(fullfile(mainPath,'CSMap',sprintf('*.mat'))); % clean CCEPS
        %  matList =  dir(fullfile(mainPath,'CSMap_All',sprintf('*.mat'))); % all Stims

        underIdx = regexp(fname,'_'); % pulling just ptID out of matlab string (take everything before _)
        patientID = fname(1:underIdx-1); % only ptID

        ptIdx = contains({matList.name},[patientID '_']);
        if ~sum(ptIdx)
            ptIdx = contains({matList.name},[patientID '_']);
        end
        load(fullfile(matList(ptIdx).folder,matList(ptIdx).name));

        %loading electrode data
        load(fullfile(mainPath,patientID,'Imaging','Registered','Electrodes.mat'));

        %loading glass brain data
        load(fullfile(mainPath,patientID,'Imaging','Registered','Surfaces.mat'));

        % get channel deets.
        inclChans = MicroElec;
        nChans = length(inclChans);

        % get microwire bundle labels
        if contains(patientID,'b')
            patientID = regexprep(patientID,'[b]',''); % removing b from ptID to get the imaging folders.
        else
        end
        [uLabels,~,SOZ] = microLabelsCCEPS(patientID(end-5:end)); % then get imaging folders

        % timing parameters.
        pre = 2; % replacing 2 to see if adaptive kernel is messing up firing rates pre-stim (didn't change anything)
        post = 3;
        post5 = 5;

        % get stim times
        trialTimes = SI30./3e4;  % stimtrials/30000
        nTrials = length(SI30);

        % align unit times
        % looping over microelectrodes
        for ch = 1:nChans
            % looping over number of units on each microelectrode
            nUnits = length(unique(ChanUnitTimestamp(inclChans(ch).*ones(size(ChanUnitTimestamp,1),1)==ChanUnitTimestamp(:,1),2)));
            for un = 1:nUnits
                fprintf('\nprocessing and plotting for unit %d of %d',un,nUnits)

                % kernel width for firing rates (negative value for adaptive measure).
                kernelWidth = -100 ./1000; % change to 100 from 50 to reduce noise.
                binWidth = 50;

                % getting unit times for the current channel and unit.
                unitTimes = ChanUnitTimestamp(ChanUnitTimestamp(:,1)==inclChans(ch) & ChanUnitTimestamp(:,2)==un,3); % in seconds

                if ~isempty(unitTimes)

                    % incrementing unit counts.
                    unitCount = unitCount+1;

                    %% ~~~~~~~~~~~~~~~~~~ cue aligned spikes ~~~~~~~~~~~~~~~~~~ %%
                    % looping over trials
                    for tt = 1:nTrials
                        % putting the data in a structure
                        spikes.channel(ch).unit(un).trial(tt).times = unitTimes(unitTimes>trialTimes(tt)-pre & unitTimes<trialTimes(tt)+post)...
                            - repmat(trialTimes(tt)-pre,length(unitTimes(unitTimes>trialTimes(tt)-pre & unitTimes<trialTimes(tt)+post)),1);
                        % and a cell for raster plot
                        spT{tt,1} = spikes.channel(ch).unit(un).trial(tt).times';
                        % calculating firing rates for each trial.
                        [~,R(tt,:),t] = psthBins(spikes.channel(ch).unit(un).trial(tt).times, binWidth, 1, 1, pre+post);
                    end % looping over trials

                    % define time frame
                    tSec = linspace(-pre,post,length(t));
                    tSec5 = linspace(-pre,post5,length(t));

                    %% ~~~~~~~~~~~~~~~~~~  specify contacts ~~~~~~~~~~~~~~~~~~ %%
                    if (inclChans(ch) - min(inclChans))<9
                        chanLabel = uLabels{1};
                        isSOZ = SOZ(1);
                    elseif (inclChans(ch) - min(inclChans))>=9 && (inclChans(ch) - min(inclChans))<17
                        chanLabel = uLabels{2};
                        isSOZ = SOZ(2);
                    elseif (inclChans(ch) - min(inclChans))>=17 && (inclChans(ch) - min(inclChans))<25
                        chanLabel = uLabels{3};
                        isSOZ = SOZ(3);
                    elseif (inclChans(ch) - min(inclChans))>=25 && (inclChans(ch) - min(inclChans))<33
                        chanLabel = uLabels{4};
                        isSOZ = SOZ(4);
                    else
                        fprintf('\n\nthere should not be units on this channel. Is there a mistake?\n')
                    end

                    % get nearby electrode stim for specific channel.
                    tmpMicroLabel = ChanLabels(MicroElec(ch));
                    macroLabel = [tmpMicroLabel{1}(2:end-1) '1']; % first macroLabel

                    % get electrode placement (gray or white matter)
                    electrodeMatter = ElecTypeRaw{ch};

                    % electrode locations.
                    trodeLocsXYZ = ElecXYZRaw; % XYZ coordinates
                    trodeLocsLabels = ElecMapRaw(:,1); % electrode Labels
                    tmpnearMicros = contains(trodeLocsLabels,macroLabel); % should this be contains?
                    nearChannels = trodeLocsLabels(tmpnearMicros,1);
                    microsToRemove = contains(nearChannels, 'm', 'IgnoreCase',false); % finding micros in nearChannels (if any remove them)
                    nearMacroLabels = nearChannels(~microsToRemove); % removing micros from macroLabels if any.
                    nearMicros = strcmp(trodeLocsLabels,nearMacroLabels);
                    NaNsToRemove = contains(trodeLocsLabels, 'NaN', 'IgnoreCase',false); % dropping NANs from electrodes.
                    trodeLocsXYZ = trodeLocsXYZ(~NaNsToRemove,:); % removing NaNs from trodeLocsXYZ
                    trodeLocsLabels = trodeLocsLabels(~NaNsToRemove); % removing NaNs from trodeLocsLabels

                    EucStr = zeros(size(SEStr)); % initialize Euclidean Dist Vector that is the same size as SEStr.

                    % euclidean distances from nearest macro ...
                    D_euclidean = sqrt(abs(sum((trodeLocsXYZ.^2) - repmat(trodeLocsXYZ(nearMicros,:),length(trodeLocsXYZ),1),2)));
                    for k = 1:length(SEStr)
                        EucStr(k) = D_euclidean(strcmp(trodeLocsLabels,SEStr{k}));
                    end

                    nearChanLogical = EucStr < 40; % Euclidean distance below 40m is our theshold
                    [sortedEuc_distance, sortedEuc_Idx] = sort(EucStr);

                    % loop to deal with empty nearChansLogical
                    if sum(nearChanLogical) == 0
                        nearMacros_EucD = {};
                    elseif sum(nearChanLogical) > 0
                        nearMacros_EucD = SEStr{nearChanLogical,1};
                    end

                    nearStims = nearChanLogical';

                    if isSOZ
                        SOZstims = contains([SEStr{:,1}],macroLabel);
                    else
                        SOZstims = false(1,length([SEStr{:,1}]));
                    end

                    % get labels for clincial SOZ
                    tmpLabels = GridLabelMap(ismember(GridChanMap,ClinicalSOZ));

                    % defining windows of interest for pre and post stim FRs
                    tWinPre = [-1.11 -0.01]; % second
                    tWinPost = [0.01 1.11]; % second
                    tWinPost3 = [0.01 2.6]; % 2.5 second window.
                    tWinPost25 = [0.01 2.51];
                    tWinPreStimFR = [-0.1 0.0]; % post stim. 100ms - 0ms.
                    tWinPostStimFR = [0.0 0.1]; % post stim. 0ms - 100ms.

                    % [========== Saving to SPES structure ==========]

                    SPESstruct(unitCount).patientID = patientID;
                    SPESstruct(unitCount).chanLabel = chanLabel;
                    SPESstruct(unitCount).electrodeMatter = electrodeMatter;
                    %     SPESstruct(unitCount).logical_IN = logical_IN(unitCount);
                    %     SPESstruct(unitCount).KMeanTwo = KMeanTwo;
                    %    SPESstruct(unitCount).KMeanThree = KMeanThree;
                    SPESstruct(unitCount).unitTimes = unitTimes; %(*30000)
                    SPESstruct(unitCount).stimTimes = trialTimes;
                    SPESstruct(unitCount).sortEucIdx = sortedEuc_Idx; % distance length of stim 
                    SPESstruct(unitCount).sortedEuc_distance = sortedEuc_distance;
                    SPESstruct(unitCount).tmpMicroLabel = tmpMicroLabel;
                    SPESstruct(unitCount).macroLabel = macroLabel;
                    SPESstruct(unitCount).isSOZ = isSOZ;
                    SPESstruct(unitCount).SOZstims = SOZstims;
                    SPESstruct(unitCount).nearStims = nearStims;
                    SPESstruct(unitCount).D_euclidean = D_euclidean;
                    SPESstruct(unitCount).nearMacros_EucD = nearMacros_EucD;

                    %% ~~~~~~~~~~~~~~~~~~~ Mean Firing Rates ~~~~~~~~~~~~~~~~~~~~ %%
                    % mean firing rates and SEs for all/SOZ trials
                    RtimeAll = mean(R); % mean FR over all trials for each time point
                    EtimeAll = std(R(nTrials,:))./sqrt(sum(nTrials));
                    % Near
                    if ~isempty(nearStims) % sum(RpreNearz) > 0
                        RtimeNear = mean(R(nearStims,:)); % mean FR over near trials for each time point
                        EtimeNear = std(R(nearStims,:))./sqrt(sum(nearStims));
                    elseif isempty(nearStims) % um(RpreNearz) < 0
                        RtimeNear = NaN;
                        EtimeNear = NaN;
                    end
                    % Far
                    if isempty(nearStims)
                        RtimeFar = mean(R);
                        EtimeFar = std(R(nTrials,:))./sqrt(sum(nTrials));
                    elseif ~isempty(nearStims)
                        RtimeFar = mean(R(~nearStims,:)); % mean FR over far trials for each time point
                        EtimeFar = std(R(~nearStims,:))./sqrt(sum(~nearStims));
                    end

                    % [========== Saving to SPES structure ==========]

                    SPESstruct(unitCount).R = R;
                    SPESstruct(unitCount).spikeWaveforms = NEV.Data.Spikes.Waveform(:,ChanUnitTimestamp(:,1)==inclChans(ch) & ChanUnitTimestamp(:,2)==un)./4;
                    SPESstruct(unitCount).RtimeAll = RtimeAll;
                    SPESstruct(unitCount).RtimeNear = RtimeNear;
                    SPESstruct(unitCount).RtimeFar = RtimeFar;

                    %% ~~~~~~~~~~~~~~ Mean Firing Rates across trials using psth chronux ~~~~~~~~~~~~ %%

                    %kernelWidth_Euc = 200 ./1000;
                    [RtrialsAll,TtrialsAll,EtrialsAll] = psth(spikes.channel(ch).unit(un).trial, kernelWidth, 'n', [0 pre+post]); % ALL
                    tsecA = TtrialsAll-repmat(pre,1,length(TtrialsAll));
                    [RtrialsFar,TtrialsFar,EtrialsFar] = psth(spikes.channel(ch).unit(un).trial(~nearStims), kernelWidth, 'n', [0 pre+post]); % FAR
                    tsecF = TtrialsFar-repmat(pre,1,length(TtrialsFar));
                    if  sum(nearStims) > 0
                        [RtrialsNear,TtrialsNear,EtrialsNear] = psth(spikes.channel(ch).unit(un).trial(nearStims), kernelWidth, 'n', [0 pre+post]); % NEAR
                        tsecN = TtrialsNear-repmat(pre,1,length(TtrialsNear));
                        if isempty(RtrialsNear)
                            RtrialsNear = NaN;
                            TtrialsNear = NaN;
                            EtrialsNear = NaN;
                        elseif ~isempty(RtrialsNear)
                        end
                    elseif sum(nearStims) == 0
                        RtrialsNear = NaN;
                        TtrialsNear = NaN;
                        EtrialsNear = NaN;
                    end

                    % [========== Saving to SPES structure ==========]

                    SPESstruct(unitCount).RtrialsAll = RtrialsAll;
                    SPESstruct(unitCount).TtrialsAll = TtrialsAll;
                    SPESstruct(unitCount).EtrialsAll = EtrialsAll;
                    SPESstruct(unitCount).RtrialsFar = RtrialsFar;
                    SPESstruct(unitCount).TtrialsFar = TtrialsFar;
                    SPESstruct(unitCount).EtrialsFar = EtrialsFar;
                    SPESstruct(unitCount).RtrialsNear = RtrialsNear;
                    SPESstruct(unitCount).TtrialsNear = TtrialsNear;
                    SPESstruct(unitCount).EtrialsNear = EtrialsNear;

                    %% ~~~~~~~~~~~~~~ Mean Firing Rates across stimulations using psth chronux ~~~~~~~~~~~~ %%

                    % Convert to a cell array of character vectors
                    SEStr_char = cellfun(@(x) char(x), SEStr, 'UniformOutput', false);
                    SEStr_string = string(SEStr_char);
                    [uniqueOccurrences,ida,idc] = unique(SEStr_string); % [unique macros, idc: count of unique macros, idx:
                    a_counts = accumarray(idc,1);
                    value_counts = [uniqueOccurrences, a_counts]; % string array
                    value_counts = cellstr(value_counts); % cell array

                    [~,TstimAll,~] = psth(spikes.channel(ch).unit(un).trial, kernelWidth, 'n', [0 pre+post]);
                    % Iterate over unique SEStr entries
                    for i = length(uniqueOccurrences):-1:1
                        % Find indices where the current unique entry occurs in SEStr
                        idx = find(strcmp(SEStr_char, uniqueOccurrences{i}));
                        EucD_unique(i) = unique(EucStr(idc == i));
                        % Extract spike times corresponding to the current unique entry
                        spikeTimes{i} = spT(idx);
                        try
                            [RstimAll(i,:),TstimAll,EstimAll(i,:)] = psth(spikes.channel(ch).unit(un).trial(idx), kernelWidth, 'n', [0 pre+post]);

                        catch
                            RstimAll(i,:) = zeros(size(TstimAll));
                            EstimAll(i,:) = zeros(size(TstimAll));

                        end
                    end

                    % sort by Euclidean distance
                    [ascendingOrder,sortIdcs] = sort(EucD_unique);

                    plotEuc = true;

                    sortIdcs_local = sortIdcs(ascendingOrder < 15); % near stim sites: less than 15mm
                    sortIdcs_distant = sortIdcs(ascendingOrder > 15); % far stim sites: more than 15mm
                    sortIdcs_below = sortIdcs(ascendingOrder <= 39);
                    sortIdcs_above = sortIdcs(ascendingOrder > 39);
                                       
                    if plotEuc
                        % plot
                        figure(15)
                        subplot(4,3,[1:2,4:5])
                        imagesc(tSec,ascendingOrder,RstimAll(sortIdcs,:),[0 (3*max(max(RstimAll))/4)])
                        axis xy
                        xlabel('time')
                        ylabel('electrode by distance')
                        xlim([-1.5 2.5])

                        colorbar(gca, 'southoutside')
                        % make a vector
                        EucFR_vector = max(RstimAll(sortIdcs,tsecA>0.05 & tsecA<0.3),[],2);
                        % use find change points
                        [FRdistance_changepts, residual_changepts] = findchangepts(EucFR_vector,Statistic="mean"); % 'MaxNumChanges',20
                        initial_params = [max(EucFR_vector),0, NaN, NaN] % max, min, midpoint, slope
                        [param, stat] = sigm_fit(ascendingOrder',EucFR_vector,[],[],0)
                        r = EucFR_vector-stat.ypred;
                        SSE = norm(r,2)^2;
                        rsquare = 1-SSE/norm(EucFR_vector-mean(EucFR_vector)).^2
                        title(sprintf('FR by distance. ResErr: %0.2f, x50: %0.2f, rsq: %0.2f', residual_changepts, param(:,3), rsquare))
                        % use wavelet
                        % FRdistance_changepts = modwt(EucFR_vector,'haar',1)
                        if ~isempty(FRdistance_changepts)
                            line([-1.5 2.5], [ascendingOrder(FRdistance_changepts) ascendingOrder(FRdistance_changepts)],'color','w','linewidth',0.75,'linestyle','--')
                        elseif isempty(FRdistance_changepts)
                        end
                        % add in plot with Euc (y) and FR (x) to clearly show change point line
                        subplot(4,3,[3,6])
                        plot(EucFR_vector,ascendingOrder)
                        xlabel('FR')
                        ylabel('Euc distance')
                        colorbar(gca, 'southoutside')
                        axis tight
                        % Below 30mm euclidean distance FRs
                        subplot(4,3,7:8)
                        hold on
                        line([0 0],[0 5],'color','w','linewidth',0.75,'linestyle','--')
                        patch([tsecA fliplr(tsecA)],[mean(RstimAll(ascendingOrder <= param(3),:),1)+mean(EstimAll(ascendingOrder <= param(3),:),1) fliplr(mean(RstimAll(ascendingOrder <= param(3),:),1)-mean(EstimAll(ascendingOrder <= param(3),:),1))],[0.7 0.7 0.7],'facealpha',0.5,'edgecolor','none')
                        plot(tsecA,mean(RstimAll(ascendingOrder <= param(3),:)),'color', [1 1 0]) % yellow for evoked firing....
                        patch([tsecA fliplr(tsecA)],[mean(RstimAll(ascendingOrder > param(3),:),1)+mean(EstimAll(ascendingOrder > param(3),:),1) fliplr(mean(RstimAll(ascendingOrder > param(3),:),1)-mean(EstimAll(ascendingOrder > param(3),:),1))],[0.7 0.7 0.7],'facealpha',0.5,'edgecolor','none')
                        plot(tsecA,mean(RstimAll(sortIdcs(ascendingOrder > param(3)),:)),'color', 'blue') % blue for less evoked firing...
                        % text(1.5,3, sprintf('sigmod: %0.2f', param(3),'Color','red'))
                        hold off
                        xlim([-1.5 2.5])
                        xlabel('time (s)')
                        ylabel('firing rate (Hz)')
                        set(gca,'Color','k')
                        % regression model
                        subplot(4,3,10)
                        if sum(RstimAll) == 0
                            text(0,0,'RStimAll is empty')
                        else
                            euclidean_lm = fitlm(EucD_unique,mean(RstimAll(:,tSec>tWinPost(1) & tSec < tWinPost(2)),2), "Exclude", mean(RstimAll(:,tSec>tWinPost(1) & tSec < tWinPost(2)),2) == 0); % one second time window post-stim
                            plot(euclidean_lm)
                            title(sprintf('micro: %s\npVal= %0.5f, Fstat= %0.5f', chanLabel,euclidean_lm.ModelFitVsNullModel.Pvalue,euclidean_lm.ModelFitVsNullModel.Fstat))
                            xlabel('distance')
                            ylabel('mean FR 1second post')
                            axis square
                            subplot(4,3,11)
                            euclidean3_lm = fitlm(EucD_unique,mean(RstimAll(:,tSec>tWinPost(1) & tSec < tWinPost3(2)),2), "Exclude", mean(RstimAll(:,tSec>tWinPost(1) & tSec < tWinPost3(2)),2) == 0); % three second time window post-stim
                            plot(euclidean3_lm)
                            title(sprintf('micro: %s\npVal= %0.5f, Fstat= %0.5f', chanLabel,euclidean3_lm.ModelFitVsNullModel.Pvalue,euclidean3_lm.ModelFitVsNullModel.Fstat))
                            xlabel('distance')
                            ylabel('mean FR (3seconds post)')
                            axis square
                        end
                        halfMaximize(figure(15), 'page')

                        % saving figure
                        saveas(15,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\EuclideanDistance\',sprintf('EuclideanDistance_FRchanges_%s_ch%d_un%d.pdf',patientID,ch,un)))
                        close(15)
                        
                        % plot (GLASS BRAIN) 
                        figure(16)
                        subplot(3,2,1)
                        imagesc(tSec,ascendingOrder,RstimAll(sortIdcs,:),[0 (3*max(max(RstimAll))/4)])
                        axis xy
                        xlabel('time')
                        ylabel('electrode by distance')
                        xlim([-1.5 2.5])
                        colorbar(gca, 'southoutside')
                        % make a vector
                        EucFR_vector = max(RstimAll(sortIdcs,tsecA>0.05 & tsecA<0.3),[],2);
                        % use find change points
                        [FRdistance_changepts, residual_changepts] = findchangepts(EucFR_vector,Statistic="mean"); % 'MaxNumChanges',20
                        initial_params = [max(EucFR_vector),0, NaN, NaN] % max, min, midpoint, slope
                        [param, stat] = sigm_fit(ascendingOrder',EucFR_vector,[],[],0)
                        r = EucFR_vector-stat.ypred;
                        SSE = norm(r,2)^2;
                        rsquare = 1-SSE/norm(EucFR_vector-mean(EucFR_vector)).^2
                        title(sprintf('FR by distance. ResErr: %0.2f, x50: %0.2f, rsq: %0.2f', residual_changepts, param(:,3), rsquare))
                        if ~isempty(FRdistance_changepts)
                            line([-1.5 2.5], [ascendingOrder(FRdistance_changepts) ascendingOrder(FRdistance_changepts)],'color','w','linewidth',0.75,'linestyle','--')
                        elseif isempty(FRdistance_changepts)
                        end
                        % Below 30mm euclidean distance FRs
                        subplot(3,2,4)
                        hold on
                        line([0 0],[0 5],'color','w','linewidth',0.75,'linestyle','--')
                        patch([tsecA fliplr(tsecA)],[mean(RstimAll(ascendingOrder <= param(3),:),1)+mean(EstimAll(ascendingOrder <= param(3),:),1) fliplr(mean(RstimAll(ascendingOrder <= param(3),:),1)-mean(EstimAll(ascendingOrder <= param(3),:),1))],[0.7 0.7 0.7],'facealpha',0.5,'edgecolor','none')
                        plot(tsecA,mean(RstimAll(ascendingOrder <= param(3),:)),'color', [1 1 0]) % yellow for evoked firing....
                        patch([tsecA fliplr(tsecA)],[mean(RstimAll(ascendingOrder > param(3),:),1)+mean(EstimAll(ascendingOrder > param(3),:),1) fliplr(mean(RstimAll(ascendingOrder > param(3),:),1)-mean(EstimAll(ascendingOrder > param(3),:),1))],[0.7 0.7 0.7],'facealpha',0.5,'edgecolor','none')
                        plot(tsecA,mean(RstimAll(sortIdcs(ascendingOrder > param(3)),:)),'color', 'blue') % blue for less evoked firing...
                        % text(1.5,3, sprintf('sigmod: %0.2f', param(3),'Color','red'))
                        hold off
                        xlim([-1.5 2.5])
                        xlabel('time (s)')
                        ylabel('firing rate (Hz)')
                        set(gca,'Color','k')
                        subplot(3,2,[2:3,5:6])
                        
                      %   % glass brain figure
                      %   figure(5)
                      %   hold on
                      %   s1 =  scatter3(ElecXYZRaw(D_euclidean>param(3),1), ElecXYZRaw(D_euclidean>param(3),2), ElecXYZRaw(D_euclidean>param(3),3), 30,'b', 'filled','MarkerEdgeColor','k');
                      %   s2 = scatter3(ElecXYZRaw(D_euclidean<=param(3),1), ElecXYZRaw(D_euclidean<=param(3),2), ElecXYZRaw(D_euclidean<=param(3),3), 30, 'yellow', 'filled','MarkerEdgeColor','k');
                      %   s3 = scatter3(ElecXYZRaw(MicroElecRaw(un:un+7),1), ElecXYZRaw(MicroElecRaw(un:un+7),2), ElecXYZRaw(MicroElecRaw(un:un+7),3), 30, 'MarkerEdgeColor','k','MarkerFaceColor','r');
                      %   p = patch('Faces', BrainSurfRaw.faces, 'Vertices', BrainSurfRaw.vertices, 'edgecolor', 'none', 'facecolor', 'flat', 'facealpha', .2);
                      %   facecolor = repmat([1 1 1 ],length(BrainSurfRaw.faces),1);
                      %   set(p, 'FaceVertexCData', facecolor);
                      %   camlight
                      %   lighting gouraud
                      % %  view(3)
                      %   halfMaximize(figure(5), 'page')
                      %   % saving figure
                      %   saveas(5,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\Glass_Brain\',sprintf('GlassBrainExampleAxial2_%s_ch%d_un%d.pdf',patientID,ch,un)))

                        % sigmodial plot
                        figure(456)
                        sigm_fit(ascendingOrder',EucFR_vector,[],[],1)
                        % saving figure
                        saveas(456,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\EuclideanDistance\Sigmodial\',sprintf('Sigmodial_%s_ch%d_un%d.pdf',patientID,ch,un)))
                        close(456)

                        % what are the mean FRs for each time bin? For stimlations that evoked near firing, how long did that firing rate last for? (0-1000ms?)
                        FR_below = mean(RstimAll(sortIdcs_below,tsecA>0 & tsecA<1)); % mean firing rate for one second post stimulation.
                        meanFR_belowBL = mean(mean(RstimAll(sortIdcs_below,tsecA>-1 & tsecA<0))); % calculate mean baseline firing one second post stim.
                        logical_FRbelow = FR_below >= meanFR_belowBL; % post stim firing rate
                        timeWindow_FRbelow = tsecA(tsecA>0 & tsecA<1); % get time window of times that FR occured until....
                        duration_FRBelow = max(timeWindow_FRbelow(logical_FRbelow)); % find max time that FR above the baseline threshold lasts for!
                        meanFR_belowEvoked = mean(FR_below(logical_FRbelow));

                        % Find evoked above distances too
                        FR_above = mean(RstimAll(sortIdcs_above,tsecA>0 & tsecA<1)); % mean firing rate for one second post stimulation.
                        meanFR_aboveBL = mean(mean(RstimAll(sortIdcs_above,tsecA>-1 & tsecA<0))); % calculate mean baseline firing one second post stim.
                        logical_FRabove = FR_above >= meanFR_aboveBL; % post stim firing rate
                        timeWindow_FRabove = tsecA(tsecA>0 & tsecA<1); % get time window of times that FR occured until....
                        duration_FRAbove = max(timeWindow_FRabove(logical_FRabove)); % find max time that FR above the baseline threshold lasts for!
                        meanFR_aboveEvoked = mean(FR_above(logical_FRabove));

                        % What regions of the brain do the longer evoked FRs come from?

                        %EucFR_model = fitglme(EucFR_table,'FR ~ Euclidean_Distance','link','logit','Distribution','InverseGaussian')

                        SPESstruct(unitCount).euclidean_lm = euclidean_lm;
                        SPESstruct(unitCount).euclidean3_lm = euclidean3_lm;
                        SPESstruct(unitCount).euclidean_lm_PValue = euclidean_lm.ModelFitVsNullModel.Pvalue;
                        SPESstruct(unitCount).euclidean3_lm_PValue = euclidean3_lm.ModelFitVsNullModel.Pvalue;
                        SPESstruct(unitCount).duration_FRBelow = duration_FRBelow;
                        SPESstruct(unitCount).meanFR_belowEvoked = meanFR_belowEvoked;
                        SPESstruct(unitCount).FRdistance_changepts = FRdistance_changepts;
                        SPESstruct(unitCount).residual_changepts =  residual_changepts;
                        SPESstruct(unitCount).rsquare = rsquare;
                        SPESstruct(unitCount).sigmod_distance = param(:,3);
                        SPESstruct(unitCount).sigmod_slope = param(:,4);
                        SPESstruct(unitCount).duration_FRAbove = duration_FRAbove;
                        SPESstruct(unitCount).meanFR_aboveEvoked = meanFR_aboveEvoked;
                        SPESstruct(unitCount).ascendingOrder = ascendingOrder;
                        SPESstruct(unitCount).EucFR_vector = EucFR_vector;
                        SPESstruct(unitCount).sortIdcs_below = sortIdcs_below;
                        SPESstruct(unitCount).sortIdcs_above = sortIdcs_above;
                        SPESstruct(unitCount).sortIdcs_above = sortIdcs_above;
                        SPESstruct(unitCount).sortIdcs = sortIdcs;

                        clear ascendingOrder sortIdcs sortIdcs_below sortIdcs_above RstimAll EucD_unique

                    elseif ~plotEuc
                    end

                    %% ~~~~~~~~~~~~~~ ANALYSIS 1: IS THE UNIT BEING MODULATED? ~~~~~~~~~~~~~~~~~ %

                    % Calculate: mean z-scored firing rates, with a poisson stand deviation, between pre and post time windows.
                    % Outcomes: (a) Increased FR (amplification), (b) Decreased FR(suppression), (c) no modulation

                    % Step 1: Create a trials (rows) x time (columns)
                    % matrix for firing rates over whole time window (psthBins)
                    RprepostAll = (R(:,tSec>tWinPre(1) & tSec<tWinPost3(2))); % using post3 to cast a wide net.

                    poisson = false;
                    if poisson
                        % Step 2: poisson fit for data to calculate lambda and poisson standard deviation
                        [lambdahatAll, lambdaciAll] = poissfit(histcounts(mean(RprepostAll))); % ALL
                        [~,varianceAll] = poisstat(lambdahatAll);
                        sdAll = sqrt(varianceAll);

                    elseif ~poisson

                        sdAll = std(mean(RprepostAll));

                    end

                    % Step 4: Create z-scored psthBins
                    RtrialsAllz = (RprepostAll - mean(RprepostAll)/sdAll);

                    % Step 5: Z-Scored over pre and post time windows
                    RpreAllz = RtrialsAllz(tSec>tWinPre(1) & tSec<tWinPre(2));
                    RpostAllz = RtrialsAllz(tSec>tWinPost(1) & tSec<tWinPost(2));
                  
                    % Step 6: Differences in firings rates over trials and time.
                    RdiffAllz = mean(RpreAllz) - mean(RpostAllz); % difference between zscore pre and post
                  
                    % Step 7 : Statistically test difference between zScored Pre and Post FRs averaged for each trial.
                    [hAllz,pAllz,ciAllz, statsAllz] = ttest(RpreAllz,RpostAllz); % ttest over trials.
                  
                    % Computing FR change without zscore.
                    RpreAll = RtrialsAll(tSec>tWinPre(1) & tSec<tWinPre(2));
                    RpostAll = RtrialsAll(tSec>tWinPost(1) & tSec<tWinPost(2));
                    RdiffAll = mean(RpreAll) - mean(RpostAll);
                    [hAll,pAll,ciAll, statsAll] = ttest(RpreAll,RpostAll); % ttest over trials.

                    % [========== Saving to SPES structure ==========]

                    SPESstruct(unitCount).RtrialsAllz = RtrialsAllz;
                    SPESstruct(unitCount).RprepostAll = RprepostAll;
                    SPESstruct(unitCount).sdAll = sdAll;
                    SPESstruct(unitCount).RpreAllz = RpreAllz;
                    SPESstruct(unitCount).RpostAllz = RpostAllz;               
                    SPESstruct(unitCount).RdiffAllz = RdiffAllz;                  
                    SPESstruct(unitCount).pAllz = pAllz;
                    SPESstruct(unitCount).hAllz = hAllz;
                    SPESstruct(unitCount).ciAllz = ciAllz;
                    SPESstruct(unitCount).statsAllz = statsAllz;                    
                    SPESstruct(unitCount).RpreAll = RpreAll;
                    SPESstruct(unitCount).RpostAll = RpostAll;
                    SPESstruct(unitCount).pAll = pAll;
                    SPESstruct(unitCount).hAll = hAll;
                    SPESstruct(unitCount).ciAll = ciAll;
                    SPESstruct(unitCount).statsAll = statsAll;
                    SPESstruct(unitCount).RdiffAll = RdiffAll;

                    %% ///////////////////////////////////// TRIAL AVERAGED FIRING RATES: UNIT CHARACTERISTICS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ %%

                    % SECTION: ~~~~~~~ TRIAL-AVERAGE FIRING RATES ~~~~~~~ %
                    %% ~~~~~~~~~ (1a) BASELINE FREQUENCY (Hz):The average of the instantaneous firing rate for the 1s preceding stimulation  ~~~~~~~~~~~~~~ %

                    BLfreqAll = mean(RtrialsAll(tsecA>tWinPre(1) & tsecA<tWinPre(2))); % Average FR
                    BLfreqAllz = mean(RtrialsAllz(tsecA>tWinPre(1) & tsecA<tWinPre(2))); % Z-Scored FR

                    % [========== Saving to SPES structure ==========]
                    SPESstruct(unitCount).BLfreqAll = BLfreqAll;
                    SPESstruct(unitCount).BLfreqAllz = BLfreqAllz;

                    %% ~~~~~~~~~ (1b) POSTSTIM FREQUENCY (Hz):The average of the instantaneous firing rate for the 1s after stimulation  ~~~~~~~~~~~~~~ %

                    PSfreqAll = mean(RtrialsAll(tsecA>tWinPost(1) & tsecA<tWinPost(2))); % Average FR
                    PSfreqAllz = mean(RtrialsAllz(tsecA>tWinPost(1) & tsecA<tWinPost(2))); % Z-Scored FR
                    [Maximum_FRAll, I_Maximum_FRAll] = max(RtrialsAll(tsecA>tWinPost(1) & tsecA<tWinPost(2))); % Max firing rate across post time window

                    % [========== Saving to SPES structure ==========]
                    SPESstruct(unitCount).PSfreqAll = PSfreqAll;
                    SPESstruct(unitCount).PSfreqAllz = PSfreqAllz;
                    SPESstruct(unitCount).Maximum_FRAll = Maximum_FRAll;
                    SPESstruct(unitCount).I_Maximum_FRAll = I_Maximum_FRAll;

                    %% ~~~~~~~~~ (2) AMPLITUDE OF THE SUPPRESSION: The lowest instantaneous firing rate during suppression. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
                    %  Calculated by finding the minimum FR in the postWindow (0.01-1.11) which should contain the lowest suppression amplitude.
                    
                    [Minimum_SuppressionAll, I_Minimum_SuppressionAll] = min(RtrialsAll(tsecA>tWinPost(1) & tsecA<tWinPost25(2))); % average firing rates across post time window, then find minimum FR
                    [Minimum_SuppressionAllz, I_Minimum_SuppressionAllz] = min(RtrialsAllz(tsecA>tWinPost(1) & tsecA<tWinPost25(2))); % z-scored  average firing rates across post time window, then find minimum FR

                    % Finding minimum Amp for stims below 40mm.... 
                    [Minimum_SuppressionAll_Near, I_Minimum_SuppressionAll_Near] = min(RtrialsNear(tsecA>tWinPost(1) & tsecA<tWinPost25(2))); % average firing rates across post time window, then find minimum FR
                    [Minimum_SuppressionAll_Far, I_Minimum_SuppressionAll_Far] = min(RtrialsFar(tsecA>tWinPost(1) & tsecA<tWinPost25(2))); % average firing rates across post time window, then find minimum FR

                    % [========== Saving to SPES structure ==========]
                    SPESstruct(unitCount).Minimum_SuppressionAll = Minimum_SuppressionAll;
                    SPESstruct(unitCount).I_Minimum_SuppressionAll = I_Minimum_SuppressionAll;
                    SPESstruct(unitCount).Minimum_SuppressionAllz = Minimum_SuppressionAllz;
                    SPESstruct(unitCount).I_Minimum_SuppressionAllz = I_Minimum_SuppressionAllz;
                    SPESstruct(unitCount).Minimum_SuppressionAll_Near = Minimum_SuppressionAll_Near;
                    SPESstruct(unitCount).I_Minimum_SuppressionAll_Near = I_Minimum_SuppressionAll_Near;
                    SPESstruct(unitCount).Minimum_SuppressionAll_Far = Minimum_SuppressionAll_Far;
                    SPESstruct(unitCount).I_Minimum_SuppressionAll_Far = I_Minimum_SuppressionAll_Far;

                    %% ~~~~~~~~~ (3) LATENCY TO REACH SUPPRESSION AMPLITUDE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
                    % This is partially calculated above in (2) by getting the minimum FR and then here, finding the timepoint the minFR occurs at, after stim.

                    SuppLatencyAll = tsecA(min(find(tsecA>tWinPost(1) & tsecA<tWinPost25(2))+I_Minimum_SuppressionAll)); % average across all firing rates: time to minimum suppression
                    SuppLatencyAllz = tsecA(min(find(tsecA>tWinPost(1) & tsecA<tWinPost25(2))+I_Minimum_SuppressionAllz)); % z-scored firing rates:  time to minimum suppression
                    SuppLatencyAll_Near = tsecA(min(find(tsecA>tWinPost(1) & tsecA<tWinPost25(2))+I_Minimum_SuppressionAll_Near)); % average across all firing rates: time to minimum suppression
                    SuppLatencyAll_Far = tsecA(min(find(tsecA>tWinPost(1) & tsecA<tWinPost25(2))+I_Minimum_SuppressionAll_Far)); % average across all firing rates: time to minimum suppression

                    % [========== Saving to SPES structure ==========]
                    SPESstruct(unitCount).suppLatencyAll = SuppLatencyAll;
                    SPESstruct(unitCount).SuppLatencyAllz = SuppLatencyAllz;
                    SPESstruct(unitCount).suppLatency_Near = SuppLatencyAll_Near;
                    SPESstruct(unitCount).suppLatency_Far = SuppLatencyAll_Far;

                    %% ~~~~~~~~~ (4) DURATION OF SUPPRESSION: Difference in time since the instantaneous firing rate crossed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
                    %  downwards and upwards the 2 SD lower threshold below the baseline frequency.
                    RpreAllSD = std(RtrialsAll(tsecA>tWinPre(1) & tsecA<tWinPre(2))); % Step 1: calculate one STD of mean firing rate over time
                    RpreAllSD_2 = RpreAllSD*2; % Step 2: calculate two STD
                    RthresholdAll = BLfreqAll-RpreAllSD_2; % Step 3: Caluculating FR suppression threshold (finding threshold for upper limit of suppression (start of suppression FR)).
                    suppressionWin = tsecA(tsecA>tWinPost25(1) & tsecA<tWinPost25(2)); % Step 4: Create a suppression window to examine FR

                    % Step 5a: Duration of Suppression for ALL
                    suppressionAllLogical = RtrialsAll(tsecA>tWinPost25(1) & tsecA<tWinPost25(2)) <= RthresholdAll; % mean FR btwn postwin less than threshold
                    suppressionDurationAll = suppressionWin(suppressionAllLogical); % index by window to find when supression happens
                    suppressionAllIdx = find(~suppressionAllLogical,1,'first'); % find first 0 in suppression array
                    %suppressionDurationAll(suppressionAllIdx:end) = []; %
                    %delete everything after the first zero (after suppression)
                    suppressionDurationTimeAll = max(suppressionDurationAll) - min(suppressionDurationAll); % Duration of Suppression in Time

                    suppressionDurationTimeAllz = zscore(suppressionDurationTimeAll);
                    RthresholdAllz = zscore(RthresholdAll);

                     %% ~~~~~~~~~ (4) RATE OF RECOVERY: Difference in time between minimum suppression an upper threshold ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
                    suppressionMinimumLogical = RtrialsAll(tsecA>tWinPost25(1) & tsecA<tWinPost25(2)) == Minimum_SuppressionAll; % find logical index of point of minimum suppression
                    recoveryPoint_tmp = find(suppressionAllLogical,1,'last'); % threshold point.. calculate the time when the threshold is crossed for the end of suppression... 
                    recoveryPointLogical = false(1, length(suppressionAllLogical)); % Create a new logical array of the same size
                    recoveryPointLogical(recoveryPoint_tmp) = true; % logical index for end time.
                    recoveryRate_Time = suppressionWin(recoveryPointLogical) - suppressionWin(suppressionMinimumLogical);
                    recoveryRate_Timez = zscore(recoveryRate_Time);

                    % [========== Saving to SPES structure ==========]
                    SPESstruct(unitCount).RpreAllSD = RpreAllSD;
                    SPESstruct(unitCount).RpreAllSD_3 = RpreAllSD_2;
                    SPESstruct(unitCount).durationThresholdAll = RthresholdAll;
                    SPESstruct(unitCount).suppressionAllLogical = suppressionAllLogical;
                    SPESstruct(unitCount).suppressionDurationTimeAll = suppressionDurationTimeAll;
                    SPESstruct(unitCount).suppressionDurationTimeAllz = suppressionDurationTimeAllz;
                    SPESstruct(unitCount).durationThresholdAllz = RthresholdAllz;
                    
                    % recovery rate saving
                    SPESstruct(unitCount).suppressionMinimumLogical = suppressionMinimumLogical;
                    SPESstruct(unitCount).recoveryPointLogical = recoveryPointLogical;
                    SPESstruct(unitCount).recoveryRate_Time = recoveryRate_Time;
                    SPESstruct(unitCount).recoveryRate_Timez = recoveryRate_Timez;

                    %% Area Under the Curve - threshold to threshold %%
                    % calculate the area between the two points on firing
                    % rate where the threshold was crossed.
                  
                    tmpFR = -RtrialsAll(tsecA>tWinPost25(1) & tsecA<tWinPost25(2)) + RthresholdAll; % firing rates below than suppression threshold
                    logical_AUC = tmpFR > 0; % logical vector for FR below threshold
                    tmpT = tsecA(tsecA>tWinPost25(1) & tsecA<tWinPost25(2)); % create time vector to plot AUC
                    if sum(logical_AUC) > 1
                        AUC = trapz(tmpT(logical_AUC), tmpFR(logical_AUC)); % calculating the Area Under the Curve (suppression)

                        figure(111) % plotting suppression.
                        hold on
                        area(tmpT(logical_AUC),tmpFR(logical_AUC))
                        text(0.2, 2,sprintf('Area = %.2f', AUC))
                        hold off
                        saveas(111,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AUC\',sprintf('suppressionAUC_%s_ch%d_un%d.pdf',patientID,ch,un)))
                        close(111)

                    else
                    end

                    SPESstruct(unitCount).AUC = AUC; % if AUC is 0 then there was no suppression. # of AUC should match number of sig suppressed units
                    SPESstruct(unitCount).tmpFR = tmpFR;
                    SPESstruct(unitCount).tmpT = tmpT;
                    SPESstruct(unitCount).logical_AUC = logical_AUC;

                  
                    %% Exponential Decay Function
%                     
%                     % Define the exponential decay function
%                     % R(t) = R_baseline - (R_baseline - R_min) * exp(-t / tau)
%                     expDecayModel = @(params, t) params(1) - (params(1) - params(2)) .* exp(-t / params(3));
% 
%                     %params(1) is the baseline firing rate.
%                     %params(2) is the minimum firing rate during suppression.
%                     %params(3) is the time constant 𝜏 τ, which determines the rate of suppression and recovery.
% 
%                     % Initial guesses for the parameters: [R_baseline, R_min, tau]
%                     initialParams = [BLfreqAll, Minimum_SuppressionAll, 1];
% 
%                     % Fit the model to the data
%                     options = optimoptions('lsqcurvefit', 'Display', 'off');
%                     [estimatedParams, ~] = lsqcurvefit(@(params, t) expDecayModel(params, t), initialParams, tsecA>tWinPost(1) & tsecA<tWinPost3(2), RtrialsAll, [], []);
% 
%                     % Plot original data
%                     figure(888);
%                     plot(tsecA>tWinPost(1) & tsecA<tWinPost3(2), RtrialsAll(tsecA>tWinPost(1) & tsecA<tWinPost3(2)), 'bo', 'MarkerFaceColor', 'b'); % Original data
%                     % Plot fitted curve
%                     fittedFiringRate = expDecayModel(estimatedParams, tsecA>tWinPost(1) & tsecA<tWinPost3(2));
%                     hold on;
%                     plot(tsecA>tWinPost(1) & tsecA<tWinPost3(2), fittedFiringRate, 'r-', 'LineWidth', 2); % Fitted curve
%                     % Add labels and legend
%                     xlabel('Time (s)');
%                     ylabel('Firing Rate (Hz)');
%                     legend('Data', 'Fitted Curve');
%                     title('Exponential Decay Fit to Firing Rate Data');
%                     maximize(888)
% 
%                     % Display estimated parameters
%                     disp('Estimated Parameters:');
%                     disp(['Baseline Firing Rate (R_baseline): ', num2str(estimatedParams(1))]);
%                     disp(['Minimum Firing Rate (R_min): ', num2str(estimatedParams(2))]);
%                     disp(['Time Constant (tau): ', num2str(estimatedParams(3))]);
% 
%                     %Time Constant (tau): Determines how quickly the firing rate decreases and recovers. A smaller ?
%                     %τ indicates a faster suppression and recovery, while a larger
%                     %τ indicates a slower process.


                    %% PLOTTING FIRING RATE RASTERS AND UNITS %%
                    if genPlots
                        % fplot figure for each unit
                        figure(un)

                        % plot any details.
                        subplot(3,3,1:2)
                        hold on
                        text(0,1.3,patientID,'fontweight','bold')
                        text(0,1.2,sprintf('microelectrode in %s, nearest macro: %s',chanLabel,macroLabel))
                        text(0,1.1,sprintf('macro Matter: %s', electrodeMatter))
                        %                         % logical index information
                        %                         if logical_IN(unitCount)
                        %                             text(0,1,sprintf('(classification: interneuron)'))
                        %                         elseif ~logical_IN(unitCount)
                        %                             text(0,1,sprintf('(classification: principal cell)'))
                        %                         end
                        % KMeans 2 index
                        %                         if KMeanTwo(unitCount) == 1
                        %                             text(0,0.9,sprintf('KMeansTwo Cluster: RED'), 'Color', 'r')
                        %                         elseif KMeanTwo(unitCount) == 2
                        %                             text(0,0.9,sprintf('KMeansTwo Cluster: BLUE'), 'Color', 'b')
                        %                         end
                        %                         % KMeans 3 index
                        %                         if KMeanThree(unitCount) == 1
                        %                             text(0,0.8,sprintf('KMeansThree Cluster: RED'), 'Color', 'r')
                        %                         elseif KMeanThree(unitCount) == 2
                        %                             text(0,0.8,sprintf('KMeansThree Cluster: BLUE'), 'Color', 'b')
                        %                         elseif KMeanThree(unitCount) == 3
                        %                             text(0,0.8,sprintf('KMeansThree Cluster: GREEN'), 'Color', 'g')
                        %                         end
                        %                         % isolated Cluster
                        %                         if isolatedCluster(unitCount) == 1
                        %                             text(0,0.7,sprintf('Isolated Cluster: TRUE'), 'Color', 'g')
                        %                         elseif isolatedCluster(unitCount) == 0
                        %                             text(0,0.7,sprintf('Isolated Cluster: FALSE'), 'Color', 'r')

                        %                         end

                        % Comparisons Stats
                        if pAllz > 0.05
                            text(0,0.6,sprintf('(All Pre/Post, p= ns)'))
                        elseif pAllz < 0.05
                            text(0,0.6,sprintf('(All Pre/Post, p= %s)',pAllz))
                        end
                        if pNearz > 0.05
                            text(0,0.5,sprintf('(Near Pre/Post, p= ns)'))
                        elseif pNearz < 0.05
                            text(0,0.5,sprintf('(Near Pre/Post, p= %s)',pNearz))
                        end
                        if pFarz > 0.05
                            text(0,0.4,sprintf('(Far Pre/Post, p= ns)'))
                        elseif pFarz < 0.05
                            text(0,0.4,sprintf('(Far Pre/Post, p= %s)',pFarz))
                        end
                        if pInstantz > 0.05
                            text(0,0.3,sprintf('(stimFR Pre/Post, p= ns)'))
                        elseif pInstantz < 0.05
                            text(0,0.3,sprintf('(stimFR Pre/Post, p= %s)',pInstantz))
                        end
                        hold off
                        axis off

                        % plot spike waveform
                        subplot(3,3,3)
                        hold on
                        plot(NEV.Data.Spikes.Waveform(:,ChanUnitTimestamp(:,1)==inclChans(ch) & ChanUnitTimestamp(:,2)==un)./4,'linewidth',0.1,'color',[0.651 0.651 0.651]) %lighter gray
                        plot(mean(NEV.Data.Spikes.Waveform(:,ChanUnitTimestamp(:,1)==inclChans(ch) & ChanUnitTimestamp(:,2)==un),2)./4,'linewidth',1,'color',[0.0196 0.2902 0.1647]) % dark green
                        hold off
                        axis tight square
                        xlabel('samples')
                        ylabel('waveform amplitude (\muV)')

                        %% plotting firing for stim near the microelectrodes
                        tmp = cellfun(@transpose, squeeze(struct2cell(spikes.channel(ch).unit(un).trial)), 'UniformOutput', false);
                        if nansum(nearStims)~=0 %changed from zero to isnan
                            %sum(nearStims)>5 || length(tmp(nearStims))>5 || ~isempty(wRpreNear)
                            % raster for nearby electrodes
                            subplot(3,3,4)
                            plotSpikeRaster(tmp(nearStims),'PlotType','vertline','rasterWindowOffset',-pre);
                            line([0 0],[0 1000],'color','k','linewidth',0.5,'linestyle','--')
                            clear tmp
                            axis xy square
                            ylabel('stimulations')
                            title(['stim near ' chanLabel])

                            % psth for nearby electrodes
                            subplot(3,3,7)
                            hold on
                            line([0 0],[min(RtrialsNear+EtrialsNear) max(RtrialsNear+EtrialsNear)],'color','k','linewidth',0.5,'linestyle','--')
                            patch([tsecN fliplr(tsecN)],[RtrialsNear+EtrialsNear fliplr(RtrialsNear-EtrialsNear)],[0.5 0.5 0.5],'facealpha',0.5,'edgecolor','none')
                            plot(tsecN,RtrialsNear,'color', [0.6350 0.0780 0.1840]) % dark red
                            hold off
                            xlim([-1.5 2.5])
                            axis square
                            xlabel('time (s)')
                            ylabel('firing rate (Hz)')
                        else
                            subplot(3,3,4)
                            title(['no stim near ' chanLabel])
                            axis off
                        end

                        % raster for distant electrodes
                        subplot(3,3,6)
                        %  MarkerFormat.MarkerSize = 1;
                        tmp = cellfun(@transpose, squeeze(struct2cell(spikes.channel(ch).unit(un).trial)), 'UniformOutput', false);
                        if sum(nearStims) == 0
                            plotSpikeRaster(tmp(~SOZstims),'PlotType','vertline','rasterWindowOffset',-pre);
                        elseif sum(nearStims) > 0
                            plotSpikeRaster(tmp(~nearStims & ~SOZstims),'PlotType','vertline','rasterWindowOffset',-pre);
                        end
                        line([0 0],[0 1000],'color','k','linewidth',0.5,'linestyle','--')
                        clear tmp
                        xlim([-1.5 2.5])
                        axis xy square
                        ylabel('stimulations')
                        title('stim elsewhere')

                        % psth for distant electrodes
                        subplot(3,3,9)
                        hold on
                        line([0 0],[min(RtrialsFar+EtrialsFar) max(RtrialsFar+EtrialsFar)],'color','k','linewidth',0.5,'linestyle','--') % doesnt work ..
                        patch([tsecF fliplr(tsecF)],[RtrialsFar+EtrialsFar fliplr(RtrialsFar-EtrialsFar)],[0.5 0.5 0.5],'facealpha',0.5,'edgecolor','none')
                        plot(tsecF,RtrialsFar,'color', [0.6350 0.0780 0.1840]) % dark red
                        hold off
                        xlim([-1.5 2.5])
                        axis square
                        xlabel('time (s)')
                        ylabel('firing rate (Hz)')

                        % save figures
                        %                halfMaximze(un,'left')
                        saveas(un,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\FiringRates\',sprintf('CCEPfiring_%s_ch%d_un%d.pdf',patientID,ch,un)))
                        close(un)
                    end
                end
            end
        end
        clear R spT nearStims SOZstims sumidxAll idxAll RtrialsNearz RtrialsFarz RtrialsAllz k D_euclidean trodeLocsXYZ nearMacros_EucD SEStr EucStr RstimAll tmpT tmpFR logical_AUC 
    end
end

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SUMMARY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%

% overall session and unit counts
fprintf('\n\n%d sessions/subjects\n',length(nevList));
fprintf('\n%d total units\n', length([SPESstruct(:)]));

keyboard

% Saving SPESstruct to main path.
save('\\155.100.91.44\D\Data\Rhiannon\CCEPS\SPES_UNITS\SPESstruct.mat','SPESstruct');

%keyboard % uncomment for new pt additions

%% ~~~~~~~~~~~~~~ ACROSS SUBJECTS SUMMARY OF UNIT CHARACTERISTICS ~~~~~~~~~~~~~~~~~ %%

SigAll = [SPESstruct(:).hAllz] == 1; % with zscore ttests.
sum(SigAll)

Sig_SPESstruct = [SPESstruct(SigAll)]; % subsetting SPES to only have Sig Units.
sig_logical_IN = [Sig_SPESstruct(:).logical_IN];

% Gray vs White Matter
electrodeMatterCell={SPESstruct.electrodeMatter}; % turn back to cell
electrodeMatterCell'; % transpose
grayMatterUnits = find(strcmp(electrodeMatterCell,'Gray')==1); % unit positions (75)
whiteMatterUnits = find(strcmp(electrodeMatterCell,'White')==1); % unit positions (110)
% # of significant gray and white matter units
SigGray = [SPESstruct(grayMatterUnits).hAll] == 1;
sum(SigGray)
SigWhite = [SPESstruct(whiteMatterUnits).hAll] == 1;
sum(SigWhite)
% props of gray and white contacts
propGray = (sum(SigGray)/length(grayMatterUnits))*100;
propWhite = (sum(SigWhite)/length(whiteMatterUnits))*100;

% subsetting SPESstruct to have only Gray or White Matter
Gray_SPESstruct = [SPESstruct(grayMatterUnits)];
White_SPESstruct = [SPESstruct(whiteMatterUnits)];

% logical Indexes for matterStructs
SigAll_Gray = [Gray_SPESstruct(:).hAll] == 1;
SigAll_White = [White_SPESstruct(:).hAll] == 1;

% number of near/far electrodes in white/gray matter
farGrayCount = sum([Gray_SPESstruct(:).hFar]== 0) + sum([Gray_SPESstruct(:).hFar]== 1);
farWhiteCount = sum([White_SPESstruct(:).hFar]== 0) + sum([White_SPESstruct(:).hFar]== 1);
nearGrayCount = sum([Gray_SPESstruct(:).hNear]== 0) + sum([Gray_SPESstruct(:).hNear]== 1);
nearWhiteCount  = sum([White_SPESstruct(:).hNear]== 0) + sum([White_SPESstruct(:).hNear]== 1);


%% VENN DIAGRAMS FOR SIGNIFICANT RESPONSES:
% ~~~~~~~~~~~~~~~~~ (1) All Neurons ~~~~~~~~~~~~~~~~~ %

Near = sum(SigNear)
Far = sum(SigFar)
Instant = sum(SigInstant)

% (1) Finding the intersections of Near Far arrays
SigStim_NearFar = sum(SigStimsFN) % summing across columns to see which contacts response for both N and F (=2).
Intersection_NearFar = sum(SigStim_NearFar(:) == 2) % total # of contacts that are both N and F.
fprintf("Intersection NearFar: ")
disp(Intersection_NearFar);
fprintf("\n");

% (2) Finding the intersections of Instantaneous and Far arrays
SigStim_FarInstant = sum(SigStimsFI) % summing across columns to see which contacts response for both I and F (=2).
Intersection_FarInstant = sum(SigStim_NearFar(:) == 2) % total # of contacts that are both I and F.
fprintf("Intersection FarInstant: ")
disp(Intersection_FarInstant);
fprintf("\n");

% (3) Finding the intersections of Instantaneous and Near arrays
SigStim_NearInstant = sum(SigStimsNI) % summing across columns to see which contacts response for both N and I (=2).
Intersection_NearInstant = sum(SigStim_NearInstant(:) == 2) % total # of contacts that are both N and I.
fprintf("Intersection NearInstant: ")
disp(Intersection_NearInstant);
fprintf("\n");

% (4) Finding the intersections of Instantaneous, Far, and Near arrays
SigStim_FarNearInstant = sum(SigStimsFNI) % summing across columns to see which contacts response for both F and N and I (= 3). doenst work if zero
Intersection_FarNearInstant = sum(SigStim_FarNearInstant(:) == 3) % total # of contacts that are both N and I.
fprintf("Intersection FarNearInstant: ")
disp(Intersection_FarNearInstant);
fprintf("\n");
clc;

clf;
figure(1238)
Plotting_Interval = 0.01;
Angles_In_Radians = (0: Plotting_Interval: 2*pi);
Circle_Plot = @(X_Offset,Y_Offset,Radius) plot(X_Offset + Radius*cos(Angles_In_Radians),Y_Offset + Radius*sin(Angles_In_Radians));

hold on
%Plotting the 3 circles%
X_Offset_Near = 0; Y_Offset_Near = 2; Radius_Near = 3;
Circle_Near = Circle_Plot(X_Offset_Near,Y_Offset_Near,Radius_Near);
fill(Circle_Near.XData, Circle_Near.YData,'r','FaceAlpha',0.2,'LineWidth',1);

X_Offset_Far = -2; Y_Offset_Far = -2; Radius_Far = 3;
Circle_Far = Circle_Plot(X_Offset_Far,Y_Offset_Far,Radius_Far);
fill(Circle_Far.XData, Circle_Far.YData,'g','FaceAlpha',0.2,'LineWidth',1);

X_Offset_Instant = 2; Y_Offset_Instant = -2; Radius_Instant = 3;
Circle_Plot(X_Offset_Instant,Y_Offset_Instant,Radius_Instant);
Circle_Instant = Circle_Plot(X_Offset_Instant,Y_Offset_Instant,Radius_Instant);
fill(Circle_Instant.XData, Circle_Instant.YData,'b','FaceAlpha',0.2,'LineWidth',1);
title("Venn Diagram");

%Writing all the labels%
Near_Label = strjoin(string(Near));
text(X_Offset_Near,Y_Offset_Near,Near_Label,'color','r');

Far_Label = strjoin(string(Far));
text(X_Offset_Far,Y_Offset_Far,Far_Label,'color','g');

Instant_Label = strjoin(string(Instant));
text(X_Offset_Instant,Y_Offset_Instant,Instant_Label,'color','b');

NearFar_Label = strjoin(string(Intersection_NearFar));
text(-1.2,0,NearFar_Label);

FarInstant_Label = strjoin(string(Intersection_FarInstant));
text(0,-2,FarInstant_Label);

NearInstant_Label = strjoin(string(Intersection_NearInstant));
text(1.2,0,NearInstant_Label);

FarNearInstant_Label = strjoin(string(Intersection_FarNearInstant));
text(0,0,FarNearInstant_Label);

%Setting the labels to be relative to the centres%
set(findall(gcf,'type','text'),'HorizontalAlignment','center');
axis equal
axis off

% saving venn diagram figure
saveas(1238,fullfile('D:\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\',sprintf('Unit_Venn_Diagram_AcrossSubs.pdf')))
close(1238)
