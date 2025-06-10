
%% ~~~~~~~~~~~~~~~~~~~~~~~~   READ ME:  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%

% This code uses a waveform classification (dependent of function
% waveformMetrics) to classify waveforms as interneurons or principal
% cells. We used three cluster methods: (1) GMM, (2) K-Means, (3) metric
% thresholds to identify the neuron type. 

% This code also produces figures for the manuscript and supplementary as
% well as additional figures that did not go into the script. Figures that
% are included in the paper are indicated. 

%% load SPESstruct data
load('\\155.100.91.44\D\Data\Rhiannon\CCEPS\SPES_UNITS\SPESstruct.mat')

% the number of units should be the length of SPESstruct, unless we have
% empty fields.

%% loop over units and run Ed's analyses. 
for un = 1:length(SPESstruct)
    meanWaveforms(un,:) = mean(SPESstruct(un).spikeWaveforms,2)';
end

% zscore waveforms to try
zWFs = zscore(meanWaveforms')';
%plot(zWFs')

%% now running ED's Code:

% meanWaveforms should be an [n x m] matrix of n total mean waveforms
% across m data samples. I feed in 401 data points for each usually, with
% the trough at the 200th data point, because the function up-samples and
% re-finds the trough, so I give it lots of leeway. Would be fine to give
% much less, but will need to update the keyPoint input if the trough isn't
% at 200.
[metrics,outputWaves] = waveformMetrics(meanWaveforms);

%% First option: fit a GMM to the three main waveform metrics:
% (depending on the data, play around with which inputs)
inputData = [metrics.troughToPeak; metrics.FWHM; metrics.asymmetry]';
gm = fitgmdist(inputData,2,'Replicates',1,'SharedCovariance',true);
probs = gm.posterior(inputData);
% work out which component had the shorter duration waveforms:
muFWHM = [nanmean(metrics.FWHM(probs(:,1) > 0.5)) nanmean(metrics.FWHM(probs(:,2) > 0.5))];
[~,isIN] = min(muFWHM);
FSprobs = probs(:,isIN);

figure(1)
scatter(metrics.troughToPeak, metrics.FWHM, 60, FSprobs,'filled')
xlabel('Trough to peak duration (ms)')
ylabel('FWHM (ms)')
cb = colorbar;
ylabel(cb,'FS interneuron probability')
colormap cool
grid on

% saving figure
saveas(1,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\WaveFormClassification\',sprintf('GMM_WF_Metrics_AcrossSubs.pdf')))
close(1)

%% Alternatively, maybe just fitting a GMM to the first N PCA scores on the
% z-scored waveforms does better, depending on how clean the data are:
% (This approach works better than the metrics on the data I've been using)
[~,pc,~,~,expl] = pca(zscore(outputWaves')');
% percentage variance explained required (if not working, try lowering
% this, since a high amount of variance in the waveforms is likely to be
% within cell-types rather than across groups, meaning including them for
% this is just incorporating inherent waveform noise rather than
% differences between cell types)
explReq = 85; % threshold for explained variance how mayn pcs need to  explain 85% and basing GMM on that.
n = find(cumsum(expl) >= explReq,1);
inputData = [metrics.troughToPeak; metrics.FWHM; metrics.asymmetry; pc(:,1:n)']'; % first three PCs. 
gm = fitgmdist(inputData,2,'Replicates',100,'SharedCovariance',true); % changed from 100
%gm = fitgmdist(pc(:,1:n),2,'Replicates',100,'SharedCovariance',true); % if we want to run the GM on the PC components only.
probs = gm.posterior(inputData);
% work out which component had the shorter duration waveforms:
muFWHM = [nanmean(metrics.FWHM(probs(:,1) > 0.5)) nanmean(metrics.FWHM(probs(:,2) > 0.5))];
[~,isIN] = min(muFWHM);
FSprobs = probs(:,isIN);

figure(2)
scatter(metrics.troughToPeak, metrics.FWHM, 60, FSprobs,'filled')
xlabel('Trough to peak duration (ms)')
ylabel('FWHM (ms)')
cb = colorbar;
ylabel(cb,'FS interneuron probability')
colormap cool
grid on

% saving figure
saveas(2,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\WaveFormClassification\',sprintf('GMM_WFPCA_Zscored_WFs_AcrossSubs.pdf')))
close(2)

% additional figures for visualizing WF metrics.
% Trough To Peak vs. Asymmetry.
figure(10)
scatter(metrics.troughToPeak, metrics.asymmetry, 60, FSprobs,'filled')
xlabel('Trough to peak duration (ms)')
ylabel('Asymmetry')
cb = colorbar;
ylabel(cb,'FS interneuron probability')
colormap cool
grid on
% saving figure
saveas(10,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\WaveFormClassification\',sprintf('GMM_WFPCA_Zscored_TTP_ASY_AcrossSubs.pdf')))
close(10)

% FWHM vs. Asymmetry
figure(11)
scatter(metrics.FWHM, metrics.asymmetry, 60, FSprobs,'filled')
xlabel('Asymmetry')
ylabel('FWHM (ms)')
cb = colorbar;
ylabel(cb,'FS interneuron probability')
colormap cool
grid on

% saving figure
saveas(11,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\WaveFormClassification\',sprintf('GMM_WFPCA_Zscored_ASY_FWHM_AcrossSubs.pdf')))
close(11)

% TroughToPeak vs. FWHM vs. Asymmetry
figure(12)
scatter3(metrics.troughToPeak, metrics.FWHM, metrics.asymmetry, 60, FSprobs,'filled')
xlabel('Trough To Peak duration (ms)')
ylabel('FWHM (ms)')
zlabel('Asymmetry')
cb = colorbar;
ylabel(cb,'FS interneuron probability')
colormap cool
grid on

% saving figure
saveas(12,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\WaveFormClassification\',sprintf('GMM_WFPCA_Zscored_3D_AcrossSubs.pdf')))
close(12)

%% 3rd option: worth trying a GMM on just the trough to peak and the log of the FWHM values in the meantime? 
%logFWHM = log(metrics.FWHM)
%gm = fitgmdist(pc(:,1:n),2,'Replicates',100,'SharedCovariance',true);
%inputData = [metrics.troughToPeak; metrics.FWHM; metrics.asymmetry]';
%gm = fitgmdist(inputData,2,'Replicates',100,'SharedCovariance',true);

% That third cell-intrinsic firing metric also (sometimes) picks up non-FS INs 
% so it could be messing up the separation on the others, and FWHM distributions generally tend towards log-normal, 
% so probably worth making that transformation prior to fitting the GMM and seeing what happens?

% might need a few more units for this. 
% intertrial coherence / shifts from pre to post?

%% ALternative Option: K-MEANS

 % k-means calculation
[idxTwo,CTwo] = kmeans(inputData,2); %returns distances from each point to every centroid in the n-by-k matrix D.
[idxThree,CThree] = kmeans(inputData,3); %returns distances from each point to every centroid in the n-by-k matrix D.

% mean metrics for k2 clusters
meanTTP_k2clusterRed = mean(inputData(idxTwo==1,1));
meanFWHM_k2clusterRed = mean(inputData(idxTwo==1,2));
meanTTP_k2clusterBlue = mean(inputData(idxTwo==2,1));
meanFWHM_k2clusterBlue = mean(inputData(idxTwo==2,2));

% mean metrics for k3 clusters
meanTTP_k3clusterRed = mean(inputData(idxThree==1,1));
meanFWHM_k3clusterRed = mean(inputData(idxThree==1,2));
meanTTP_k3clusterBlue = mean(inputData(idxThree==2,1));
meanFWHM_k3clusterBlue = mean(inputData(idxThree==2,2));
meanTTP_k3clusterGreen = mean(inputData(idxThree==3,1));
meanFWHM_k3clusterGreen = mean(inputData(idxThree==3,2));

% plots for k-means
figure(3);
plot(inputData(idxTwo==1,1),inputData(idxTwo==1,2),'r.','MarkerSize',20)
hold on
plot(inputData(idxTwo==2,1),inputData(idxTwo==2,2),'b.','MarkerSize',20)
plot(CTwo(:,1),CTwo(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
title('K2: Cluster Assignments and Centroids','FontSize',14)
hold off

KMeanTwo = idxTwo; % 1 = red, 2 = blue

% saving figure
saveas(3,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\WaveFormClassification\',sprintf('KMeansTwo_WFClassification_AcrossSubs.pdf')))
close(3)

noWfs = false;
if noWfs
    figure(4);
    plot(inputData(idxThree==1,1),inputData(idxThree==1,2),'r.','MarkerSize',20)
    hold on
    plot(inputData(idxThree==2,1),inputData(idxThree==2,2),'b.','MarkerSize',20)
    plot(inputData(idxThree==3,1),inputData(idxThree==3,2),'g.','MarkerSize',20)
    plot(CThree(:,1),CThree(:,2),'kx',...
        'MarkerSize',15,'LineWidth',3)
    legend('Cluster 1','Cluster 2', 'Cluster 3','Centroids',...
        'Location','NW')
    title('K3: Cluster Assignments and Centroids','FontSize',14)
    hold off

    KMeanThree = idxThree; % 1 = red, 2 = blue, 3 = green

    % saving figure
    saveas(4,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\WaveFormClassification\',sprintf('KMeansThree_WFClassification_AcrossSubs.pdf')))
    close(4)
else
    % to plot the mean waveforms across clusters. 
    figure(4);
    subplot(2,3,2)
    plot(inputData(idxThree==1,1),inputData(idxThree==1,2),'r.','MarkerSize',20)
    hold on
    plot(inputData(idxThree==2,1),inputData(idxThree==2,2),'b.','MarkerSize',20)
    plot(inputData(idxThree==3,1),inputData(idxThree==3,2),'g.','MarkerSize',20)
    plot(CThree(:,1),CThree(:,2),'kx',...
        'MarkerSize',15,'LineWidth',3)
    legend('Cluster 1','Cluster 2', 'Cluster 3','Centroids',...
        'Location','NW')
    title('K3: Cluster Assignments and Centroids','FontSize',14)
    xlabel('Trough-To-Peak (ms)','FontSize',12)
    ylabel('FWHM (ms)','FontSize',12)
    hold off
    axis tight square

    KMeanThree = idxThree; % 1 = red, 2 = blue, 3 = green
    % red cluster
    subplot(2,3,4)
    plot(meanWaveforms(idxThree==1,:)','r-')
    axis tight square
    xlabel('samples')
    numObservationsRed = sum(idxThree == 1);
    title(sprintf('# of Waveforms: %d\nMean TroughToPeak = %.2f ms\nMean FWHM = %.2f ms',numObservationsRed, meanTTP_k3clusterRed, meanFWHM_k3clusterRed));
    % blue cluster
    subplot(2,3,5)
    plot(meanWaveforms(idxThree==2,:)','b-')
    axis tight square
    xlabel('samples')
    numObservationsBlue = sum(idxThree == 2);
    title(sprintf('# of Waveforms: %d\nMean TroughToPeak = %.2f ms\nMean FWHM = %.2f ms', numObservationsBlue,meanTTP_k3clusterBlue, meanFWHM_k3clusterBlue));
    % green cluster
    subplot(2,3,6)
    plot(meanWaveforms(idxThree==3,:)','g-')
    axis tight square
    xlabel('samples')
    numObservationsGreen = sum(idxThree == 3);
    title(sprintf('# of Waveforms: %d\nMean TroughToPeak = %.2f ms\nMean FWHM = %.2f ms', numObservationsGreen,meanTTP_k3clusterGreen, meanFWHM_k3clusterGreen));

    % saving figure
    maximize(4)
    saveas(4,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\WaveFormClassification\',sprintf('KMeansThree_WFs_AcrossSubs.pdf')))
    close(4)

end

%% Supplementary Plot: Clusters for different data segregations %%

figure(5)
% plotting wfs by patient
patientID = {SPESstruct.patientID}; % extract ptID from structure
patientID'; % transpose
[patientNumber, ~, icPatientNumber] = unique(patientID);

subplot(2,2,1)
scatter(metrics.troughToPeak, metrics.FWHM, 80, icPatientNumber, "filled", 'MarkerEdgeColor', 'k')
colormap colorcube
xlabel('Trough to peak duration (ms)','FontSize',12)
ylabel('FWHM (ms)','FontSize',12)
title('Indiviudal Patient Responses','FontSize',14)
axis square
grid on

% plot by region
unitLocs = {SPESstruct(:).chanLabel}; % get unit locations
unitLocs = string(unitLocs); % changing to string
AMY = "Amygdala"; % finding locations with Amygdala
AMYLocs = contains(unitLocs(1,:), AMY);
OFC = ("OFC"| "Orbital" | "vmPFC"); % finding locations with OFC
OFCLocs = contains(unitLocs(1,:), OFC);
HIPP = "Hippocampus"; % finding locations with HIPP
HIPPLocs = contains(unitLocs(1,:), HIPP);
CING = "Cingulate"; % finding locations with Cingulate
CINGLocs = contains(unitLocs(1,:), CING);
pCING = "Posterior"; % finding locations with Cingulate
pCINGLocs = contains(unitLocs(1,:), pCING);
Regions = [AMYLocs; OFCLocs; HIPPLocs; CINGLocs; pCINGLocs];
RegionLabels = {'Amygdala', 'OFC', 'Hippocampus', 'Cingulate', 'Posterior Cingulate'};
RegionColors = lines(length(RegionLabels)); % Generate different colors for each region

subplot(2,2,2)
hold on
for i = 1:size(Regions, 1)
    scatter(metrics.troughToPeak(Regions(i,:)), metrics.FWHM(Regions(i,:)), 80, RegionColors(i,:), "filled",'MarkerEdgeColor', 'k','DisplayName', RegionLabels{i})
end
hold off
% labels
xlabel('Trough to peak duration (ms)', 'FontSize', 12)
ylabel('FWHM (ms)','FontSize',12)
grid on
legend('show')
title('Regional responses','FontSize',14)
axis square

% Plotting significant vs nonsignificant units
sigAll = [SPESstruct.hAll]; % extract significant unit from structure
sigAll = logical(sigAll)'; % transpose
[uniquesigAll, ~, icSigAll] = unique(sigAll);
Significant = logical(icSigAll == 2);
NonSignificant = logical(icSigAll == 1);
Units = [Significant; NonSignificant];
UnitLabels = {'Significant', 'Non-Significant'};

subplot(2,2,3)
hold on
scatter(metrics.troughToPeak(icSigAll==2), metrics.FWHM(icSigAll==2), 80, 'b', "filled", 'MarkerEdgeColor', 'k')
scatter(metrics.troughToPeak(icSigAll==1), metrics.FWHM(icSigAll==1), 80, 'k', "filled", 'MarkerEdgeColor', 'k')
hold off
xlabel('Trough to peak duration (ms)','FontSize',12)
ylabel('FWHM (ms)','FontSize',12)
title('Sig or NS Unit Responses','FontSize',14)
legend('Significant Response', 'No Response'); 
grid on
axis square

% plot by suppressed vs excited FRs
dataAll = [SPESstruct.RdiffAll];
suppressionAll = dataAll < 0; % 1 = suppression, 0 = excitatory
[suppressionAll, ~, icsuppressionAll] = unique(suppressionAll); % 2 = suppression, 1 = excitatory

subplot(2,2,4) % 
hold on;
scatter(metrics.troughToPeak(icsuppressionAll == 2), metrics.FWHM(icsuppressionAll == 2), 80, 'o', 'filled', 'MarkerEdgeColor', 'k'); % Blue for suppressed units
scatter(metrics.troughToPeak(icsuppressionAll == 1), metrics.FWHM(icsuppressionAll == 1), 80, 'red', 'filled', 'MarkerEdgeColor', 'k'); % Orange for excitory
scatter(metrics.troughToPeak(NonSignificant), metrics.FWHM(NonSignificant), 80, 'k', 'filled', 'MarkerEdgeColor', 'k'); % Black for Nonsignifcant responses
hold off;
xlabel('Trough to peak duration (ms)','FontSize',12);
ylabel('FWHM (ms)','FontSize',12);
title('Suppressed or Enhanced Unit Clusters','FontSize',14);
grid on; 
legend('Suppressed Response', 'Enhanced Response', 'No Response'); 
axis square
maximize(figure(5))
% saving figure
saveas(5,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\WaveFormClassification\',sprintf('byPtRegSigMod_WFs_AcrossSubs.pdf')))
close(5)

%% Figure for SOZ cluster %%

figure(51)
% plotting wfs by patient
isSOZ = [SPESstruct.isSOZ]; % extract ptID from structure
isSOZ'; % transpose

hold on
scatter(metrics.troughToPeak(isSOZ==1), metrics.FWHM(isSOZ==1), 60, 'yellow', "filled", 'MarkerEdgeColor', 'k')
scatter(metrics.troughToPeak(isSOZ==0), metrics.FWHM(isSOZ==0), 60, 'k', "filled", 'MarkerEdgeColor', 'k')
hold off
xlabel('Trough to peak duration (ms)','FontSize',12)
ylabel('FWHM (ms)','FontSize',12)
title('SOZ unit microrecording','FontSize',14)
legend('SOZ unit', 'other unit'); 
grid on
axis square
saveas(51,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\WaveFormClassification\',sprintf('bySOZ_WFs_AcrossSubs.pdf')))
close(51)

%% For "wide" interneurons, I remove the units that have been deemed FS by
% the above approaches, then use fit_ACG from the CellExplorer toolbox to
% calculate their "tau rise" component in the exponentials, and fit a
% 2-component GMM to that.
% For example, I usually have my single units stored in my "NeuroClass"
% objects meaning I'd do this to calculate the autocorrelations, though any
% approach for that would work:
lags = -50:0.5:50;
ac = NaN(length(lags),length(SPESstruct));
for u = 1:length(SPESstruct)
    ac(:,u) = autocorr(SPESstruct(u).unitTimes,lags);
    ac(:,u) = ac(:,u)./length(SPESstruct(u).unitTimes)./0.0005;
    % Note that the original CCG approach used 'norm','rate' so need to
    % convert from 'count' to 'rate'. This does that correctly, matching
    % their calculated rates, but it makes me suspicious that their rate is
    % not "spikes per second"
end
% Now fit the exponentials:
params = fit_ACG(ac,false);

% And fit the GMM. I found acg_tau_rise was the most useful, because it's
% the one that captures whether the unit fires with a steady increase in
% probability over the first 10s of milliseconds, which is stereotypical of
% interneurons. The burst one should do a good job of finding the subset of
% excitatory pyramidal cells that burst-fire.
gmAC = fitgmdist(params.acg_tau_rise',2,'Replicates',100);
wideProbs = gmAC.posterior(params.acg_tau_rise');
muACG = [nanmean(params.acg_tau_rise(wideProbs(:,1) > 0.5)) nanmean(params.acg_tau_rise(wideProbs(:,2) > 0.5))];
[~,isW] = max(wideProbs); % find which group had the higher tau_acg, and is therefore more likely to be the wide interneurons
wideProbs = wideProbs(isW,:); % [20230814EHS] switched the indexing dimension to fix an error, but not sure if the answer is correct. 

%% Isolating small fast response cluster %%

isolatedCluster = metrics.troughToPeak <= 0.35 & metrics.FWHM <= 0.35; % isolating just the fast FWHM and fast TTP units

% Plotting isolated cluster
% Define custom colors
sage_green = [119, 171, 86] / 255; % RGB values for sage green
muted_purple = [79, 49, 170] / 255; % RGB values for muted purple

figure(9)
hold on; 
scatter(metrics.troughToPeak(isolatedCluster == 1), metrics.FWHM(isolatedCluster == 1), 70, sage_green, 'filled', 'MarkerEdgeColor', 'k'); % Red for FS interneuron
scatter(metrics.troughToPeak(isolatedCluster == 0), metrics.FWHM(isolatedCluster == 0), 70, muted_purple, 'filled', 'MarkerEdgeColor', 'k'); % Black for Principal Cell
xlabel('Trough to peak duration (ms)','FontSize',12);
ylabel('FWHM (ms)','FontSize',12);
title('Cell Type Clusters','FontSize',14);
grid on; 
legend('FS interneuron', 'Principal Cell'); 
hold off; 

% saving figure
saveas(9,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\WaveFormClassification\',sprintf('byIsolatedCLuster_Zscored_WFs_AcrossSubs.pdf')))
close(9) 

 %% ~~~~~~~~~~~~~~~~~ STATISTICS FOR IN, PCS and SOZs ~~~~~~~~~~~~~~~~~~~~~~~ %%
 
 % ~~~~~~~~ Quanitifying Responses ~~~~~~~~~~~ %

 % count and prop of isolated FS Interneuron cluster.
 INcount = sum(isolatedCluster); % 37
 PCcount = sum(~isolatedCluster); % 191
 SOZcount = sum(isSOZ); %27
 inSOZcount = sum(isolatedCluster & isSOZ) % 5
 pcSOZcount = sum(~isolatedCluster & isSOZ) % 22
 INnonSOZcount = sum(isolatedCluster & ~isSOZ) % 32
 PCnonSOZcount = sum(~isolatedCluster & ~isSOZ) % 169

 % logical indexes
 inSOZlogical = logical(isolatedCluster & isSOZ) % 5
 pcSOZlogical = logical(~isolatedCluster & isSOZ) % 22
 INnonSOZlogical = logical(isolatedCluster & ~isSOZ) % 32
 PCnonSOZlogical = logical(~isolatedCluster & ~isSOZ) % 169

 INperc = INcount/length(isolatedCluster)*100; % 16.23%
 PCperc = PCcount/length(~isolatedCluster)*100; % 83.77%
 SOZprop = SOZcount/length(isSOZ)*100; % 11.84% 

 % proportion test for SOZ in IN and PCs
 % statistics.
 INprop = sum(isolatedCluster & isSOZ); % 13.51% (5/37)
 PCprop = sum(~isolatedCluster & isSOZ); % 11.51% (22/191)  
 totals = [sum(isolatedCluster) sum(~isolatedCluster)];
 [hsig,psig,chisig,dgsig] = prop_test([INprop PCprop], totals,false); % not significant differences in proportions on SOZ INs and PCs.

% Mean and SD for each Cluster (a couple of NaNs because of inverted waveforms)
meanTTP_clusterPurple = mean(inputData(isolatedCluster==0,1), 'omitnan'); % mean TTP PC
stdTTP_clusterPurple = std(inputData(isolatedCluster==0,1), 'omitnan'); % std TTP PC
meanFWHM_clusterPurple = mean(inputData(isolatedCluster==0,2), 'omitnan'); % mean FWHM PC
stdFWHM_clusterPurple = std(inputData(isolatedCluster==0,2), 'omitnan'); % std FWHM PC
meanASYM_clusterPurple = mean(inputData(isolatedCluster==0,3), 'omitnan'); % mean ASYM PC
stdASYM_clusterPurple = std(inputData(isolatedCluster==0,3), 'omitnan'); % std ASYM PC
meanTTP_clusterGreen = mean(inputData(isolatedCluster==1,1), 'omitnan'); % mean TTP IN
stdTTP_clusterGreen = std(inputData(isolatedCluster==1,1), 'omitnan'); % std TTP IN
meanFWHM_clusterGreen = mean(inputData(isolatedCluster==1,2), 'omitnan'); % mean FWHM IN
stdFWHM_clusterGreen = std(inputData(isolatedCluster==1,2), 'omitnan'); % std FWHM IN
meanASYM_clusterGreen = mean(inputData(isolatedCluster==1,3), 'omitnan'); % mean ASYM IN
stdASYM_clusterGreen = std(inputData(isolatedCluster==1,3), 'omitnan'); % std ASYM IN

% Mean and SD for SOZ (a couple of NaNs because of inverted waveforms)
meanTTP_SOZ = mean(inputData(isSOZ,1), 'omitnan'); % mean TTP SOZ
stdTTP_SOZ = std(inputData(isSOZ,1), 'omitnan'); % std TTP SOZ
meanFWHM_SOZ = mean(inputData(isSOZ,2), 'omitnan'); % mean FWHM SOZ
stdFWHM_SOZ = std(inputData(isSOZ,2), 'omitnan'); % std FWHM SOZ
meanASYM_SOZ = mean(inputData(isSOZ,3), 'omitnan'); % mean ASYM SOZ
stdASYM_SOZ = std(inputData(isSOZ,3), 'omitnan'); % std ASYM SOZ
% nonSOZ
meanTTP_nonSOZ = mean(inputData(~isSOZ,1), 'omitnan'); % mean TTP nonSOZ
stdTTP_nonSOZ = std(inputData(~isSOZ,1), 'omitnan'); % std TTP nonSOZ
meanFWHM_nonSOZ = mean(inputData(~isSOZ,2), 'omitnan'); % mean FWHM nonSOZ
stdFWHM_nonSOZ = std(inputData(~isSOZ,2), 'omitnan'); % std FWHM nonSOZ
meanASYM_nonSOZ = mean(inputData(~isSOZ,3), 'omitnan'); % mean ASYM nonSOZ
stdASYM_nonSOZ = std(inputData(~isSOZ,3), 'omitnan'); % std ASYM nonSOZ

% Cell Type Matrix
allTTP_clusterPurple = inputData(isolatedCluster==0,1); % TTP PC
allTTP_clusterPurple = allTTP_clusterPurple(~isnan(allTTP_clusterPurple));
allFWHM_clusterPurple = inputData(isolatedCluster==0,2); % FWHM PC
allFWHM_clusterPurple = allFWHM_clusterPurple(~isnan(allFWHM_clusterPurple));
allASYM_clusterPurple = inputData(isolatedCluster==0,3); % ASYM PC
allASYM_clusterPurple = allASYM_clusterPurple(~isnan(allASYM_clusterPurple));
allTTP_clusterGreen = inputData(isolatedCluster==1,1); % TTP IN
allTTP_clusterGreen = allTTP_clusterGreen(~isnan(allTTP_clusterGreen));
allFWHM_clusterGreen = inputData(isolatedCluster==1,2); % FWHM IN
allFWHM_clusterGreen = allFWHM_clusterGreen(~isnan(allFWHM_clusterGreen));
allASYM_clusterGreen = inputData(isolatedCluster==1,3); % ASYM IN
allASYM_clusterGreen = allASYM_clusterGreen(~isnan(allASYM_clusterGreen));

% SOZ Matrix (ALL)
allTTP_SOZ = inputData(isSOZ,1); % TTP SOZ
allTTP_SOZ = allTTP_SOZ(~isnan(allTTP_SOZ));
allFWHM_SOZ = inputData(isSOZ,2); % FWHM SOZ
allFWHM_SOZ = allFWHM_SOZ(~isnan(allFWHM_SOZ));
allASYM_SOZ = inputData(isSOZ,3); % ASYM SOZ
allASYM_SOZ = allASYM_SOZ(~isnan(allASYM_SOZ));
% non-SOZ (ALL)
allTTP_nonSOZ = inputData(~isSOZ,1); % TTP nonSOZ
allTTP_nonSOZ = allTTP_nonSOZ(~isnan(allTTP_nonSOZ));
allFWHM_nonSOZ = inputData(~isSOZ,2); % FWHM nonSOZ
allFWHM_nonSOZ = allFWHM_nonSOZ(~isnan(allFWHM_nonSOZ));
allASYM_nonSOZ = inputData(~isSOZ,3); % ASYM nonSOZ
allASYM_nonSOZ = allASYM_nonSOZ(~isnan(allASYM_nonSOZ));

% SOZ Matrix (just PC)
allTTP_SOZpc = inputData(isSOZ & ~isolatedCluster,1); % TTP SOZ
allTTP_SOZpc = allTTP_SOZpc(~isnan(allTTP_SOZpc));
allFWHM_SOZpc = inputData(isSOZ & ~isolatedCluster,2); % FWHM SOZ
allFWHM_SOZpc = allFWHM_SOZpc(~isnan(allFWHM_SOZpc));
allASYM_SOZpc = inputData(isSOZ & ~isolatedCluster,3); % ASYM SOZ
allASYM_SOZpc = allASYM_SOZpc(~isnan(allASYM_SOZpc));
% non-SOZ (just pc)
allTTP_nonSOZpc = inputData(~isSOZ & ~isolatedCluster,1); % TTP nonSOZ
allTTP_nonSOZpc = allTTP_nonSOZpc(~isnan(allTTP_nonSOZpc));
allFWHM_nonSOZpc = inputData(~isSOZ & ~isolatedCluster,2); % FWHM nonSOZ
allFWHM_nonSOZpc = allFWHM_nonSOZpc(~isnan(allFWHM_nonSOZpc));
allASYM_nonSOcZpc = inputData(~isSOZ & ~isolatedCluster,3); % ASYM nonSOZ
allASYM_nonSOZpc = allFWHM_nonSOZpc(~isnan(allFWHM_nonSOZpc));

% statistics between mean WF metrics
[TTPp,TTPh,TTPstats] = ranksum(inputData(isolatedCluster==0,1),inputData(~isolatedCluster==0,1)); % Trough-to-peak(SIG)
[FWHMp,FWHMh, FWHMstats] = ranksum(inputData(isolatedCluster==0,2),inputData(~isolatedCluster==0,2)); % Full-width(SIG)
[ASYMp,ASYMh,ASYMstats] = ranksum(inputData(isolatedCluster==0,3),inputData(~isolatedCluster==0,3)); % Asymmetry(SIG)

% statistics between SOZ or nonSOZ metrics
[TTP_SOZp,TTP_SOZh, TTP_SOZstats] = ranksum(inputData(isSOZ,1),inputData(~isSOZ,1)); % Trough-to-peak (SIG)
[FWHM_SOZp,FWHM_SOZh, FWHM_SOZstats] = ranksum(inputData(isSOZ,2),inputData(~isSOZ,2)); % Full-width (Non-SIG)
[ASYM_SOZp,ASYM_SOZh, ASYM_SOZstats] = ranksum(inputData(isSOZ,3),inputData(~isSOZ,3)); % Asymmetry (Non-SIG)

% statistics between SOZ or nonSOZ Principal Cell metrics
[TTP_pcSOZp,TTP_pcSOZh, TTP_pcSOZstats] = ranksum(inputData(isSOZ & ~isolatedCluster,1),inputData(~isSOZ & ~isolatedCluster,1)); % Trough-to-peak (SIG)
[FWHM_pcSOZp,FWHM_pcSOZh, FWHM_pcSOZstats] = ranksum(inputData(isSOZ & ~isolatedCluster,2),inputData(~isSOZ & ~isolatedCluster,2)); % Full-width (Non-SIG)
[ASYM_pcSOZp,ASYM_pcSOZh, ASYM_pcSOZstats] = ranksum(inputData(isSOZ & ~isolatedCluster,3),inputData(~isSOZ & ~isolatedCluster,3)); % Asymmetry (Non-SIG)

% TABLE TO SHOW ISOLATED CLUSTER METRICS
data = { meanTTP_clusterPurple, stdTTP_clusterPurple, meanFWHM_clusterPurple, stdFWHM_clusterPurple, meanASYM_clusterPurple, stdASYM_clusterPurple;
    meanTTP_clusterGreen, stdTTP_clusterGreen, meanFWHM_clusterGreen, stdFWHM_clusterGreen, meanASYM_clusterGreen, stdASYM_clusterGreen};
colNames = {'Trough To Peak M (ms)', 'Trough To Peak SD (ms)', 'FWHM M (ms)', 'FWHM Sd (ms)', 'Asymmetry M', 'Asymmetry SD'};
rowNames = {'Principal Cell', 'Interneuron'};
clusterTable = cell2table(data, 'VariableNames', colNames, 'RowNames', rowNames);
disp(clusterTable);

% SCATTER AND BARPLOT OF METRICS
meanMetrics = [meanTTP_clusterPurple, meanFWHM_clusterPurple, meanASYM_clusterPurple; meanTTP_clusterGreen, meanFWHM_clusterGreen, meanASYM_clusterGreen];
sdMetrics = [stdTTP_clusterPurple, stdFWHM_clusterPurple, stdASYM_clusterPurple; stdTTP_clusterGreen,stdFWHM_clusterGreen,stdASYM_clusterGreen]; 

% Define colors
sage_green = [119, 171, 86] / 255; % RGB values for sage green
muted_purple = [79, 49, 170] / 255; % RGB values for muted purple
boxColors = [muted_purple; sage_green];

isolatedCluster = isolatedCluster';

%% ~~~~~~~~~~~~~ MANUSCRIPT FIGURE 2: WAVEFORM CLASSIFICATION ~~~~~~~~~~~~
figure(10)
% Cluster Scatter Plot
subplot(2,4,1);
hold on;
scatter(metrics.troughToPeak(isolatedCluster == 1), metrics.FWHM(isolatedCluster == 1), 70, sage_green, 'filled', 'MarkerEdgeColor', 'k'); 
scatter(metrics.troughToPeak(isolatedCluster == 0), metrics.FWHM(isolatedCluster == 0), 70, muted_purple, 'filled', 'MarkerEdgeColor', 'k');
xlabel('Trough to peak duration (ms)','FontSize',12);
ylabel('FWHM (ms)','FontSize',12);
title('Cell Type Clusters','FontSize',14);
grid on;
legend('Interneuron', 'Principal Cell');
axis square
hold off;
% Principal Cell Cluster
subplot(2,4,2)
plot(meanWaveforms(isolatedCluster==0,:)', 'color', muted_purple)
axis tight square
xlabel('samples')
numObservationsPurple = sum(isolatedCluster==0);
title(sprintf('PC WF count: %s', num2str(numObservationsPurple)));
axis square
% Interneuron Cluster
subplot(2,4,3)
plot(meanWaveforms(isolatedCluster==1,:)', 'color', sage_green)
axis tight square
xlabel('samples')
numObservationsGreen = sum(isolatedCluster==1);
title(sprintf('IN WF count: %s', num2str(numObservationsGreen)));
axis square
% Subplot 1 for Trough-To-Peak (ms)
subplot(2,4,5);
% Create table for data
isolatedClusterMetric = double(isolatedCluster);
metricTable = table(isolatedClusterMetric, inputData(:,1), inputData(:,2),inputData(:,3), 'VariableNames',{'isolatedCluster', 'TTP', 'FWHM', 'ASYM'});
metricTable.isolatedCluster = grp2idx(metricTable.isolatedCluster);
% plot
hc = boxchart(metricTable.isolatedCluster, metricTable.TTP); % group by isolated Cluster
hc.BoxFaceColor = muted_purple;
%hc(2).BoxFaceColor = sage_green;
hold on
% overlay the scatter plots
for n=1:max(unique(metricTable.isolatedCluster))
    hs = scatter(ones(sum(metricTable.isolatedCluster==n),1) + n-1, metricTable.TTP(metricTable.isolatedCluster == n),"filled",'jitter','on','JitterAmount',0.1);
    hs.MarkerFaceAlpha = 0.7;
    hs.MarkerEdgeColor = 'k';
    hs.MarkerFaceColor = muted_purple;
    hs.MarkerFaceColor = sage_green;
end
xticks([1 2])
xticklabels(["Principal Cells", "Interneurons"])
ylabel('Trough-To-Peak (ms)', 'FontSize', 12);
ylim([0.0 max(metricTable.TTP)+0.1])
title('Trough-To-Peak', 'FontSize', 14);

axis square;
% Subplot 2 for FWHM (ms)
subplot(2,4,6);
% plot
hc = boxchart(metricTable.isolatedCluster, metricTable.FWHM); % group by isolated Cluster
hc.BoxFaceColor = muted_purple;
%hc(2).BoxFaceColor = sage_green;
hold on
% overlay the scatter plots
for n=1:max(unique(metricTable.isolatedCluster))
    hs = scatter(ones(sum(metricTable.isolatedCluster==n),1) + n-1, metricTable.FWHM(metricTable.isolatedCluster == n),"filled",'jitter','on','JitterAmount',0.1);
    hs.MarkerFaceAlpha = 0.8;
    hs.MarkerEdgeColor = 'k';
    hs.MarkerFaceColor = muted_purple;
    hs.MarkerFaceColor = sage_green;
end
xticks([1 2])
xticklabels(["Principal Cells", "Interneurons"])
ylabel('Full Width Half Max (ms)', 'FontSize', 12);
ylim([0.0 max(metricTable.FWHM)+0.1])
title('Full Width Half Max', 'FontSize', 14);
axis square;
% Subplot 3 for Asymmetry (if keeping ratio)
% subplot(2,4,7);
% % Plot the scatter plot
% scatter(find(~isolatedCluster), metrics.asymmetry(~isolatedCluster), 70, 'filled', 'MarkerFaceColor', muted_purple, 'MarkerEdgeColor', 'k');
% hold on;
% scatter(find(isolatedCluster), metrics.asymmetry(isolatedCluster), 70, 'filled', 'MarkerFaceColor', sage_green, 'MarkerEdgeColor', 'k');
% % Add a dashed black line at y-axis 1 across x-axis
% xlim([0 250])
% plot(xlim, ones(size(xlim)), 'k--', 'LineWidth', 1.5);
% % Add labels
% ylabel('WF Asymmetry (Ratio)', 'FontSize', 12);
% xlabel('Unit #')
% title('Asymmetry','FontSize',14);
% axis square
% grid on
% hold off;

% If switching Asymmetry to boxplot
subplot(2,4,7);
hc = boxchart(metricTable.isolatedCluster, metricTable.ASYM); % group by isolated Cluster
hc.BoxFaceColor = muted_purple;
%hc(2).BoxFaceColor = sage_green;
hold on
% overlay the scatter plots
for n=1:max(unique(metricTable.isolatedCluster))
    hs = scatter(ones(sum(metricTable.isolatedCluster==n),1) + n-1, metricTable.ASYM(metricTable.isolatedCluster == n),"filled",'jitter','on','JitterAmount',0.1);
    hs.MarkerFaceAlpha = 0.8;
    hs.MarkerEdgeColor = 'k';
    hs.MarkerFaceColor = muted_purple;
    hs.MarkerFaceColor = sage_green;
end
xticks([1 2])
xticklabels(["Principal Cells", "Interneurons"])
ylabel('Asymmetry (ms)', 'FontSize', 12);
ylim([0.0 max(metricTable.ASYM)+0.1])
title('Asymmetry', 'FontSize', 14);
axis square;


% sub-figure: plotting cell type and region responses.
AMYLocsIN = sum(AMYLocs(isolatedCluster))/sum(AMYLocs)*100;
AMYLocsPC = sum(AMYLocs(~isolatedCluster))/sum(AMYLocs)*100;
OFCLocsIN = sum(OFCLocs(isolatedCluster))/sum(OFCLocs)*100;
OFCLocsPC = sum(OFCLocs(~isolatedCluster))/sum(OFCLocs)*100;
HIPPLocsIN = sum(HIPPLocs(isolatedCluster))/sum(HIPPLocs)*100;
HIPPLocsPC = sum(HIPPLocs(~isolatedCluster))/sum(HIPPLocs)*100;
CINGLocsIN = sum(CINGLocs(isolatedCluster))/sum(CINGLocs)*100;
CINGLocsPC = sum(CINGLocs(~isolatedCluster))/sum(CINGLocs)*100;
pCINGLocsIN = sum(pCINGLocs(isolatedCluster))/sum(pCINGLocs)*100;
pCINGLocsPC = sum(pCINGLocs(~isolatedCluster))/sum(pCINGLocs)*100;

CellTypeRegions = [AMYLocsIN, AMYLocsPC; OFCLocsIN, OFCLocsPC; HIPPLocsIN, HIPPLocsPC; CINGLocsIN CINGLocsPC; pCINGLocsIN, pCINGLocsPC]; % dont need

% Create bar plot with specified colors
subplot(2,4,4)
b = bar(CellTypeRegions, 'stacked', 'FaceColor', sage_green, 'BarWidth', 1); % Transpose the data to have regions as columns
b(2).FaceColor = muted_purple;
b(1).BarWidth = .8;
b(2).BarWidth = .8;
% Add xtick labels
RegionLabels = {'Amygdala', 'OFC & vmPFC', 'Hippocampus', 'Cingulate', 'Posterior Cingulate'};
xticks(1:5)
xticklabels(RegionLabels)
ylim([0 100])
% Add labels and title
xlabel('Regions', 'FontSize', 12)
ylabel('Cell Type', 'FontSize', 12)
title('Cell Type Regions', 'FontSize', 14)
axis square

% sub-figure: raw baseline frequency scores.
% cell type baseline frequency
mean_BLfreq_IN_All = mean([SPESstruct(isolatedCluster).BLfreqAll], "omitnan");
std_BLfreq_IN_All = std([SPESstruct(isolatedCluster).BLfreqAll], "omitnan");
mean_BLfreq_PC_All = mean([SPESstruct(~isolatedCluster).BLfreqAll], "omitnan");
std_BLfreq_PC_All = std([SPESstruct(~isolatedCluster).BLfreqAll], "omitnan");

% statistics bewteen BL freq characteristics
[BLfreq_p,BLfreq_h, BLfreq_stats] = ranksum([SPESstruct(isolatedCluster).BLfreqAll],[SPESstruct(~isolatedCluster).BLfreqAll]); % BL freq

h1_IN = [SPESstruct(isolatedCluster).BLfreqAll];
h2_PC = [SPESstruct(~isolatedCluster).BLfreqAll];

subplot(2,4,8);
hold on
title('BL Freq IN/PC Responses')
histogram(h2_PC,'facealpha', 1,'FaceColor',muted_purple,'edgecolor','k', binWidth=2);
hold on;
histogram(h1_IN,'facealpha', 1,'FaceColor',sage_green,'edgecolor','k', binWidth=2);
ylabel('IN vs PC Units FR count')
xlabel('baseline window FR (Hz)')
text(40,12,sprintf('BLFreq:\nIN= %0.2fHz\nPC= %0.2fHz', mean_BLfreq_IN_All,mean_BLfreq_PC_All))
axis square
hold off
maximize(figure(10))
% save figure

saveas(10,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\WaveFormClassification\',sprintf('meanMetrics_IsolatedCluster_AcrossSubs.pdf')))
close(10)

% k2 Cluster
meanTTP_k2clusterRed = mean(inputData(idxTwo==1,1), 'omitnan'); % mean TTP PC
stdTTP_k2clusterRed = std(inputData(idxTwo==1,1), 'omitnan'); % std TTP PC
meanFWHM_k2clusterRed = mean(inputData(idxTwo==1,2), 'omitnan'); % mean FWHM PC
stdFWHM_k2clusterRed = std(inputData(idxTwo==1,2), 'omitnan'); % std FWHM PC
meanASYM_k2clusterRed = mean(inputData(idxTwo==1,3), 'omitnan'); % mean ASYM PC
stdASYM_k2clusterRed = std(inputData(idxTwo==1,3), 'omitnan'); % std ASYM PC
meanTTP_k2clusterBlue = mean(inputData(idxTwo==2,1), 'omitnan'); % mean TTP IN
stdTTP_k2clusterBlue = std(inputData(idxTwo==2,1), 'omitnan'); % std TTP IN
meanFWHM_k2clusterBlue = mean(inputData(idxTwo==2,2), 'omitnan'); % mean FWHM IN
stdFWHM_k2clusterBlue = std(inputData(idxTwo==2,2), 'omitnan'); % std FWHM IN
meanASYM_k2clusterBlue = mean(inputData(idxTwo==2,3), 'omitnan'); % mean ASYM IN
stdASYM_k2clusterBlue = std(inputData(idxTwo==2,3), 'omitnan'); % std ASYM IN

% mean metrics for k3 clusters
meanTTP_k3clusterOrange = mean(inputData(idxThree==1,1), 'omitnan');
stdTTP_k3clusterOrange = std(inputData(idxThree==1,1), 'omitnan');
meanFWHM_k3clusterOrange = mean(inputData(idxThree==1,2), 'omitnan');
stdFWHM_k3clusterOrange = std(inputData(idxThree==1,2), 'omitnan');
meanASYM_k3clusterOrange = mean(inputData(idxThree==1,3), 'omitnan');
stdASYM_k3clusterOrange = std(inputData(idxThree==1,3), 'omitnan');
meanTTP_k3clusterYellow = mean(inputData(idxThree==2,1), 'omitnan');
stdTTP_k3clusterYellow = std(inputData(idxThree==2,1), 'omitnan');
meanFWHM_k3clusterYellow = mean(inputData(idxThree==2,2), 'omitnan');
stdFWHM_k3clusterYellow = std(inputData(idxThree==2,2), 'omitnan');
meanASYM_k3clusterYellow = mean(inputData(idxThree==2,3), 'omitnan');
stdASYM_k3clusterYellow = std(inputData(idxThree==2,3), 'omitnan');
meanTTP_k3clusterTerra = mean(inputData(idxThree==3,1), 'omitnan');
stdTTP_k3clusterTerra = std(inputData(idxThree==3,1), 'omitnan');
meanFWHM_k3clusterTerra = mean(inputData(idxThree==3,2), 'omitnan');
stdFWHM_k3clusterTerra = std(inputData(idxThree==3,2), 'omitnan');
meanASYM_k3clusterTerra = mean(inputData(idxThree==3,3), 'omitnan');
stdASYM_k3clusterTerra = std(inputData(idxThree==3,3), 'omitnan');

meanMetricsK = [meanTTP_k3clusterOrange, meanFWHM_k3clusterOrange, meanASYM_k3clusterOrange; meanTTP_k3clusterYellow, meanFWHM_k3clusterYellow, meanASYM_k3clusterYellow; meanTTP_k3clusterTerra, meanFWHM_k3clusterTerra, meanASYM_k3clusterTerra];
sdMetricsK = [stdTTP_k3clusterOrange, stdFWHM_k3clusterOrange, stdASYM_k3clusterOrange; stdTTP_k3clusterYellow,stdFWHM_k3clusterYellow,stdASYM_k3clusterYellow; stdTTP_k3clusterTerra,stdFWHM_k3clusterTerra,stdASYM_k3clusterTerra]; 

 %% ~~~~~SUPPLEMENTARY FIGURE 2: WAVEFORM CLASSIFICATION K3 MEANS ~~~~~~

Burnt_Orange = [1 0.333 0];
Mustard_Yellow = [1 0.8 0];
Terracotta = [0.8 0.467 0.133];

figure(11)
% Cluster Scatter Plot
subplot(2,4,1);
hold on;
scatter(metrics.troughToPeak(idxThree == 1), metrics.FWHM(idxThree == 1), 70, Burnt_Orange, 'filled', 'MarkerEdgeColor', 'k'); % 
scatter(metrics.troughToPeak(idxThree == 2), metrics.FWHM(idxThree == 2), 70, Mustard_Yellow, 'filled', 'MarkerEdgeColor', 'k'); % 
scatter(metrics.troughToPeak(idxThree == 3), metrics.FWHM(idxThree == 3), 70, Terracotta, 'filled', 'MarkerEdgeColor', 'k'); % 
xlabel('Trough to peak duration (ms)','FontSize',12);
ylabel('FWHM (ms)','FontSize',12);
title('Three K-Clusters','FontSize',14);
grid on;
legend('Cluster 1', 'Cluster 2', 'Cluster 3');
axis square
hold off;
% K1 Cluster
subplot(2,4,2)
plot(meanWaveforms(idxThree==1,:)', 'color', Burnt_Orange)
axis tight square
xlabel('samples')
numObservationsk1 = sum(idxThree==1);
title(sprintf('K1 WF count: %s', num2str(numObservationsk1)));
axis square
% K2 Cluster
subplot(2,4,3)
plot(meanWaveforms(idxThree==2,:)', 'color', Mustard_Yellow)
axis tight square
xlabel('samples')
numObservationsk2 = sum(idxThree==2);
title(sprintf('K2 WF count: %s', num2str(numObservationsk2)));
axis square
% K3 Cluster
subplot(2,4,4)
plot(meanWaveforms(idxThree==3,:)', 'color', Terracotta)
axis tight square
xlabel('samples')
numObservationsk3 = sum(idxThree==3);
title(sprintf('K3 WF count: %s', num2str(numObservationsk3)));
axis square
% Subplot 1 for Trough-To-Peak (ms)
subplot(2,4,5);
% Create table for data
kClusterMetric = idxThree;
KmetricTable = table(kClusterMetric, inputData(:,1), inputData(:,2),inputData(:,3), 'VariableNames',{'kCluster', 'TTP', 'FWHM', 'ASYM'});
KmetricTable.kCluster = grp2idx(KmetricTable.kCluster);
% plot
hc = boxchart(KmetricTable.kCluster, KmetricTable.TTP); % group by kCluster
hc.BoxFaceColor = Burnt_Orange;
%hc(2).BoxFaceColor = Mustard_Yellow;
%hc(2).BoxFaceColor = Terracotta;
hold on
% overlay the scatter plots
for n=1:max(unique(KmetricTable.kCluster))
    hs = scatter(ones(sum(KmetricTable.kCluster==n),1) + n-1, KmetricTable.TTP(KmetricTable.kCluster == n),"filled",'jitter','on','JitterAmount',0.1);
    hs.MarkerFaceAlpha = 0.7;
    hs.MarkerEdgeColor = 'k';
    hs.MarkerFaceColor = Burnt_Orange;
    hs.MarkerFaceColor = Mustard_Yellow;
    hs.MarkerFaceColor = Terracotta;
end
xticks([1 2 3])
xticklabels(["K-Cluster 1", "K-Cluster 2", "K-Cluster 3"])
ylim([0.0 max(KmetricTable.TTP)+0.1])
ylabel('Trough-To-Peak (ms)', 'FontSize', 12);
title('Trough-To-Peak', 'FontSize', 14);
axis square;
% Subplot 2 for FWHM (ms)
subplot(2,4,6);
hc = boxchart(KmetricTable.kCluster, KmetricTable.FWHM); % group by kCluster
hc.BoxFaceColor = Burnt_Orange;
%hc(2).BoxFaceColor = Mustard_Yellow;
%hc(2).BoxFaceColor = Terracotta;
hold on
% overlay the scatter plots
for n=1:max(unique(KmetricTable.kCluster))
    hs = scatter(ones(sum(KmetricTable.kCluster==n),1) + n-1, KmetricTable.FWHM(KmetricTable.kCluster == n),"filled",'jitter','on','JitterAmount',0.1);
    hs.MarkerFaceAlpha = 0.7;
    hs.MarkerEdgeColor = 'k';
    hs.MarkerFaceColor = Burnt_Orange;
    hs.MarkerFaceColor = Mustard_Yellow;
    hs.MarkerFaceColor = Terracotta;
end
xticks([1 2 3])
xticklabels(["K-Cluster 1", "K-Cluster 2", "K-Cluster 3"])
ylim([0.0 max(KmetricTable.FWHM)+0.1])
ylabel('Full Width Half Max (ms)', 'FontSize', 12);
title('Full Width Half Max', 'FontSize', 14);
axis square;
% % Subplot 3 for Asymmetry
% subplot(2,4,7);
% % Plot the scatter plot
% scatter(find(idxThree == 1), metrics.asymmetry(idxThree == 1), 70, 'filled', 'MarkerFaceColor', Burnt_Orange, 'MarkerEdgeColor', 'k');
% hold on;
% scatter(find(idxThree == 2), metrics.asymmetry(idxThree == 2), 70, 'filled', 'MarkerFaceColor', Mustard_Yellow, 'MarkerEdgeColor', 'k');
% scatter(find(idxThree == 3), metrics.asymmetry(idxThree == 3), 70, 'filled', 'MarkerFaceColor', Terracotta, 'MarkerEdgeColor', 'k');
% % Add a dashed black line at y-axis 1 across x-axis
% xlim([0 250]);
% plot(xlim, ones(size(xlim)), 'k--', 'LineWidth', 1.5);
% % Add labels
% ylabel('WF Asymmetry (Ratio)', 'FontSize', 12);
% xlabel('Unit #')
% title('Asymmetry','FontSize',14);
% axis square
% grid on
% hold off;

% Asymmetry if boxplots
subplot(2,4,7);
hc = boxchart(KmetricTable.kCluster, KmetricTable.ASYM); % group by kCluster
hc.BoxFaceColor = Burnt_Orange;
%hc(2).BoxFaceColor = Mustard_Yellow;
%hc(2).BoxFaceColor = Terracotta;
hold on
% overlay the scatter plots
for n=1:max(unique(KmetricTable.kCluster))
    hs = scatter(ones(sum(KmetricTable.kCluster==n),1) + n-1, KmetricTable.ASYM(KmetricTable.kCluster == n),"filled",'jitter','on','JitterAmount',0.1);
    hs.MarkerFaceAlpha = 0.7;
    hs.MarkerEdgeColor = 'k';
    hs.MarkerFaceColor = Burnt_Orange;
    hs.MarkerFaceColor = Mustard_Yellow;
    hs.MarkerFaceColor = Terracotta;
end
xticks([1 2 3])
xticklabels(["K-Cluster 1", "K-Cluster 2", "K-Cluster 3"])
ylim([0.0 max(KmetricTable.ASYM)+0.1])
ylabel('Asymmetry (Ratio)', 'FontSize', 12);
title('Asymmetry', 'FontSize', 14);
axis square;

% sub-figure: raw baseline frequency scores.
% cell type baseline frequency
mean_BLfreq_k1_All = mean([SPESstruct(idxThree == 1).BLfreqAll], "omitnan");
std_BLfreq_k1_All = std([SPESstruct(idxThree == 1).BLfreqAll], "omitnan");
mean_BLfreq_k2_All = mean([SPESstruct(idxThree == 2).BLfreqAll], "omitnan");
std_BLfreq_k2_All = std([SPESstruct(idxThree == 2).BLfreqAll], "omitnan");
mean_BLfreq_k3_All = mean([SPESstruct(idxThree == 3).BLfreqAll], "omitnan");
std_BLfreq_k3_All = std([SPESstruct(idxThree == 3).BLfreqAll], "omitnan");

h1_k1 = [SPESstruct(idxThree == 1).BLfreqAll];
h2_k2 = [SPESstruct(idxThree == 2).BLfreqAll];
h3_k3 = [SPESstruct(idxThree == 3).BLfreqAll];

subplot(2,4,8);
hold on
title('BL Freq Cluster Responses')
histogram(h2_k2,'facealpha', 1,'FaceColor',Mustard_Yellow,'edgecolor','k', binWidth=2);
hold on;
histogram(h1_k1,'facealpha', 1,'FaceColor',Burnt_Orange,'edgecolor','k', binWidth=2);
histogram(h3_k3,'facealpha', 1,'FaceColor',Terracotta,'edgecolor','k', binWidth=2);
ylabel('Cluster Units FR count')
xlabel('baseline window FR (Hz)')
text(40,12,sprintf('BLFreq:\nK1= %0.2fHz\nK2= %0.2fHz\nK3= %0.2fHz', mean_BLfreq_k1_All,mean_BLfreq_k2_All,mean_BLfreq_k3_All))
axis square
hold off
maximize(figure(11))
% saving figure 
saveas(11,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\WaveFormClassification\',sprintf('meanMetrics_K3Cluster_AcrossSubs.pdf')))
close(11)

%%  ~~~~~~~~~ Manuscript Figure: SOZ units ~~~~~~~~~~~~ %%
nonSOZmean = mean(meanWaveforms(~isSOZ,:));
nonSOZn = sum(~isSOZ);
nonSOZstd = std(meanWaveforms(~isSOZ,:), 0, 1);
nonSOZse = nonSOZstd/sqrt(nonSOZn);

nonSOZse = std(meanWaveforms(~isSOZ))/ sqrt(length(meanWaveforms(~isSOZ)));
SOZmean = mean(meanWaveforms(isSOZ,:));
SOZse = std(meanWaveforms(isSOZ,:))/ sqrt(length(meanWaveforms(isSOZ,:)));

% Create patch variables
x_nonSOZvalues = 1:size(mean(meanWaveforms(~isSOZ,:)), 2); % x-coordinates for non-SOz
x_nonSOZpatch = [x_nonSOZvalues fliplr(x_nonSOZvalues)]; % create the x-coordinates for the patch
y_nonSOZpatch = [nonSOZmean + nonSOZse fliplr(nonSOZmean - nonSOZse)]; % create the y-coordinates for the patch
x_SOZvalues = 1:size(mean(meanWaveforms(isSOZ,:)), 2); % x-coordinates for SOz
x_SOZpatch = [x_SOZvalues fliplr(x_SOZvalues)]; % create the x-coordinates for the patch
y_SOZpatch = [SOZmean + SOZse fliplr(SOZmean - SOZse)]; % create the y-coordinates for the patch

figure(8)
% non-SOZ Cluster
hold on
patch(x_nonSOZpatch, y_nonSOZpatch, [0.5 0.5 0.5], 'facealpha', 0.5, 'edgecolor', 'none');
plot(mean(meanWaveforms(~isSOZ,:)), 'color', 'k','LineWidth',2.0)
patch(x_SOZpatch, y_SOZpatch, [0.5 0.5 0.5], 'facealpha', 0.5, 'edgecolor', 'none');
plot(mean(meanWaveforms(isSOZ,:)), 'color', 'yellow','LineWidth',2.0)
hold off
maximize(8)
axis tight square
xlabel('samples')
numObservationsnonSOZ = sum(~isSOZ);
numObservationsSOZ = sum(isSOZ);
title(sprintf('nonSOZ WF: %s, SOZ WF: %s', num2str(numObservationsnonSOZ), num2str(numObservationsSOZ)));
close(8)

figure(18)
% Cluster Scatter Plot
subplot(2,4,1);
hold on;
scatter(metrics.troughToPeak(isSOZ), metrics.FWHM(isSOZ), 70, 'yellow', 'filled', 'MarkerEdgeColor', 'k'); 
scatter(metrics.troughToPeak(~isSOZ), metrics.FWHM(~isSOZ), 70, 'k', 'filled', 'MarkerEdgeColor', 'k');
xlabel('Trough to peak duration (ms)','FontSize',12);
ylabel('FWHM (ms)','FontSize',12);
title('SOZ Clusters','FontSize',14);
grid on;
legend('SOZ', 'Other');
axis square
hold off;
% non-SOZ Cluster
subplot(2,4,2)
plot(meanWaveforms(~isSOZ,:)', 'color', 'k')
axis tight square
xlabel('samples')
numObservationsnonSOZ = sum(~isSOZ);
title(sprintf('nonSOZ WF count: %s', num2str(numObservationsnonSOZ)));
axis square
% SOZ Cluster
subplot(2,4,3)
plot(mean(meanWaveforms(isSOZ,:))', 'color', 'yellow')
axis tight square
xlabel('samples')
numObservationsSOZ = sum(isSOZ);
title(sprintf('SOZ WF count: %s', num2str(numObservationsSOZ)));
axis square
% Subplot 1 for Trough-To-Peak (ms)
subplot(2,4,5);
% Create table for data
SOZMetric = double(isSOZ);
SOZmetricTable = table(SOZMetric', inputData(:,1), inputData(:,2),inputData(:,3), 'VariableNames',{'isSOZ', 'TTP', 'FWHM', 'ASYM'});
SOZmetricTable.isSOZ = grp2idx(SOZmetricTable.isSOZ);
% plot
hc = boxchart(SOZmetricTable.isSOZ, SOZmetricTable.TTP); % group by isolated Cluster
hc.BoxFaceColor = 'k';
%hc(2).BoxFaceColor = 'yellow';
hold on
% overlay the scatter plots
for n=1:max(unique(SOZmetricTable.isSOZ))
    hs = scatter(ones(sum(SOZmetricTable.isSOZ==n),1) + n-1, SOZmetricTable.TTP(SOZmetricTable.isSOZ == n),"filled",'jitter','on','JitterAmount',0.1);
    hs.MarkerFaceAlpha = 0.7;
    hs.MarkerEdgeColor = 'k';
    hs.MarkerFaceColor = 'k';
    hs.MarkerFaceColor = 'yellow';
end
xticks([1 2])
xticklabels(["other", "SOZ"])
ylabel('Trough-To-Peak (ms)', 'FontSize', 12);
ylim([0.0 max(SOZmetricTable.TTP)+0.1])
title('Trough-To-Peak', 'FontSize', 14);
axis square;
% Subplot 2 for FWHM (ms)
subplot(2,4,6);
% plot
hc = boxchart(SOZmetricTable.isSOZ, SOZmetricTable.FWHM); % group by isolated Cluster
hc.BoxFaceColor = 'k';
%hc(2).BoxFaceColor = 'yellow';
hold on
% overlay the scatter plots
for n=1:max(unique(SOZmetricTable.isSOZ))
    hs = scatter(ones(sum(SOZmetricTable.isSOZ==n),1) + n-1, SOZmetricTable.FWHM(SOZmetricTable.isSOZ == n),"filled",'jitter','on','JitterAmount',0.1);
    hs.MarkerFaceAlpha = 0.8;
    hs.MarkerEdgeColor = 'k';
    hs.MarkerFaceColor = 'k';
    hs.MarkerFaceColor = 'yellow';
end
xticks([1 2])
xticklabels(["other", "SOZ"])
ylabel('Full Width Half Max (ms)', 'FontSize', 12);
ylim([0.0 max(SOZmetricTable.FWHM)+0.1])
title('Full Width Half Max', 'FontSize', 14);
axis square;
% % Subplot 3 for Asymmetry
% subplot(2,4,7);
% % Plot the scatter plot
% scatter(find(~isSOZ), metrics.asymmetry(~isSOZ), 70, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
% hold on;
% scatter(find(isSOZ), metrics.asymmetry(isSOZ), 70, 'filled', 'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'k');
% % Add a dashed black line at y-axis 1 across x-axis
% xlim([0 250]);
% plot(xlim, ones(size(xlim)), 'k--', 'LineWidth', 1.5);
% % Add labels
% ylabel('WF Asymmetry (Ratio)', 'FontSize', 12);
% xlabel('Unit #')
% title('Asymmetry','FontSize',14);
% axis square
% grid on
% hold off;

% If doing asymmetry as boxplots.
subplot(2,4,7);
% plot
hc = boxchart(SOZmetricTable.isSOZ, SOZmetricTable.ASYM); % group by isolated Cluster
hc.BoxFaceColor = 'k';
%hc(2).BoxFaceColor = 'yellow';
hold on
% overlay the scatter plots
for n=1:max(unique(SOZmetricTable.isSOZ))
    hs = scatter(ones(sum(SOZmetricTable.isSOZ==n),1) + n-1, SOZmetricTable.ASYM(SOZmetricTable.isSOZ == n),"filled",'jitter','on','JitterAmount',0.1);
    hs.MarkerFaceAlpha = 0.8;
    hs.MarkerEdgeColor = 'k';
    hs.MarkerFaceColor = 'k';
    hs.MarkerFaceColor = 'yellow';
end
xticks([1 2])
xticklabels(["other", "SOZ"])
ylabel('Asymmetry (Ratio)', 'FontSize', 12);
ylim([0.0 max(SOZmetricTable.ASYM)+0.1])
title('Asymmetry', 'FontSize', 14);
axis square;

% sub-figure: plotting cell type and region responses.
AMYLocsSOZ = sum(AMYLocs(isSOZ))/sum(AMYLocs)*100
AMYLocsNonSOZ = sum(AMYLocs(~isSOZ))/sum(AMYLocs)*100
OFCLocsSOZ = sum(OFCLocs(isSOZ))/sum(OFCLocs)*100
OFCLocsNonSOZ = sum(OFCLocs(~isSOZ))/sum(OFCLocs)*100
HIPPLocsSOZ = sum(HIPPLocs(isSOZ))/sum(HIPPLocs)*100
HIPPLocsNonSOZ = sum(HIPPLocs(~isSOZ))/sum(HIPPLocs)*100
CINGLocsSOZ = sum(CINGLocs(isSOZ))/sum(CINGLocs)*100
CINGLocsNonSOZ = sum(CINGLocs(~isSOZ))/sum(CINGLocs)*100
pCINGLocsSOZ = sum(pCINGLocs(isSOZ))/sum(pCINGLocs)*100
pCINGLocsNonSOZ = sum(pCINGLocs(~isSOZ))/sum(pCINGLocs)*100

SOZRegions = [AMYLocsSOZ, AMYLocsNonSOZ; OFCLocsSOZ, OFCLocsNonSOZ; HIPPLocsSOZ, HIPPLocsNonSOZ; CINGLocsSOZ, CINGLocsNonSOZ; pCINGLocsSOZ, pCINGLocsNonSOZ]; % dont need

% Create bar plot with specified colors
subplot(2,4,4)
b = bar(SOZRegions, 'stacked', 'FaceColor', 'yellow', 'BarWidth', 1); % Transpose the data to have regions as columns
b(2).FaceColor = 'k';
b(1).BarWidth = .8;
b(2).BarWidth = .8;
% Add xtick labels
RegionLabels = {'Amygdala', 'OFC & vmPFC', 'Hippocampus', 'Cingulate', 'Posterior Cingulate'};
xticks(1:5)
xticklabels(RegionLabels)
ylim([0 100])
% Add labels and title
xlabel('Regions', 'FontSize', 12)
ylabel('SOZ', 'FontSize', 12)
title('SOZ Regions', 'FontSize', 14)
axis square

% sub-figure: raw baseline frequency scores.
% cell type baseline frequency
mean_BLfreq_SOZ_All = mean([SPESstruct(isSOZ).BLfreqAll], "omitnan");
std_BLfreq_SOZ_All = std([SPESstruct(isSOZ).BLfreqAll], "omitnan");
mean_BLfreq_nonSOZ_All = mean([SPESstruct(~isSOZ).BLfreqAll], "omitnan");
std_BLfreq_nonSOZ_All = std([SPESstruct(~isSOZ).BLfreqAll], "omitnan");

% statistics bewteen BL freq characteristics
[BLfreq_SOZh,BLfreq_SOZp,BLfreq_SOZci, BLfreq_SOZstats] = ttest2([SPESstruct(isSOZ).BLfreqAll],[SPESstruct(~isSOZ).BLfreqAll],'Tail','both', 'varType', 'unequal'); % BL freq
h1_SOZ = [SPESstruct(isSOZ).BLfreqAll];
h2_nonSOZ = [SPESstruct(~isSOZ).BLfreqAll];

subplot(2,4,8);
hold on
title('BL Freq SOZ/nonSOZ Responses')
histogram(h2_nonSOZ,'facealpha', 1,'FaceColor','k','edgecolor','k', binWidth=2);
hold on;
histogram(h1_SOZ,'facealpha', 1,'FaceColor','yellow','edgecolor','k', binWidth=2);
ylabel('SOZ vs other Units FR count')
xlabel('baseline window FR (Hz)')
text(40,12,sprintf('BLFreq:\nSOZ= %0.2fHz\nother= %0.2fHz', mean_BLfreq_SOZ_All,mean_BLfreq_nonSOZ_All))
axis square
hold off
maximize(figure(18))
% save figure
saveas(18,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\WaveFormClassification\',sprintf('meanMetrics_SOZ_AcrossSubs.pdf')))
close(18)

%%~~~~~~~~~ Neuronal Feature Characterization Metrics ~~~~~~~~~~~~~~~~~  %%
% We will only analyze  suppression units to look at suppression amplitude etc. 

sigUnits = [SPESstruct(:).hAll] == 1; % significantly modulated units = 174
sigINs = sigUnits(isolatedCluster); % number of significant INs = 33
sigPCs = sigUnits(~isolatedCluster);  % number of significant PCs = 141 
sigSOZ = isSOZ(sigUnits);  % number of significant SOZs = 20

% logical index for significant cell types and SOZ
isolatedClusterSig = isolatedCluster(sigUnits) % all significant units, 1 = INs, 0 = PCs
isolatedClusterSOZ = isolatedCluster(isSOZ) % all SOZ units, 1 = INs, 0 = PCs

sigSPESstruct = [SPESstruct(sigUnits)]; % Dropping units that werent modulated.
dataSigUnits = [sigSPESstruct(:).RdiffAll];
suppressionSigUnits = dataSigUnits < 0; % how many units showed suppression characteristics? = 159
sigsupp_SPESstruct = [sigSPESstruct(suppressionSigUnits)]; % Sig Units (174)
excitationSigUnits = dataSigUnits > 0; % how many units showed excitation characteristics? = 15.
sigsexci_SPESstruct = [sigSPESstruct(excitationSigUnits)]; % Sig Units (15)

p1 = 159 / 174;  % 0.9138
p2 = 0.6;
cohens_h = 2 * (asin(sqrt(p1)) - asin(sqrt(p2)));


suppressedVals = dataSigUnits(suppressionSigUnits);
nonSuppressedVals = dataSigUnits(~suppressionSigUnits);

mean1 = mean(suppressedVals);
mean2 = mean(nonSuppressedVals);
stdPooled = sqrt((var(suppressedVals) + var(nonSuppressedVals)) / 2);

cohens_d = (mean1 - mean2) / stdPooled;


% isolatedCluster and Suppressed logical index
isolatedClusterSigSupp = isolatedClusterSig(suppressionSigUnits);  % all significant suppressed units, 1 = INs, 0 = PCs
isolatedClusterSigExci = isolatedClusterSig(excitationSigUnits);
% SOZ and suppressed logical index
SOZsigSupp = sigSOZ(suppressionSigUnits); % all SOZ that are suppressed 
SOZsigExci = sigSOZ(excitationSigUnits);
% INs non and SOZ, and suppressed
isolatedClusterSOZSigSupp = isolatedClusterSigSupp(SOZsigSupp); % SOZ units that are sig suppressed, 1 = SOZINs
isolatedClusternonSOZSigSupp = isolatedClusterSigSupp(~SOZsigSupp); % SOZ units that are sig suppressed, 1 = nonSOZINs
% PCs non and SOZ, and suppressed
nonisolatedClusterSOZSigSupp = isolatedClusterSigSupp(SOZsigSupp) == 0; % SOZ units that are sig suppressed, 1 = SOZPCs
nonisolatedClusternonSOZSigSupp = isolatedClusterSigSupp(~SOZsigSupp) ==0; % SOZ units that are sig suppressed, 1 = nonSOZPCs

 %% Proportion testing: between cell types and SOZs

 % (1) No differnece in the # of INs and PCs that were suppressed/excited
 INprop_supp = sum(isolatedClusterSigSupp); % (31/33)
 PCprop_supp = sum(~isolatedClusterSigSupp); % (128/141)
 INprop_exci = sum(isolatedClusterSigExci); % (2/33)
 PCprop_exci = sum(~isolatedClusterSigExci); % (13/141)
 totals = [sum(isolatedClusterSig) sum(~isolatedClusterSig)];
 [hsig_supp,psig_supp,chisig_supp,dgsig_supp] = prop_test([INprop_supp PCprop_supp], totals,false); % not significant differences in proportions on SOZ INs and PCs.
 [hsig_exci,psig_exci,chisig_exci,dgsig_exci] = prop_test([INprop_exci PCprop_exci], totals,false); % not significant differences in proportions on SOZ INs and PCs.

% (2) No difference in SOZ or nonSOZ suppressed
 SOZprop_supp = sum(SOZsigSupp); % (5)
 nonSOZprop_supp = sum(~SOZsigSupp); %  (13)
 SOZtotals = [sum(sigSOZ) sum(~sigSOZ)]; % incorrect.. what is the total of interest here?
 [hsig_SOZsupp,psig_SOZsupp,chisig_SOZsupp,dgsig_SOZsupp] = prop_test([SOZprop_supp nonSOZprop_supp], SOZtotals,false); % not significant differences in proportions on SOZ INs and PCs.

 % (3) No differnece in the # of SOZ INs and SOZ PCs that were suppressed/excited
 INSOZprop_supp = sum(isolatedClusterSOZSigSupp); % (5)
 PCSOZprop_supp = sum(nonisolatedClusterSOZSigSupp); %  (13)
 SOZtotals = [sum(isolatedClusterSigSupp) sum(~isolatedClusterSigSupp)]; % incorrect.. what is the total of interest here?
 [hsig_celltypeSOZsupp,psig_celltypeSOZsupp,chisig_celltypeSOZsupp,dgsig_celltypeSOZsupp] = prop_test([INSOZprop_supp PCSOZprop_supp], SOZtotals,false); % not significant differences in proportions on SOZ INs and PCs.

%% ~~~~~~~ 1A. All Suppressed: Average firing rate responses for suppressed units ~~~~~~~~~ %% 

% sig supp units: baseline frequency (-1.1 - -0.1s) (raw)
mean_BLfreq_SigSupp = mean([sigsupp_SPESstruct(:).BLfreqAll], "omitnan");
std_BLfreq_SigSupp = std([sigsupp_SPESstruct(:).BLfreqAll], "omitnan");
%sig supp units: post-stim frequency  (0.1 - 1.1s) (raw)
mean_PSfreq_SigSupp = mean([sigsupp_SPESstruct(:).PSfreqAll], "omitnan");
std_PSfreq_SigSupp = std([sigsupp_SPESstruct(:).PSfreqAll], "omitnan");
% sig supp units: FR Hz change pre to post (raw)
mean_FR_SigSupp = mean([sigsupp_SPESstruct(:).RdiffAll], "omitnan");
std_FR_SigSupp = std([sigsupp_SPESstruct(:).RdiffAll], "omitnan");
% sig supp units: significant suppression amplitude (raw)
mean_suppAmp_SigSupp = mean([sigsupp_SPESstruct(:).Minimum_SuppressionAll], "omitnan");
std_suppAmp_SigSupp = std([sigsupp_SPESstruct(:).Minimum_SuppressionAll], "omitnan");
% sig supp units: significant suppression latency
mean_suppLatency_SigSupp = mean([sigsupp_SPESstruct(:).suppLatencyAll], "omitnan");
std_suppLatency_SigSupp = std([sigsupp_SPESstruct(:).suppLatencyAll], "omitnan");
% sig supp units: significant suppression duration
mean_suppDuration_SigSupp = mean([sigsupp_SPESstruct(:).suppressionDurationTimeAll], "omitnan");
std_suppDuration_SigSupp = std([sigsupp_SPESstruct(:).suppressionDurationTimeAll], "omitnan");
% sig supp units: significant suppression threshold
mean_suppThreshold_SigSupp = mean([sigsupp_SPESstruct(:).durationThresholdAll], "omitnan");
std_suppThreshold_SigSupp = std([sigsupp_SPESstruct(:).durationThresholdAll], "omitnan");

% sig supp units: significant suppression threshold
mean_recoveryRate_SigSupp = mean([sigsupp_SPESstruct(:).recoveryRate_Time], "omitnan");
std_recoveryRate_SigSupp = std([sigsupp_SPESstruct(:).recoveryRate_Time], "omitnan");

mean_recoveryRate_SigSupp = mean([sigsupp_SPESstruct(isolatedClusterSigSupp).recoveryRate_Time], "omitnan");
std_recoveryRate_SigSupp = std([sigsupp_SPESstruct(isolatedClusterSigSupp).recoveryRate_Time], "omitnan");

% changing empty cells from supp duration to zeros...
% If you want to assign it back to the structure
for i = 1:length(sigsupp_SPESstruct)
    if isempty(sigsupp_SPESstruct(i).suppressionDurationTimeAll)
        sigsupp_SPESstruct(i).suppressionDurationTimeAll = 0;
    end
end

 

% scatter (all: amplitude vs latency)
scatter([sigsupp_SPESstruct(:).Minimum_SuppressionAll], [sigsupp_SPESstruct(:).suppLatencyAll])
mdl = fitlm([sigsupp_SPESstruct(:).Minimum_SuppressionAll], [sigsupp_SPESstruct(:).suppLatencyAll])
plot(mdl)
% scatter (all: amplitude vs duration) (SIG).
scatter([sigsupp_SPESstruct(:).Minimum_SuppressionAll], [sigsupp_SPESstruct(:).suppressionDurationTimeAll])
mdl = fitlm([sigsupp_SPESstruct(:).Minimum_SuppressionAll], [sigsupp_SPESstruct(:).suppressionDurationTimeAll], 'Exclude', [sigsupp_SPESstruct(:).suppressionDurationTimeAll] == 0)
plot(mdl)
% scatter (all: duration vs latency)
scatter([sigsupp_SPESstruct(:).suppLatencyAll], [sigsupp_SPESstruct(:).suppressionDurationTimeAll])
mdl = fitlm([sigsupp_SPESstruct(:).suppLatencyAll], [sigsupp_SPESstruct(:).suppressionDurationTimeAll], 'Exclude', [sigsupp_SPESstruct(:).suppressionDurationTimeAll] == 0)
plot(mdl)

% scatter (INs: amplitude vs duration)
mdl = fitlm([sigsupp_SPESstruct(isolatedClusterSigSupp).Minimum_SuppressionAll], [sigsupp_SPESstruct(isolatedClusterSigSupp).suppressionDurationTimeAll], 'Exclude', [sigsupp_SPESstruct(isolatedClusterSigSupp).suppressionDurationTimeAll] == 0)
% scatter (INs: amplitude vs latency)
mdl = fitlm([sigsupp_SPESstruct(isolatedClusterSigSupp).Minimum_SuppressionAll], [sigsupp_SPESstruct(isolatedClusterSigSupp).suppLatencyAll])
% scatter (INs: duration vs latency)
mdl = fitlm([sigsupp_SPESstruct(isolatedClusterSigSupp).suppLatencyAll], [sigsupp_SPESstruct(isolatedClusterSigSupp).suppressionDurationTimeAll], 'Exclude', [sigsupp_SPESstruct(isolatedClusterSigSupp).suppressionDurationTimeAll] == 0)
% scatter (PCs: amplitude vs duration) SIG
mdl = fitlm([sigsupp_SPESstruct(~isolatedClusterSigSupp).Minimum_SuppressionAll], [sigsupp_SPESstruct(~isolatedClusterSigSupp).suppressionDurationTimeAll], 'Exclude', [sigsupp_SPESstruct(~isolatedClusterSigSupp).suppressionDurationTimeAll] == 0)
% scatter (PCs: amplitude vs latency)
mdl = fitlm([sigsupp_SPESstruct(~isolatedClusterSigSupp).Minimum_SuppressionAll], [sigsupp_SPESstruct(~isolatedClusterSigSupp).suppLatencyAll])
% scatter (PCs: duration vs latency)
mdl = fitlm([sigsupp_SPESstruct(~isolatedClusterSigSupp).suppLatencyAll], [sigsupp_SPESstruct(~isolatedClusterSigSupp).suppressionDurationTimeAll], 'Exclude', [sigsupp_SPESstruct(~isolatedClusterSigSupp).suppressionDurationTimeAll] == 0)
% scatter (SOZs: amplitude vs duration)
mdl = fitlm([sigsupp_SPESstruct(SOZsigSupp).Minimum_SuppressionAll], [sigsupp_SPESstruct(SOZsigSupp).suppressionDurationTimeAll], 'Exclude', [sigsupp_SPESstruct(SOZsigSupp).suppressionDurationTimeAll] == 0)
% scatter (SOZs: amplitude vs latency)
mdl = fitlm([sigsupp_SPESstruct(SOZsigSupp).Minimum_SuppressionAll], [sigsupp_SPESstruct(SOZsigSupp).suppLatencyAll])
% scatter (SOZs: duration vs latency)
mdl = fitlm([sigsupp_SPESstruct(SOZsigSupp).suppLatencyAll], [sigsupp_SPESstruct(SOZsigSupp).suppressionDurationTimeAll], 'Exclude', [sigsupp_SPESstruct(SOZsigSupp).suppressionDurationTimeAll] == 0)
% scatter (nonSOZs: amplitude vs duration) SIG
mdl = fitlm([sigsupp_SPESstruct(~SOZsigSupp).Minimum_SuppressionAll], [sigsupp_SPESstruct(~SOZsigSupp).suppressionDurationTimeAll], 'Exclude', [sigsupp_SPESstruct(~SOZsigSupp).suppressionDurationTimeAll] == 0)
% scatter (nonSOZs: amplitude vs latency) SIG
mdl = fitlm([sigsupp_SPESstruct(~SOZsigSupp).Minimum_SuppressionAll], [sigsupp_SPESstruct(~SOZsigSupp).suppLatencyAll])
% scatter (nonSOZs: duration vs latency) NS
mdl = fitlm([sigsupp_SPESstruct(~SOZsigSupp).suppLatencyAll], [sigsupp_SPESstruct(~SOZsigSupp).suppressionDurationTimeAll], 'Exclude', [sigsupp_SPESstruct(~SOZsigSupp).suppressionDurationTimeAll] == 0)


%mdl = fitlm([sigsupp_SPESstruct(~SOZsigSupp).recoveryRate_Time], [sigsupp_SPESstruct(~SOZsigSupp).suppressionDurationTimeAll], 'Exclude', [sigsupp_SPESstruct(~SOZsigSupp).suppressionDurationTimeAll] == 0)


%% Supplementary Regression Figure
% Scatter Plot: Cell Type: Suppression Duration vs. Suppression Latency 
% Create figure and scatter plots
figure(123);
subplot(2,3,1);
hold on;
% Extract data for IN
temp_y1 = [sigsupp_SPESstruct(isolatedClusterSigSupp).suppLatencyAll];
temp_x1 = [sigsupp_SPESstruct(isolatedClusterSigSupp).suppressionDurationTimeAll];
y1 = temp_y1(temp_x1 ~= 0);
x1 = temp_x1(temp_x1 ~= 0);
% Extract data for PC
temp_y2 = [sigsupp_SPESstruct(~isolatedClusterSigSupp).suppLatencyAll];
temp_x2 = [sigsupp_SPESstruct(~isolatedClusterSigSupp).suppressionDurationTimeAll];
y2 = temp_y2(temp_x2 ~= 0);
x2 = temp_x2(temp_x2 ~= 0);
scatter(x1, y1, 80, sage_green, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
scatter(x2, y2, 80, muted_purple, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
% Fit linear regression models
mdl1 = fitlm([sigsupp_SPESstruct(isolatedClusterSigSupp).suppressionDurationTimeAll], [sigsupp_SPESstruct(isolatedClusterSigSupp).suppLatencyAll], 'Exclude', [sigsupp_SPESstruct(isolatedClusterSigSupp).suppressionDurationTimeAll] == 0); % Regression model for sage green data
mdl2 = fitlm([sigsupp_SPESstruct(~isolatedClusterSigSupp).suppressionDurationTimeAll], [sigsupp_SPESstruct(~isolatedClusterSigSupp).suppLatencyAll], 'Exclude', [sigsupp_SPESstruct(~isolatedClusterSigSupp).suppressionDurationTimeAll] == 0); % Regression model for muted purple data
% Plot regression lines
xRange1 = linspace(min(x1), max(x1), 100); % Generate x values for plotting
yPred1 = predict(mdl1, xRange1'); % Predict y values from the model
plot(xRange1, yPred1, 'Color', sage_green, 'LineStyle', '-', 'LineWidth', 2);
xRange2 = linspace(min(x2), max(x2), 100); % Generate x values for plotting
yPred2 = predict(mdl2, xRange2'); % Predict y values from the model
plot(xRange2, yPred2, 'Color', muted_purple, 'LineStyle', '-', 'LineWidth', 2);
text(2, 1.1, sprintf('IN pVal = %.4f', mdl1.ModelFitVsNullModel.Pvalue))
text(2, 1, sprintf('PC pVal = %.4f', mdl2.ModelFitVsNullModel.Pvalue))
% axis properties
ylabel('Suppression Latency (s)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Suppression Duration (s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Suppression Duration vs. Suppression Latency', ...
      'FontSize', 14, 'FontWeight', 'bold');
% Add grid lines for better readability
grid on;
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'Box', 'on', 'LineWidth', 1.5);
axis square
hold off;
% Maximize figure window
set(gcf, 'Units', 'normalized', 'Position', [0 0 1 1]);

% Scatter Plot: Cell Type: Suppression Amplitude vs. Suppression Duration
% Create figure and scatter plots
subplot(2,3,2);
hold on;
% Extract data for IN
temp_y3 = [sigsupp_SPESstruct(isolatedClusterSigSupp).Minimum_SuppressionAll];
temp_x3 = [sigsupp_SPESstruct(isolatedClusterSigSupp).suppressionDurationTimeAll];
y3 = temp_y3(temp_x3 ~= 0);
x3 = temp_x3(temp_x3 ~= 0);
% Extract data for PC
temp_y4 = [sigsupp_SPESstruct(~isolatedClusterSigSupp).Minimum_SuppressionAll];
temp_x4 = [sigsupp_SPESstruct(~isolatedClusterSigSupp).suppressionDurationTimeAll];
y4 = temp_y4(temp_x4 ~= 0);
x4 = temp_x4(temp_x4 ~= 0);
scatter(x3, y3, 80, sage_green, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
scatter(x4, y4, 80, muted_purple, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
% Fit linear regression models
mdl3 = fitlm([sigsupp_SPESstruct(isolatedClusterSigSupp).suppressionDurationTimeAll], [sigsupp_SPESstruct(isolatedClusterSigSupp).Minimum_SuppressionAll], 'Exclude', [sigsupp_SPESstruct(isolatedClusterSigSupp).suppressionDurationTimeAll] == 0); % Regression model for sage green data
mdl4 = fitlm([sigsupp_SPESstruct(~isolatedClusterSigSupp).suppressionDurationTimeAll], [sigsupp_SPESstruct(~isolatedClusterSigSupp).Minimum_SuppressionAll], 'Exclude', [sigsupp_SPESstruct(~isolatedClusterSigSupp).suppressionDurationTimeAll] == 0); % Regression model for muted purple data
% Plot regression lines
xRange2 = linspace(min(x3), max(x3), 100); % Generate x values for plotting
yPred2 = predict(mdl3, xRange2'); % Predict y values from the model
plot(xRange2, yPred2, 'Color', sage_green, 'LineStyle', '-', 'LineWidth', 2);
xRange3 = linspace(min(x4), max(x4), 100); % Generate x values for plotting
yPred3 = predict(mdl4, xRange3'); % Predict y values from the model
plot(xRange3, yPred3, 'Color', muted_purple, 'LineStyle', '-', 'LineWidth', 2);
text(2, 60, sprintf('IN pVal = %.4f', mdl3.ModelFitVsNullModel.Pvalue))
text(2, 55, sprintf('PC pVal = %.4f', mdl4.ModelFitVsNullModel.Pvalue))
% axis properties
ylabel('Suppression Amplitude (spks/s)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Suppression Duration (s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Suppression Amplitude vs. Suppression Duration', ...
      'FontSize', 14, 'FontWeight', 'bold');
% Add grid lines for better readability
grid on;
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'Box', 'on', 'LineWidth', 1.5);
axis square
hold off;
% Maximize figure window
set(gcf, 'Units', 'normalized', 'Position', [0 0 1 1]);

% Scatter Plot: Cell Type: Suppression Amplitude vs. Suppression Latency
% Create figure and scatter plots
subplot(2,3,3);
hold on;
% Extract data for IN
temp_y9 = [sigsupp_SPESstruct(isolatedClusterSigSupp).Minimum_SuppressionAll];
temp_x9 = [sigsupp_SPESstruct(isolatedClusterSigSupp).suppLatencyAll];
y9 = temp_y9(temp_x9 ~= 0);
x9 = temp_x9(temp_x9 ~= 0);
% Extract data for PC
temp_y10 = [sigsupp_SPESstruct(~isolatedClusterSigSupp).Minimum_SuppressionAll];
temp_x10 = [sigsupp_SPESstruct(~isolatedClusterSigSupp).suppLatencyAll];
y10 = temp_y10(temp_x10 ~= 0);
x10 = temp_x10(temp_x10 ~= 0);
scatter(x9, y9, 80, sage_green, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
scatter(x10, y10, 80, muted_purple, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
% Fit linear regression models
mdl9 = fitlm([sigsupp_SPESstruct(isolatedClusterSigSupp).suppLatencyAll], [sigsupp_SPESstruct(isolatedClusterSigSupp).Minimum_SuppressionAll]); % Regression model for sage green data
mdl10 = fitlm([sigsupp_SPESstruct(~isolatedClusterSigSupp).suppLatencyAll], [sigsupp_SPESstruct(~isolatedClusterSigSupp).Minimum_SuppressionAll]); % Regression model for muted purple data
% Plot regression lines
xRange2 = linspace(min(x9), max(x9), 100); % Generate x values for plotting
yPred2 = predict(mdl9, xRange2'); % Predict y values from the model
plot(xRange2, yPred2, 'Color', sage_green, 'LineStyle', '-', 'LineWidth', 2);
xRange3 = linspace(min(x10), max(x10), 100); % Generate x values for plotting
yPred3 = predict(mdl10, xRange3'); % Predict y values from the model
plot(xRange3, yPred3, 'Color', muted_purple, 'LineStyle', '-', 'LineWidth', 2);
text(2, 60, sprintf('IN pVal = %.4f', mdl9.ModelFitVsNullModel.Pvalue))
text(2, 55, sprintf('PC pVal = %.4f', mdl10.ModelFitVsNullModel.Pvalue))
% axis properties
ylabel('Suppression Amplitude (spks/s)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Suppression Latency (s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Suppression Amplitude vs. Suppression Duration', ...
      'FontSize', 14, 'FontWeight', 'bold');
% Add grid lines for better readability
grid on;
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'Box', 'on', 'LineWidth', 1.5);
axis square
hold off;
% Maximize figure window
set(gcf, 'Units', 'normalized', 'Position', [0 0 1 1]);

% Scatter Plot: SOZ Suppression Amplitude vs. Suppression Latency
% Create figure and scatter plots
subplot(2,3,4);
hold on;
% Extract data for IN
temp_y7 = [sigsupp_SPESstruct(SOZsigSupp).suppLatencyAll];
temp_x7 = [sigsupp_SPESstruct(SOZsigSupp).suppressionDurationTimeAll];
y7 = temp_y7(temp_x7 ~= 0);
x7 = temp_x7(temp_x7 ~= 0);
% Extract data for PC
temp_y8 = [sigsupp_SPESstruct(~SOZsigSupp).suppLatencyAll];
temp_x8 = [sigsupp_SPESstruct(~SOZsigSupp).suppressionDurationTimeAll];
y8 = temp_y8(temp_x8 ~= 0);
x8 = temp_x8(temp_x8 ~= 0);
scatter(x8, y8, 80, 'k', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
scatter(x7, y7, 80, 'y', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5); % plot SOZ second so we can see it
% Fit linear regression models
mdl7 = fitlm([sigsupp_SPESstruct(SOZsigSupp).suppressionDurationTimeAll], [sigsupp_SPESstruct(SOZsigSupp).suppLatencyAll], 'Exclude', [sigsupp_SPESstruct(SOZsigSupp).suppressionDurationTimeAll] == 0); % Regression model for sage green data
mdl8 = fitlm([sigsupp_SPESstruct(~SOZsigSupp).suppressionDurationTimeAll], [sigsupp_SPESstruct(~SOZsigSupp).suppLatencyAll], 'Exclude', [sigsupp_SPESstruct(~SOZsigSupp).suppressionDurationTimeAll] == 0); % Regression model for muted purple data
% Plot regression lines
xRange4 = linspace(min(x7), max(x7), 100); % Generate x values for plotting
yPred4 = predict(mdl7, xRange4'); % Predict y values from the model
plot(xRange4, yPred4, 'Color', 'y', 'LineStyle', '-', 'LineWidth', 2);
xRange5 = linspace(min(x8), max(x8), 100); % Generate x values for plotting
yPred5 = predict(mdl8, xRange5'); % Predict y values from the model
plot(xRange5, yPred5, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2);
text(2.2, 1, sprintf('SOZ pVal = %.4f', mdl7.ModelFitVsNullModel.Pvalue))
text(2.2, 1.1, sprintf('nonSOZ pVal = %.4f', mdl8.ModelFitVsNullModel.Pvalue))
% axis properties
ylabel('Suppression Latency (s)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Suppression Duration (s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Suppression Latency vs. Suppression Duration', ...
      'FontSize', 14, 'FontWeight', 'bold');
% Add grid lines for better readability
grid on;
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'Box', 'on', 'LineWidth', 1.5);
axis square
hold off;
% Maximize figure window
set(gcf, 'Units', 'normalized', 'Position', [0 0 1 1]);

% Scatter Plot: SOZ Suppression Amplitude vs. Suppression Duration
% Create figure and scatter plots
subplot(2,3,5);
hold on;
% Extract data for IN
temp_y5 = [sigsupp_SPESstruct(SOZsigSupp).Minimum_SuppressionAll];
temp_x5 = [sigsupp_SPESstruct(SOZsigSupp).suppressionDurationTimeAll];
y5 = temp_y5(temp_x5 ~= 0);
x5 = temp_x5(temp_x5 ~= 0);
% Extract data for PC
temp_y6 = [sigsupp_SPESstruct(~SOZsigSupp).Minimum_SuppressionAll];
temp_x6 = [sigsupp_SPESstruct(~SOZsigSupp).suppressionDurationTimeAll];
y6 = temp_y6(temp_x6 ~= 0);
x6 = temp_x6(temp_x6 ~= 0);
scatter(x6, y6, 80, 'k', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
scatter(x5, y5, 80, 'y', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5); % plot SOZ second so we can see it
% Fit linear regression models
mdl5 = fitlm([sigsupp_SPESstruct(SOZsigSupp).suppressionDurationTimeAll], [sigsupp_SPESstruct(SOZsigSupp).Minimum_SuppressionAll], 'Exclude', [sigsupp_SPESstruct(SOZsigSupp).suppressionDurationTimeAll] == 0); % Regression model for sage green data
mdl6 = fitlm([sigsupp_SPESstruct(~SOZsigSupp).suppressionDurationTimeAll], [sigsupp_SPESstruct(~SOZsigSupp).Minimum_SuppressionAll], 'Exclude', [sigsupp_SPESstruct(~SOZsigSupp).suppressionDurationTimeAll] == 0); % Regression model for muted purple data
% Plot regression lines
xRange3 = linspace(min(x5), max(x5), 100); % Generate x values for plotting
yPred3 = predict(mdl5, xRange3'); % Predict y values from the model
plot(xRange3, yPred3, 'Color', 'y', 'LineStyle', '-', 'LineWidth', 2);
xRange4 = linspace(min(x6), max(x6), 100); % Generate x values for plotting
yPred4 = predict(mdl6, xRange4'); % Predict y values from the model
plot(xRange4, yPred4, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2);
text(2, 60, sprintf('SOZ pVal = %.4f', mdl5.ModelFitVsNullModel.Pvalue))
text(2, 55, sprintf('nonSOZ pVal = %.4f', mdl6.ModelFitVsNullModel.Pvalue))
% axis properties
ylabel('Suppression Amplitude (spks/s)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Suppression Duration (s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Suppression Amplitude vs. Suppression Duration', ...
      'FontSize', 14, 'FontWeight', 'bold');
% Add grid lines for better readability
grid on;
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'Box', 'on', 'LineWidth', 1.5);
axis square
hold off;
% Maximize figure window
set(gcf, 'Units', 'normalized', 'Position', [0 0 1 1]);

% Scatter Plot: SOZ Suppression Amplitude vs. Suppression Latency
% Create figure and scatter plots
subplot(2,3,6);
hold on;
% Extract data for IN
temp_y11 = [sigsupp_SPESstruct(SOZsigSupp).Minimum_SuppressionAll];
temp_x11 = [sigsupp_SPESstruct(SOZsigSupp).suppLatencyAll];
y11 = temp_y11(temp_x11 ~= 0);
x11 = temp_x11(temp_x11 ~= 0);
% Extract data for PC
temp_y12 = [sigsupp_SPESstruct(~SOZsigSupp).Minimum_SuppressionAll];
temp_x12 = [sigsupp_SPESstruct(~isolatedClusterSigSupp).suppLatencyAll];
y12 = temp_y12(temp_x12 ~= 0);
x12 = temp_x12(temp_x12 ~= 0);
scatter(x11, y11, 80, 'y', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
scatter(x12, y12, 80, 'k', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
% Fit linear regression models
mdl11 = fitlm([sigsupp_SPESstruct(SOZsigSupp).suppLatencyAll], [sigsupp_SPESstruct(SOZsigSupp).Minimum_SuppressionAll]); % Regression model for sage green data
mdl12 = fitlm([sigsupp_SPESstruct(~SOZsigSupp).suppLatencyAll], [sigsupp_SPESstruct(~SOZsigSupp).Minimum_SuppressionAll]); % Regression model for muted purple data
% Plot regression lines
xRange2 = linspace(min(x11), max(x11), 100); % Generate x values for plotting
yPred2 = predict(mdl11, xRange2'); % Predict y values from the model
plot(xRange2, yPred2, 'Color', 'y', 'LineStyle', '-', 'LineWidth', 2);
xRange3 = linspace(min(x12), max(x12), 100); % Generate x values for plotting
yPred3 = predict(mdl12, xRange3'); % Predict y values from the model
plot(xRange3, yPred3, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2);
text(2, 60, sprintf('SOZ pVal = %.4f', mdl11.ModelFitVsNullModel.Pvalue))
text(2, 55, sprintf('nonSOZ pVal = %.4f', mdl12.ModelFitVsNullModel.Pvalue))
% axis properties
ylabel('Suppression Amplitude (spks/s)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Suppression Latency (s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Suppression Amplitude vs. Suppression Duration', ...
      'FontSize', 14, 'FontWeight', 'bold');
% Add grid lines for better readability
grid on;
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'Box', 'on', 'LineWidth', 1.5);
axis square
hold off;
% Maximize figure window
set(gcf, 'Units', 'normalized', 'Position', [0 0 1 1]);

% saving figure 
saveas(123,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\Regression\',sprintf('Regressions_CellType_SOZ_suppressionFeats_NEW.pdf')))

 %% ~~~~~~~ 1B. Cell Type: Average firing rate responses for suppressed units ~~~~~~~~~ %%

% (1) cell type: baseline frequency (-1.1 - -0.1s) (raw)
mean_BLfreq_IN_SigSupp = mean([sigsupp_SPESstruct(isolatedClusterSigSupp).BLfreqAll], "omitnan");
std_BLfreq_IN_SigSupp = std([sigsupp_SPESstruct(isolatedClusterSigSupp).BLfreqAll], "omitnan");
mean_BLfreq_PC_SigSupp = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp).BLfreqAll], "omitnan");
std_BLfreq_PC_SigSupp = std([sigsupp_SPESstruct(~isolatedClusterSigSupp).BLfreqAll], "omitnan"); 

[BLfreq_sigsupp_p,BLfreq_sigsupp_h, BLfreq_sigsupp_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp).BLfreqAll],[sigsupp_SPESstruct(~isolatedClusterSigSupp).BLfreqAll]); % BL freq

% (2) cell type: post-stim frequency  (0.1 - 1.1s) (raw)
mean_PSfreq_IN_Sig = mean([sigsupp_SPESstruct(isolatedClusterSigSupp).PSfreqAll], "omitnan");
std_PSfreq_IN_Sig = std([sigsupp_SPESstruct(isolatedClusterSigSupp).PSfreqAll], "omitnan");
mean_PSfreq_PC_Sig = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp).PSfreqAll], "omitnan");
std_PSfreq_PC_Sig = std([sigsupp_SPESstruct(~isolatedClusterSigSupp).PSfreqAll], "omitnan");

[PSfreq_sigsupp_p,PSfreq_sigsupp_h, PSfreq_sigsupp_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp).PSfreqAll],[sigsupp_SPESstruct(~isolatedClusterSigSupp).PSfreqAll]); % PS freq

% (3) cell type: FR Hz change pre to post (raw)
mean_FR_IN_Sig = mean([sigsupp_SPESstruct(isolatedClusterSigSupp).RdiffAll], "omitnan");
std_FR_IN_Sig = std([sigsupp_SPESstruct(isolatedClusterSigSupp).RdiffAll], "omitnan");
mean_FR_PC_Sig = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp).RdiffAll], "omitnan");
std_FR_PC_Sig = std([sigsupp_SPESstruct(~isolatedClusterSigSupp).RdiffAll], "omitnan");

[Rdiff_p,Rdiff_h,Rdiff_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp).RdiffAll],[sigsupp_SPESstruct(~isolatedClusterSigSupp).RdiffAll],'Tail','both'); % suppression amplitude

% (4) cell type: significant suppression amplitude (raw)
mean_suppAmp_IN = mean([sigsupp_SPESstruct(isolatedClusterSigSupp).Minimum_SuppressionAll], "omitnan");
std_suppAmp_IN = std([sigsupp_SPESstruct(isolatedClusterSigSupp).Minimum_SuppressionAll], "omitnan");
mean_suppAmp_PC = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp).Minimum_SuppressionAll], "omitnan");
std_suppAmp_PC = std([sigsupp_SPESstruct(~isolatedClusterSigSupp).Minimum_SuppressionAll], "omitnan");

[suppAmp_p,suppAmp_h,suppAmp_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp).Minimum_SuppressionAll],[sigsupp_SPESstruct(~isolatedClusterSigSupp).Minimum_SuppressionAll],'Tail','both'); % suppression amplitude

% (5) cell type: significant suppression latency
mean_suppLatency_IN = mean([sigsupp_SPESstruct(isolatedClusterSigSupp).suppLatencyAll], "omitnan");
std_suppLatency_IN = std([sigsupp_SPESstruct(isolatedClusterSigSupp).suppLatencyAll], "omitnan");
mean_suppLatency_PC = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp).suppLatencyAll], "omitnan");
std_suppLatency_PC = std([sigsupp_SPESstruct(~isolatedClusterSigSupp).suppLatencyAll], "omitnan");

[suppLatency_p,suppLatency_h,suppLatency_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp).suppLatencyAll],[sigsupp_SPESstruct(~isolatedClusterSigSupp).suppLatencyAll],'Tail','both'); % suppression latency

% Loop through the structure and replace zeros with NaN
for i = 1:length(sigsupp_SPESstruct)
    if sigsupp_SPESstruct(i).suppressionDurationTimeAll == 0
        sigsupp_SPESstruct(i).suppressionDurationTimeAll = NaN;
    end
end

% (6) cell type: significant suppression duration
mean_suppDuration_IN = mean([sigsupp_SPESstruct(isolatedClusterSigSupp).suppressionDurationTimeAll], "omitnan");
std_suppDuration_IN = std([sigsupp_SPESstruct(isolatedClusterSigSupp).suppressionDurationTimeAll], "omitnan");
mean_suppDuration_PC = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp).suppressionDurationTimeAll], "omitnan");
std_suppDuration_PC = std([sigsupp_SPESstruct(~isolatedClusterSigSupp).suppressionDurationTimeAll], "omitnan");

[suppDuration_p,suppDuration_h,suppDuration_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp).suppressionDurationTimeAll],[sigsupp_SPESstruct(~isolatedClusterSigSupp).suppressionDurationTimeAll],'Tail','both'); % suppression Duration Time

% (7) cell type: significant suppression threshold
mean_suppThreshold_IN = mean([sigsupp_SPESstruct(isolatedClusterSigSupp).durationThresholdAll], "omitnan");
std_suppThreshold_IN = std([sigsupp_SPESstruct(isolatedClusterSigSupp).durationThresholdAll], "omitnan");
mean_suppThreshold_PC = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp).durationThresholdAll], "omitnan");
std_suppThreshold_PC = std([sigsupp_SPESstruct(~isolatedClusterSigSupp).durationThresholdAll], "omitnan");

[suppThreshold_p,suppThreshold_h,suppThreshold_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp).durationThresholdAll],[sigsupp_SPESstruct(~isolatedClusterSigSupp).durationThresholdAll],'Tail','both'); % suppression Duration Time

% (8) cell type: significant recoveryrate
mean_recoveryRate_IN = mean([sigsupp_SPESstruct(isolatedClusterSigSupp).recoveryRate_Time], "omitnan");
std_recoveryRate_IN = std([sigsupp_SPESstruct(isolatedClusterSigSupp).recoveryRate_Time], "omitnan");
mean_recoveryRate_PC = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp).recoveryRate_Time], "omitnan");
std_recoveryRate_PC = std([sigsupp_SPESstruct(~isolatedClusterSigSupp).recoveryRate_Time], "omitnan");

[recoveryRate_p,recoveryRate_h,recoveryRate_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp).recoveryRate_Time],[sigsupp_SPESstruct(~isolatedClusterSigSupp).recoveryRate_Time],'Tail','both'); % suppression recovery Time

 %% ~~~~~~~ 1C. Cell Type: Z-Scored firing rate responses for suppressed units ~~~~~~~~~ %%

% (1) cell type: baseline frequency (-1.1 - -0.1s) (SIG)
mean_BLfreqZ_IN_SigSupp = mean([sigsupp_SPESstruct(isolatedClusterSigSupp).BLfreqAllz], "omitnan");
mean_BLfreqZ_PC_SigSupp = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp).BLfreqAllz], "omitnan");

[BLfreqZ_sigsupp_p,BLfreqZ_sigsupp_h, BLfreqZ_sigsupp_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp).BLfreqAllz],[sigsupp_SPESstruct(~isolatedClusterSigSupp).BLfreqAllz]); % BL freq

%  (2) cell type: post-stim frequency  (0.1 - 1.1s) (SIG)
mean_PSfreqZ_IN_Sig = mean([sigsupp_SPESstruct(isolatedClusterSigSupp).PSfreqAllz], "omitnan");
mean_PSfreqZ_PC_Sig = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp).PSfreqAllz], "omitnan");

[PSfreqZ_sigsupp_p,PSfreqZ_sigsupp_h, PSfreqZ_sigsupp_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp).PSfreqAllz],[sigsupp_SPESstruct(~isolatedClusterSigSupp).PSfreqAllz]); % PS freq

% (3) cell type FR Hz change pre to post NS
mean_FRz_IN_Sig = mean([sigsupp_SPESstruct(isolatedClusterSigSupp).RdiffAllz], "omitnan");
mean_FRz_PC_Sig = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp).RdiffAllz], "omitnan");

[Rdiffz_p,Rdiffz_h,Rdiffz_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp).RdiffAllz],[sigsupp_SPESstruct(~isolatedClusterSigSupp).RdiffAllz],'Tail','both'); % suppression amplitude

% (4) cell type: significant suppression amplitude (SIG)
mean_suppAmpZ_IN = mean([sigsupp_SPESstruct(isolatedClusterSigSupp).Minimum_SuppressionAllz], "omitnan");
mean_suppAmpZ_PC = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp).Minimum_SuppressionAllz], "omitnan");

[suppAmpZ_p,suppAmpZ_h,suppAmpZ_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp).Minimum_SuppressionAllz],[sigsupp_SPESstruct(~isolatedClusterSigSupp).Minimum_SuppressionAllz],'Tail','both'); % suppression amplitude

% (5) cell type: significant suppression latency 

% Extract the time data
timeData = [sigsupp_SPESstruct(:).suppLatencyAll];
% Calculate the mean and standard deviation
mu = mean(timeData, 'omitnan'); % Omit NaNs if present
sigma = std(timeData, 'omitnan'); % Omit NaNs if present
% Z-score the time data
zScores = (timeData - mu) / sigma;
% If you want to assign it back to the structure
for i = 1:length(sigsupp_SPESstruct)
    sigsupp_SPESstruct(i).SuppLatencyAllz = zScores(i); % Store in a new field
end
clear timeData

mean_suppLatencyZ_IN = mean([sigsupp_SPESstruct(isolatedClusterSigSupp).SuppLatencyAllz], "omitnan");
mean_suppLatencyZ_PC = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp).SuppLatencyAllz], "omitnan");

[suppLatencyZ_p,suppLatencyZ_h,suppLatencyZ_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp).SuppLatencyAllz],[sigsupp_SPESstruct(~isolatedClusterSigSupp).SuppLatencyAllz],'Tail','both'); % suppression latency

% Extract the time data
timeData = [sigsupp_SPESstruct(:).suppressionDurationTimeAll];
% Calculate the mean and standard deviation
mu = mean(timeData, 'omitnan'); % Omit NaNs if present
sigma = std(timeData, 'omitnan'); % Omit NaNs if present
% Z-score the time data
zScores = (timeData - mu) / sigma;
% If you want to assign it back to the structure
for i = 1:length(sigsupp_SPESstruct)
    sigsupp_SPESstruct(i).suppressionDurationTimeAll_z = zScores(i); % Store in a new field
end
clear timeData

% (6) cell type: significant suppression duration
suppDurationZ_IN = mean([sigsupp_SPESstruct(isolatedClusterSigSupp).suppressionDurationTimeAll_z], "omitnan");
suppDurationZ_PC = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp).suppressionDurationTimeAll_z], "omitnan");

[suppDurationZ_p,suppDurationZ_h,suppDurationZ_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp).suppressionDurationTimeAll_z],[sigsupp_SPESstruct(~isolatedClusterSigSupp).suppressionDurationTimeAll_z],'Tail','both'); % suppression Duration Time

% Extract the threshold data
thresholdData = [sigsupp_SPESstruct(:).durationThresholdAll];
% Calculate the mean and standard deviation
mu = mean(thresholdData, 'omitnan'); % Omit NaNs if present
sigma = std(thresholdData, 'omitnan'); % Omit NaNs if present
% Z-score the time data
zScores = (thresholdData - mu) / sigma;
% If you want to assign it back to the structure
for i = 1:length(sigsupp_SPESstruct)
    sigsupp_SPESstruct(i).durationThresholdAll_z = zScores(i); % Store in a new field
end

% (7) cell type: significant suppression threshold
mean_suppThresholdz_IN = mean([sigsupp_SPESstruct(isolatedClusterSigSupp).durationThresholdAll_z], "omitnan");
mean_suppThresholdz_PC = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp).durationThresholdAll_z], "omitnan");

[suppThresholdz_p,suppThresholdz_h,suppThresholdz_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp).durationThresholdAll_z],[sigsupp_SPESstruct(~isolatedClusterSigSupp).durationThresholdAll_z],'Tail','both'); % suppression Duration Time

% (8) cell type: significant recoveryrate
% % Extract the recovery data
% recoveryData = [sigsupp_SPESstruct(:).recoveryRate_Time];
% % Calculate the mean and standard deviation
% mu = mean(recoveryData, 'omitnan'); % Omit NaNs if present
% sigma = std(recoveryData, 'omitnan'); % Omit NaNs if present
% % Z-score the time data
% zScores = (recoveryData - mu) / sigma;
% % If you want to assign it back to the structure
% for i = 1:length(sigsupp_SPESstruct)
%      sigsupp_SPESstruct(i).recoveryRate_Time_z = zScores(i); % Store in a new field
% end

%mean_recoveryRatez_IN = mean([sigsupp_SPESstruct(isolatedClusterSigSupp).recoveryRate_Time_z], "omitnan");
%mean_recoveryRatez_PC = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp).recoveryRate_Time_z], "omitnan");

%[recoveryRatez_p,recoveryRatez_h,recoveryRatez_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp).recoveryRate_Time_z],[sigsupp_SPESstruct(~isolatedClusterSigSupp).recoveryRate_Time_z],'Tail','both'); % suppression recovery Time

 %% ~~~~~~~ 2A. SOZ/non-SOZ: Average firing rate responses for suppressed units ~~~~~~~~~ %%

% (1) SOZ: baseline for  SOZ units
mean_BL_SOZ = mean([sigsupp_SPESstruct(SOZsigSupp).BLfreqAll], "omitnan");
std_BL_SOZ = std([sigsupp_SPESstruct(SOZsigSupp).BLfreqAll], "omitnan");
mean_BL_nonSOZ = mean([sigsupp_SPESstruct(~SOZsigSupp).BLfreqAll], "omitnan");
std_BL_nonSOZ = std([sigsupp_SPESstruct(~SOZsigSupp).BLfreqAll], "omitnan");

[BLSOZ_p,BLSOZ_h,BLSOZ_stats] = ranksum([sigsupp_SPESstruct(SOZsigSupp).BLfreqAll],[sigsupp_SPESstruct(~SOZsigSupp).BLfreqAll],'Tail','both'); % excitation amplitude

% (2) SOZ: post-stimulation for  SOZ units
mean_PS_SOZ = mean([sigsupp_SPESstruct(SOZsigSupp).PSfreqAll], "omitnan");
std_PS_SOZ = std([sigsupp_SPESstruct(SOZsigSupp).PSfreqAll], "omitnan");
mean_PS_nonSOZ = mean([sigsupp_SPESstruct(~SOZsigSupp).PSfreqAll], "omitnan");
std_PS_nonSOZ = std([sigsupp_SPESstruct(~SOZsigSupp).PSfreqAll], "omitnan");

[PSSOZ_p,PSSOZ_h,PSSOZ_stats] = ranksum([sigsupp_SPESstruct(SOZsigSupp).PSfreqAll],[sigsupp_SPESstruct(~SOZsigSupp).PSfreqAll],'Tail','both'); % excitation amplitude

% (3) SOZ: FR difference R DIFF
mean_Rdiff_SOZ = mean([sigsupp_SPESstruct(SOZsigSupp).RdiffAll], "omitnan");
std_Rdiff_SOZ = std([sigsupp_SPESstruct(SOZsigSupp).RdiffAll], "omitnan");
mean_Rdiff_nonSOZ = mean([sigsupp_SPESstruct(~SOZsigSupp).RdiffAll], "omitnan");
std_Rdiff_nonSOZ = std([sigsupp_SPESstruct(~SOZsigSupp).RdiffAll], "omitnan");

[RdiffSOZ_p,RdiffSOZ_h,RdiffSOZ_stats] = ranksum([sigsupp_SPESstruct(SOZsigSupp).RdiffAll],[sigsupp_SPESstruct(~SOZsigSupp).RdiffAll],'Tail','both'); % excitation amplitude

% (4) SOZ: Minimum Suppression Amplitude
mean_suppAmp_SOZ = mean([sigsupp_SPESstruct(SOZsigSupp).Minimum_SuppressionAll], "omitnan");
std_suppAmp_SOZ = std([sigsupp_SPESstruct(SOZsigSupp).Minimum_SuppressionAll], "omitnan");
mean_suppAmp_nonSOZ = mean([sigsupp_SPESstruct(~SOZsigSupp).Minimum_SuppressionAll], "omitnan");
std_suppAmp_nonSOZ = std([sigsupp_SPESstruct(~SOZsigSupp).Minimum_SuppressionAll], "omitnan");

[suppAmpSOZ_p, suppAmpSOZ_h,suppAmpSOZ_stats] = ranksum([sigsupp_SPESstruct(SOZsigSupp).Minimum_SuppressionAll],[sigsupp_SPESstruct(~SOZsigSupp).Minimum_SuppressionAll],'Tail','both'); % excitation amplitude

% (5) Suppression Latency
mean_suppLatency_SOZ = mean([sigsupp_SPESstruct(SOZsigSupp).suppLatencyAll], "omitnan");
std_suppLatency_SOZ = std([sigsupp_SPESstruct(SOZsigSupp).suppLatencyAll], "omitnan");
mean_suppLatency_nonSOZ = mean([sigsupp_SPESstruct(~SOZsigSupp).suppLatencyAll], "omitnan");
std_suppLatency_nonSOZ = std([sigsupp_SPESstruct(~SOZsigSupp).suppLatencyAll], "omitnan");

[suppLatencySOZ_p, suppLatencySOZ_h,suppLatencySOZ_stats] = ranksum([sigsupp_SPESstruct(SOZsigSupp).suppLatencyAll],[sigsupp_SPESstruct(~SOZsigSupp).suppLatencyAll],'Tail','both'); % excitation amplitude

% (6) SOZ: significant suppression duration
mean_suppDuration_SOZ = mean([sigsupp_SPESstruct(SOZsigSupp).suppressionDurationTimeAll], "omitnan");
std_suppDuration_SOZ = std([sigsupp_SPESstruct(SOZsigSupp).suppressionDurationTimeAll], "omitnan");
mean_suppDuration_nonSOZ = mean([sigsupp_SPESstruct(~SOZsigSupp).suppressionDurationTimeAll], "omitnan");
std_suppDuration_nonSOZ = std([sigsupp_SPESstruct(~SOZsigSupp).suppressionDurationTimeAll], "omitnan");

[suppDurationSOZ_p, suppDurationSOZ_h,suppDurationSOZ_stats] = ranksum([sigsupp_SPESstruct(SOZsigSupp).suppressionDurationTimeAll],[sigsupp_SPESstruct(~SOZsigSupp).suppressionDurationTimeAll],'Tail','both'); % excitation amplitude

% (7) SOZ: Suppression Threshold
mean_suppThreshold_SOZ = mean([sigsupp_SPESstruct(SOZsigSupp).durationThresholdAll], "omitnan");
std_suppThreshold_SOZ = std([sigsupp_SPESstruct(SOZsigSupp).durationThresholdAll], "omitnan");
mean_suppThreshold_nonSOZ = mean([sigsupp_SPESstruct(~SOZsigSupp).durationThresholdAll], "omitnan");
std_suppThreshold_nonSOZ = std([sigsupp_SPESstruct(~SOZsigSupp).durationThresholdAll], "omitnan");

[suppThresholdSOZ_p,suppThresholdSOZ_h,suppThresholdSOZ_stats] = ranksum([sigsupp_SPESstruct(SOZsigSupp).durationThresholdAll],[sigsupp_SPESstruct(~SOZsigSupp).durationThresholdAll],'Tail','both'); % excitation amplitude

% (8) SOZ: significant recoveryrate
mean_recoveryRate_SOZ = mean([sigsupp_SPESstruct(SOZsigSupp).recoveryRate_Time], "omitnan");
std_recoveryRate_SOZ = std([sigsupp_SPESstruct(SOZsigSupp).recoveryRate_Time], "omitnan");
mean_recoveryRate_nonSOZ = mean([sigsupp_SPESstruct(~SOZsigSupp).recoveryRate_Time], "omitnan");
std_recoveryRate_nonSOZ = std([sigsupp_SPESstruct(~SOZsigSupp).recoveryRate_Time], "omitnan");

[recoveryRateSOZ_p,recoveryRateSOZ_h,recoveryRateSOZ_stats] = ranksum([sigsupp_SPESstruct(SOZsigSupp).recoveryRate_Time],[sigsupp_SPESstruct(~SOZsigSupp).recoveryRate_Time],'Tail','both'); 

% not sig.
hold on 
hist([sigsupp_SPESstruct(~SOZsigSupp).recoveryRate_Time])
hist([sigsupp_SPESstruct(SOZsigSupp).recoveryRate_Time])
hold off

 %% ~~~~~~~ 2B. SOZ/non-SOZ: Z-Scored firing rate responses for suppressed units ~~~~~~~~~ %%

% (1) baseline  SOZ units
mean_BLz_SOZ = mean([sigsupp_SPESstruct(SOZsigSupp).BLfreqAllz], "omitnan");
mean_BLz_nonSOZ = mean([sigsupp_SPESstruct(~SOZsigSupp).BLfreqAllz], "omitnan");

[BLSOZz_p,BLSOZz_h,BLSOZz_stats] = ranksum([sigsupp_SPESstruct(SOZsigSupp).BLfreqAllz],[sigsupp_SPESstruct(~SOZsigSupp).BLfreqAllz],'Tail','both'); % excitation amplitude

% (2) SOZ: post-stimulation for  SOZ units
mean_PSz_SOZ = mean([sigsupp_SPESstruct(SOZsigSupp).PSfreqAll], "omitnan");
mean_PSz_nonSOZ = mean([sigsupp_SPESstruct(~SOZsigSupp).PSfreqAll], "omitnan");

[PSSOZz_p,PSSOZz_h,PSSOZz_stats] = ranksum([sigsupp_SPESstruct(SOZsigSupp).PSfreqAllz],[sigsupp_SPESstruct(~SOZsigSupp).PSfreqAllz],'Tail','both'); % excitation amplitude

% (3) SOZ: FR difference R DIFF
mean_Rdiffz_SOZ = mean([sigsupp_SPESstruct(SOZsigSupp).RdiffAllz], "omitnan");
mean_Rdiffz_nonSOZ = mean([sigsupp_SPESstruct(~SOZsigSupp).RdiffAllz], "omitnan");

[RdiffSOZz_p,RdiffSOZz_h,RdiffSOZz_stats] = ranksum([sigsupp_SPESstruct(SOZsigSupp).RdiffAllz],[sigsupp_SPESstruct(~SOZsigSupp).RdiffAllz],'Tail','both'); % excitation amplitude

% (4) SOZ: Minimum Suppression Amplitude
mean_suppAmpz_SOZ = mean([sigsupp_SPESstruct(SOZsigSupp).Minimum_SuppressionAllz], "omitnan");
mean_suppAmpz_nonSOZ = mean([sigsupp_SPESstruct(~SOZsigSupp).Minimum_SuppressionAllz], "omitnan");

[suppAmpSOZz_p,suppAmpSOZz_h,suppAmpSOZz_stats] = ranksum([sigsupp_SPESstruct(SOZsigSupp).Minimum_SuppressionAllz],[sigsupp_SPESstruct(~SOZsigSupp).Minimum_SuppressionAllz],'Tail','both'); % excitation amplitude

% (5) SOZ: Suppression Latency 
mean_suppLatencyz_SOZ = mean([sigsupp_SPESstruct(SOZsigSupp).SuppLatencyAllz]);
mean_suppLatencyz_nonSOZ = mean([sigsupp_SPESstruct(~SOZsigSupp).SuppLatencyAllz]);

[suppLatencySOZz_p,suppLatencySOZz_h,suppLatencySOZz_stats] = ranksum([sigsupp_SPESstruct(SOZsigSupp).SuppLatencyAllz],[sigsupp_SPESstruct(~SOZsigSupp).SuppLatencyAllz],'Tail','both'); % excitation amplitude

% Extract the time data (copied from above need to edit)
timeDataSOZ = [sigsupp_SPESstruct(:).suppressionDurationTimeAll];
% Calculate the mean and standard deviation
muSOZ = mean(timeDataSOZ, 'omitnan'); % Omit NaNs if present
sigmaSOZ = std(timeDataSOZ, 'omitnan'); % Omit NaNs if present
% Z-score the time data
zScoresSOZ = (timeDataSOZ - muSOZ) / sigmaSOZ;
% If you want to assign it back to the structure
for i = 1:length(sigsupp_SPESstruct)
    sigsupp_SPESstruct(i).suppressionDurationTimeAll_z = zScoresSOZ(i); % Store in a new field
end

% (7) SOZ: Suppression Duration
mean_suppDurationz_SOZ = mean([sigsupp_SPESstruct(SOZsigSupp).suppressionDurationTimeAll_z]);
mean_suppDurationz_nonSOZ = mean([sigsupp_SPESstruct(~SOZsigSupp).suppressionDurationTimeAll_z]);

[suppDurationSOZz_p,suppDurationSOZz_h,suppDurationSOZz_stats] = ranksum([sigsupp_SPESstruct(SOZsigSupp).suppressionDurationTimeAll_z],[sigsupp_SPESstruct(~SOZsigSupp).suppressionDurationTimeAll_z],'Tail','both'); % excitation amplitude

% (7) SOZ: Suppression Threshold
mean_suppThresholdz_SOZ = mean([sigsupp_SPESstruct(SOZsigSupp).durationThresholdAll_z]);
mean_suppThresholdz_nonSOZ = mean([sigsupp_SPESstruct(~SOZsigSupp).durationThresholdAll_z]);

[suppThresholdSOZz_p,suppThresholdSOZz_h,suppThresholdSOZz_stats] = ranksum([sigsupp_SPESstruct(SOZsigSupp).durationThresholdAll_z],[sigsupp_SPESstruct(~SOZsigSupp).durationThresholdAll_z],'Tail','both'); % excitation amplitude

%(8) AUC
AUC_SOZz = zscore([sigsupp_SPESstruct(SOZsigSupp).AUC])
AUC_nonSOZz = zscore([sigsupp_SPESstruct(~SOZsigSupp).AUC])

[AUCSOZz_p,AUCSOZz_h,AUCSOZz_stats] = ranksum(AUC_SOZz,AUC_nonSOZz,'Tail','both'); % excitation amplitude

% (9) Rate of Recovery (SIG)
mean_recoveryRate_SOZz = mean([sigsupp_SPESstruct(SOZsigSupp).recoveryRate_Time_z], "omitnan");
mean_recoveryRate_nonSOZz = mean([sigsupp_SPESstruct(~SOZsigSupp).recoveryRate_Time_z], "omitnan");

[recoveryRateSOZz_p,recoveryRateSOZz_h,recoveryRateSOZz_stats] = ranksum([sigsupp_SPESstruct(SOZsigSupp).recoveryRate_Time],[sigsupp_SPESstruct(~SOZsigSupp).recoveryRate_Time_z],'Tail','both'); 


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 3A. SOZ/CellType: Average firing rate responses for suppressed units ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%

% (1) SOZ btwn Celltypes: Baseline Frequency 
mean_BLfreq_IN_SOZSigSupp = mean([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).BLfreqAll], "omitnan");
std_BLfreq_IN_SOZSigSupp = std([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).BLfreqAll], "omitnan");
mean_BLfreq_PC_SOZSigSupp = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).BLfreqAll], "omitnan");
std_BLfreq_PC_SOZSigSupp = std([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).BLfreqAll], "omitnan"); 

[BLfreq_SOZsigsupp_p,BLfreq_SOZsigsupp_h, BLfreq_SOZsigsupp_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).BLfreqAll],[sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).BLfreqAll]); % BL freq

% (2) SOZ btwn Celltypes: Post-Stim Frequency
mean_PSfreq_IN_SOZSig = mean([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).PSfreqAll], "omitnan");
std_PSfreq_IN_SOZSig = std([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).PSfreqAll], "omitnan");
mean_PSfreq_PC_SOZSig = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).PSfreqAll], "omitnan");
std_PSfreq_PC_SOZSig = std([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).PSfreqAll], "omitnan"); 

[PSfreq_SOZsigsupp_p,PSfreq_SOZsigsupp_h, PSfreq_SOZsigsupp_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).PSfreqAll],[sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).PSfreqAll]); % PS freq

% (3) SOZ btwn Celltypes: R diff
mean_FR_IN_SOZSig = mean([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).RdiffAll], "omitnan");
std_FR_IN_SOZSig = std([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).RdiffAll], "omitnan");
mean_FR_PC_SOZSig = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).RdiffAll], "omitnan");
std_FR_PC_SOZSig = std([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).RdiffAll], "omitnan");

[Rdiffz_SOZp,Rdiffz_SOZh,Rdiffz_SOZstats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).RdiffAll],[sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).RdiffAll],'Tail','both'); % suppression amplitude

% (4) SOZ btwn Celltypes: Suppression Amplitude
mean_suppAmp_SOZIN = mean([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).Minimum_SuppressionAll], "omitnan");
std_suppAmp_SOZIN = std([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).Minimum_SuppressionAll], "omitnan");
mean_suppAmp_SOZPC = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).Minimum_SuppressionAll], "omitnan");
std_suppAmp_SOZPC = std([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).Minimum_SuppressionAll], "omitnan");

[suppAmp_SOZp,suppAmp_SOZh,suppAmp_SOZstats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).Minimum_SuppressionAll],[sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).Minimum_SuppressionAll],'Tail','both'); % suppression amplitude

% (5) SOZ btwn Celltypes: Suppression Latency 
mean_suppLatency_SOZIN = mean([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).suppLatencyAll], "omitnan");
std_suppLatency_SOZIN = std([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).suppLatencyAll], "omitnan");
mean_suppLatency_SOZPC = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).suppLatencyAll], "omitnan");
std_suppLatency_SOZPC = std([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).suppLatencyAll], "omitnan");

[suppLatency_SOZ_ct_p,suppLatency_SOZh,suppLatency_SOZstats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).suppLatencyAll],[sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).suppLatencyAll],'Tail','both'); % suppression latency

% (6) SOZ btwn Celltypes: Suppression Duration
mean_suppDuration_SOZIN = mean([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).suppressionDurationTimeAll], "omitnan");
std_suppDuration_SOZIN = std([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).suppressionDurationTimeAll], "omitnan");
mean_suppDuration_SOZPC = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).suppressionDurationTimeAll], "omitnan");
std_suppDuration_SOZPC = std([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).suppressionDurationTimeAll], "omitnan");

[suppDuration_SOZp,suppDuration_SOZh,suppDuration_SOZstats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).suppressionDurationTimeAll],[sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).suppressionDurationTimeAll],'Tail','both'); % suppression Duration Time

% (7) SOZ btwn Celltypes: Suppression Threshold
mean_suppDuration_SOZIN = mean([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).durationThresholdAll], "omitnan");
std_suppDuration_SOZIN = std([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).durationThresholdAll], "omitnan");
mean_suppDuration_SOZPC = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).durationThresholdAll], "omitnan");
std_suppDuration_SOZPC = std([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).durationThresholdAll], "omitnan");

[suppThreshold_SOZp,suppThreshold_SOZh,suppThreshold_SOZstats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).durationThresholdAll],[sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).durationThresholdAll],'Tail','both'); % suppression Duration Time

% (8) SOZ: significant recoveryrate
mean_recoveryRate_SOZIN = mean([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).recoveryRate_Time], "omitnan");
std_recoveryRate_SOZIN = std([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).recoveryRate_Time], "omitnan");
mean_recoveryRate_nonSOZPC = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).recoveryRate_Time], "omitnan");
std_recoveryRate_nonSOZPC = std([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).recoveryRate_Time], "omitnan");

% PCs between SOZ/nonSOZ
[recoveryRateSOZ_p,recoveryRateSOZ_h,recoveryRateSOZ_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).recoveryRate_Time],[sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).recoveryRate_Time],'Tail','both'); 
% PCs between SOZ/nonSOZ
[recoveryRateSOZ_p,recoveryRateSOZ_h,recoveryRateSOZ_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp(~SOZsigSupp)).recoveryRate_Time],[sigsupp_SPESstruct(~isolatedClusterSigSupp(~SOZsigSupp)).recoveryRate_Time],'Tail','both'); 

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 3B. SOZ/CellType: Zscored firing rate responses for suppressed units ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%

% (1) SOZ btwn Celltypes: Baseline Frequency 
mean_BLfreqZ_IN_SOZSigSupp = mean([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).BLfreqAllz], "omitnan");
mean_BLfreqZ_PC_SOZSigSupp = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).BLfreqAllz], "omitnan");

[BLfreqZ_SOZ_ct_p,BLfreqZ_SOZ_ct_h, BLfreqZ_SOZ_ct_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).BLfreqAllz],[sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).BLfreqAllz]); % BL freq

% (2) SOZ btwn Celltypes: Post-Stim Frequency
mean_PSfreqZ_IN_SOZSig = mean([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).PSfreqAllz], "omitnan");
mean_PSfreqZ_PC_SOZSig = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).PSfreqAllz], "omitnan");

[PSfreqZ_SOZ_ct_p,PSfreqZ_SOZ_ct_h, PSfreqZ_SOZ_ct_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).PSfreqAllz],[sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).PSfreqAllz]); % PS freq

% (3) SOZ btwn Celltypes: R diff
mean_FRz_IN_SOZSig = mean([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).RdiffAllz], "omitnan");
mean_FRz_PC_SOZSig = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).RdiffAllz], "omitnan");

[Rdiffz_SOZ_ct_p,Rdiffz_SOZ_ct_h,Rdiffz_SOZ_ct_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).RdiffAllz],[sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).RdiffAllz],'Tail','both'); % suppression amplitude

% (4) SOZ btwn Celltypes: Suppression Amplitude
mean_suppAmpZ_SOZIN = mean([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).Minimum_SuppressionAllz], "omitnan");
mean_suppAmpZ_SOZPC = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).Minimum_SuppressionAllz], "omitnan");

[suppAmp_SOZ_ct_p,suppAmp_SOZ_ct_h,suppAmp_SOZ_ct_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).Minimum_SuppressionAllz],[sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).Minimum_SuppressionAllz],'Tail','both'); % suppression amplitude

% (5) SOZ btwn Celltypes: Suppression Latency (if weird check the for loop)
mean_suppLatencyZ_SOZIN = mean([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).SuppLatencyAllz], "omitnan");
mean_suppLatencyZ_SOZPC = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).SuppLatencyAllz], "omitnan");

[suppLatency_SOZ_ct_p,suppLatency_SOZ_ct_h,suppLatency_SOZ_ct_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).SuppLatencyAllz],[sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).SuppLatencyAllz],'Tail','both'); % suppression latency

% (6) SOZ btwn Celltypes: Suppression Duration (if weird check the for loop)
mean_suppDurationZ_SOZIN = mean([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).suppressionDurationTimeAll_z], "omitnan");
mean_suppDurationZ_SOZPC = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).suppressionDurationTimeAll_z], "omitnan");

[suppDuration_SOZ_ct_p,suppDuration_SOZ_ct_h,suppDuration_SOZ_ct_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).suppressionDurationTimeAll_z],[sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).suppressionDurationTimeAll_z],'Tail','both'); % suppression Duration Time

% (7) SOZ btwn Celltypes: Suppression Threshold (if weird check the for loop)
mean_suppThresholdz_SOZIN = mean([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).durationThresholdAll_z], "omitnan");
mean_suppThresholdz_SOZPC = mean([sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).durationThresholdAll_z], "omitnan");

[suppThresholdz_SOZ_ct_p,suppThresholdz_SOZ_ct_h,suppThresholdz_SOZ_ct_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).durationThresholdAll_z],[sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).durationThresholdAll_z],'Tail','both'); % suppression Duration Time


% (8) recovery rate

[recoveryRateSOZz_p,recoveryRateSOZz_h,recoveryRateSOZz_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).recoveryRate_Time_z],[sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).recoveryRate_Time_z],'Tail','both'); 
[recoveryRateSOZz_p,recoveryRateSOZz_h,recoveryRateSOZz_stats] = ranksum([sigsupp_SPESstruct(isolatedClusterSigSupp(~SOZsigSupp)).recoveryRate_Time_z],[sigsupp_SPESstruct(isolatedClusterSigSupp(SOZsigSupp)).recoveryRate_Time_z],'Tail','both'); 
[recoveryRateSOZz_p,recoveryRateSOZz_h,recoveryRateSOZz_stats] = ranksum([sigsupp_SPESstruct(~isolatedClusterSigSupp(~SOZsigSupp)).recoveryRate_Time_z],[sigsupp_SPESstruct(~isolatedClusterSigSupp(SOZsigSupp)).recoveryRate_Time_z],'Tail','both'); 


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 4A. All: Averaged firing rate responses for enhanced units ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%

% (1) Enhanced: Baseline Frequency
mean_BL_FR = mean([sigsexci_SPESstruct(:).BLfreqAll], "omitnan");
std_BL_FR = std([sigsexci_SPESstruct(:).BLfreqAll], "omitnan");

% (2) Enhanced: Post-Stimulation Frequency
mean_PS_FR = mean([sigsexci_SPESstruct(:).PSfreqAll], "omitnan");
std_PS_FR = std([sigsexci_SPESstruct(:).PSfreqAll], "omitnan");

% (3) Enhanced: Maximum FR post-stim
mean_max_FR = mean([sigsexci_SPESstruct(:).Maximum_FRAll], "omitnan");
std_max_FR = std([sigsexci_SPESstruct(:).Maximum_FRAll], "omitnan");

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 4B. Enhanced btwn Cell Types: Averaged firing rate responses for enhanced units ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
keyboard
% (1) Baseline Frequency
mean_BL_FR_IN = mean([sigsexci_SPESstruct(isolatedClusterSigExci).BLfreqAll], "omitnan");
std_BL_FR_IN = std([sigsexci_SPESstruct(isolatedClusterSigExci).BLfreqAll], "omitnan");
mean_BL_FR_PC = mean([sigsexci_SPESstruct(~isolatedClusterSigExci).BLfreqAll], "omitnan");
std_BL_FR_PC = std([sigsexci_SPESstruct(~isolatedClusterSigExci).BLfreqAll], "omitnan");

[BLFR_p,BLSFR_h,BLFR_stats] = ranksum([sigsexci_SPESstruct(isolatedClusterSigExci).BLfreqAll],[sigsexci_SPESstruct(~isolatedClusterSigExci).BLfreqAll],'Tail','both'); 

% (2) Enhanced: Post-Stimulation Frequency
mean_PS_FR_IN = mean([sigsexci_SPESstruct(isolatedClusterSigExci).PSfreqAll], "omitnan");
std_PS_FR_IN = std([sigsexci_SPESstruct(isolatedClusterSigExci).PSfreqAll], "omitnan");
mean_PS_FR_PC = mean([sigsexci_SPESstruct(~isolatedClusterSigExci).PSfreqAll], "omitnan");
std_PS_FR_PC = std([sigsexci_SPESstruct(~isolatedClusterSigExci).PSfreqAll], "omitnan");

[PSFR_p,PSFR_h,PSFR_stats] = ranksum([sigSPESstruct(isolatedClusterSigExci).PSfreqAll],[sigSPESstruct(~isolatedClusterSigExci).PSfreqAll],'Tail','both');

% (3) Enhanced: Maximum FR post-stim
mean_max_FR_IN = mean([sigsexci_SPESstruct(isolatedClusterSigExci).Maximum_FRAll], "omitnan");
std_max_FR_IN = std([sigsexci_SPESstruct(isolatedClusterSigExci).Maximum_FRAll], "omitnan");
mean_max_FR_PC = mean([sigsexci_SPESstruct(~isolatedClusterSigExci).Maximum_FRAll], "omitnan");
std_max_FR_PC = std([sigsexci_SPESstruct(~isolatedClusterSigExci).Maximum_FRAll], "omitnan");

[maxFR_p,maxFR_h,maxFR_stats] = ranksum([sigSPESstruct(isolatedClusterSigExci).Maximum_FRAll],[sigSPESstruct(~isolatedClusterSigExci).Maximum_FRAll],'Tail','both');

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 4C. Enhanced btwn Cell Types: Averaged firing rate responses for enhanced units ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%

% (1) Baseline Frequency
mean_BL_FRz_IN = mean([sigsexci_SPESstruct(isolatedClusterSigExci).BLfreqAllz], "omitnan");
mean_BL_FRz_PC = mean([sigsexci_SPESstruct(~isolatedClusterSigExci).BLfreqAllz], "omitnan");

[BLFRz_p,BLSFRz_h,BLFRz_stats] = ranksum([sigsexci_SPESstruct(isolatedClusterSigExci).BLfreqAllz],[sigsexci_SPESstruct(~isolatedClusterSigExci).BLfreqAllz],'Tail','both'); 

% (2) Enhanced: Post-Stimulation Frequency
mean_PS_FRz_IN = mean([sigsexci_SPESstruct(isolatedClusterSigExci).PSfreqAllz], "omitnan");
mean_PS_FRz_PC = mean([sigsexci_SPESstruct(~isolatedClusterSigExci).PSfreqAllz], "omitnan");

[PSFRz_p,PSFRz_h,PSFRz_stats] = ranksum([sigSPESstruct(isolatedClusterSigExci).PSfreqAllz],[sigSPESstruct(~isolatedClusterSigExci).PSfreqAllz],'Tail','both');

% (3) Enhanced: Maximum FR post-stim
mean_max_FRz_IN = mean([sigsexci_SPESstruct(isolatedClusterSigExci).Maximum_FRAllz], "omitnan");
mean_max_FRz_PC = mean([sigsexci_SPESstruct(~isolatedClusterSigExci).Maximum_FRAllz], "omitnan");

[maxFRz_p,maxFRz_h,maxFRz_stats] = ranksum([sigSPESstruct(isolatedClusterSigExci).Maximum_FRAllz],[sigSPESstruct(~isolatedClusterSigExci).Maximum_FRAllz],'Tail','both');


%scatter([sigsupp_SPESstruct(:).durationThresholdAll],[sigsupp_SPESstruct(:).suppressionDurationTimeAll])


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Analysis AREA UNDER THE CURVE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%

% Cell type
AUC = [SPESstruct(:).AUC];
IN_AUC = [sigsupp_SPESstruct(isolatedClusterSigSupp).AUC];
PC_AUC = [sigsupp_SPESstruct(~isolatedClusterSigSupp).AUC];

sigAUC = [sigSPESstruct(:).AUC];

sum(IN_AUC > 0) % count of INs that have an AUC: 18
sum(PC_AUC > 0)
IN_AUC(IN_AUC == 0) = [];
PC_AUC(PC_AUC == 0) = [];

AUC(AUC == 0) = [];

sigAUCz = zscore(sigAUC);

IN_AUC_z = zscore(IN_AUC);
PC_AUC_z = zscore(PC_AUC);

meanIN_AUC = mean(IN_AUC);
meanPC_AUC = mean(PC_AUC);
stdIN_AUC = std(IN_AUC);
stdPC_AUC = std(PC_AUC);

[AUC_p, AUC_h,AUC_stats] = ranksum(IN_AUC,PC_AUC,'Tail','both'); % AUC
[AUCz_p, AUCz_h,AUCz_stats] = ranksum(IN_AUC_z,PC_AUC_z,'Tail','both'); % AUCz

suppressionAUC = [suppressionSigUnits', AUC']

[p,h,stats] = ranksum(AUC(isolatedCluster), AUC(~isolatedCluster));

figure;
% Plot the histograms
hold on;
numBins = 10;
histogram(PC_AUC, numBins, 'FaceColor', muted_purple, 'FaceAlpha', 0.7, 'EdgeColor', 'k');
histogram(IN_AUC, numBins, 'FaceColor', sage_green, 'FaceAlpha', 0.7, 'EdgeColor', 'k');
hold off;
xlabel('AUC');
ylabel('Frequency');
legend('PC AUC', 'IN AUC');
title('Histogram of AUC Values');

figure;
% Plot the histograms
hold on;
numBins = 10;
histogram(PC_AUC_z, numBins, 'FaceColor', muted_purple, 'FaceAlpha', 0.7, 'EdgeColor', 'k');
histogram(IN_AUC_z, numBins, 'FaceColor', sage_green, 'FaceAlpha', 0.7, 'EdgeColor', 'k');
hold off;
xlabel('AUCz');
ylabel('Frequency');
legend('PC AUCz', 'IN AUCz');
title('Histogram of AUC z-Values');

% SOZ
SOZ_AUC = [sigsupp_SPESstruct(SOZsigSupp).AUC];
nonSOZ_AUC = [sigsupp_SPESstruct(~SOZsigSupp).AUC];

sum(SOZ_AUC > 0) % count of INs that have an AUC: 18
sum(nonSOZ_AUC > 0)
SOZ_AUC(SOZ_AUC == 0) = [];
nonSOZ_AUC(nonSOZ_AUC == 0) = [];

SOZ_AUC_z = zscore(SOZ_AUC);
nonSOZ_AUC_z = zscore(nonSOZ_AUC);

meanSOZ_AUC = mean(SOZ_AUC)
meannonSOZ_AUC = mean(nonSOZ_AUC)

[SOZ_AUC_p, SOZ_AUC_h,SOZ_AUC_stats] = ranksum(SOZ_AUC,nonSOZ_AUC,'Tail','both'); % AUC SIG but doesnt matter 
[SOZ_AUCz_p, SOZ_AUCz_h,SOZ_AUCz_stats] = ranksum(SOZ_AUC_z,nonSOZ_AUC_z,'Tail','both'); % AUCz NS

figure;
% Plot the histograms
hold on;
numBins = 10;
histogram(nonSOZ_AUC, numBins, 'FaceColor', 'k', 'FaceAlpha', 0.6, 'EdgeColor', 'k');
histogram(SOZ_AUC, numBins, 'FaceColor', 'y', 'FaceAlpha', 0.9, 'EdgeColor', 'k');
hold off;
xlabel('AUC');
ylabel('Frequency');
legend('SOZ AUC', 'nonSOZ AUC');
title('Histogram of AUC Values');

figure;
% Plot the histograms
hold on;
numBins = 10;
histogram(nonSOZ_AUC_z, numBins, 'FaceColor', 'k', 'FaceAlpha', 0.6, 'EdgeColor', 'k');
histogram(SOZ_AUC_z, numBins, 'FaceColor', 'y', 'FaceAlpha', 0.9, 'EdgeColor', 'k');
hold off;
xlabel('AUCz');
ylabel('Frequency');
legend('SOZ AUCz', 'nonSOZ AUCz');
title('Histogram of AUC z-Values');

%% Sub-Analysis 
% Which regions evoke the greatest stimulation responses? Amperage change?

% subset Locs as a logical value for just significant  values.
sigAMYLocs = AMYLocs(sigAll) % 25
sigOFCLocs = OFCLocs(sigAll) % 72
sigHIPPLocs = HIPPLocs(sigAll) % 32
sigCINGLocs = CINGLocs(sigAll) % 45
sigpCINGLocs = pCINGLocs(sigAll) % 1

sigsuppAMYLocs = AMYLocs(suppressionSigUnits) % 25
sigsuppOFCLocs = OFCLocs(suppressionSigUnits) % 72
sigsuppHIPPLocs = HIPPLocs(suppressionSigUnits) % 32
sigsuppCINGLocs = CINGLocs(suppressionSigUnits) % 45
sigsupppCINGLocs = pCINGLocs(suppressionSigUnits) % 1

% PLOT: Barchart for suppression Sig Units
sigSuppRegions = [sum(AMYLocs) sum(sigsuppAMYLocs); sum(OFCLocs) sum(sigsuppOFCLocs); sum(HIPPLocs) sum(sigsuppHIPPLocs); sum(CINGLocs) sum(sigsuppCINGLocs); sum(pCINGLocs) sum(sigsupppCINGLocs)];
bar(sigSuppRegions)

% (1) Calculate baseline frequencies for units in each region
meanBLfreq_AMY = mean([sigSPESstruct(sigAMYLocs).BLfreqAll], "omitnan");
meanBLfreq_OFC = mean([sigSPESstruct(sigOFCLocs).BLfreqAll], "omitnan");
meanBLfreq_HIPP = mean([sigSPESstruct(sigHIPPLocs).BLfreqAll], "omitnan");
meanBLfreq_CING = mean([sigSPESstruct(sigCINGLocs).BLfreqAll], "omitnan");
meanBLfreq_pCING = mean([sigSPESstruct(sigpCINGLocs).BLfreqAll], "omitnan");

% (2) Calculate post stimulation FR for units in each region (ALL)
meanPSfreq_AMY = mean([sigSPESstruct(sigAMYLocs).PSfreqAll], "omitnan");
meanPSfreq_OFC = mean([sigSPESstruct(sigOFCLocs).PSfreqAll], "omitnan");
meanPSfreq_HIPP = mean([sigSPESstruct(sigHIPPLocs).PSfreqAll], "omitnan");
meanPSfreq_CING = mean([sigSPESstruct(sigCINGLocs).PSfreqAll], "omitnan");
meanPSfreq_pCING = mean([sigSPESstruct(sigpCINGLocs).PSfreqAll], "omitnan");

% (3) Calculate post stimulation FR for units in each region (SUPPRESSION)
meanPSsupp_AMY = mean([sigSPESstruct(sigAMYLocs & suppressionSigUnits).PSfreqAll], "omitnan");
meanPSsupp_OFC = mean([sigSPESstruct(sigOFCLocs & suppressionSigUnits).PSfreqAll], "omitnan");
meanPSsupp_HIPP = mean([sigSPESstruct(sigHIPPLocs & suppressionSigUnits).PSfreqAll], "omitnan");
meanPSsupp_CING = mean([sigSPESstruct(sigCINGLocs & suppressionSigUnits).PSfreqAll], "omitnan");
meanPSsupp_pCING = mean([sigSPESstruct(sigpCINGLocs & suppressionSigUnits).PSfreqAll], "omitnan");

% (4) Calculate difference between baseline frequency and post-stimulation suppression change for each region
suppressionAmperage_sigAMY = [sigSPESstruct(sigAMYLocs & suppressionSigUnits).BLfreqAll] - [sigSPESstruct(sigAMYLocs & suppressionSigUnits).Minimum_SuppressionAll];
suppressionAmperage_sigOFC = [sigSPESstruct(sigOFCLocs & suppressionSigUnits).BLfreqAll] - [sigSPESstruct(sigOFCLocs & suppressionSigUnits).Minimum_SuppressionAll];
suppressionAmperage_sigHIPP = [sigSPESstruct(sigHIPPLocs & suppressionSigUnits).BLfreqAll] - [sigSPESstruct(sigHIPPLocs & suppressionSigUnits).Minimum_SuppressionAll];
suppressionAmperage_sigCING = [sigSPESstruct(sigCINGLocs & suppressionSigUnits).BLfreqAll] - [sigSPESstruct(sigCINGLocs & suppressionSigUnits).Minimum_SuppressionAll];
suppressionAmperage_sigpCING = [sigSPESstruct(sigpCINGLocs & suppressionSigUnits).BLfreqAll] - [sigSPESstruct(sigpCINGLocs & suppressionSigUnits).Minimum_SuppressionAll];

mean(suppressionAmperage_sigAMY)
mean(suppressionAmperage_sigOFC)
mean(suppressionAmperage_sigHIPP)
mean(suppressionAmperage_sigCING)
mean(suppressionAmperage_sigpCING)

% colors for regions
Muted_Blue =  [0.4, 0.5, 0.6]; % amygdala
Muted_Green = [0.45, 0.55, 0.45]; % cingulate
Dark_Blue = [0.1, 0.15, 0.3]; % hippocampus
Muted_Purple = [0.55, 0.406, 0.6]; % OFC
Dark_Brown = [0.6588,0.5490,0.4902]; % posterior cingulate

% Evoked firing rate changes between different regions
figure(22)
hold on
title("Regional: Baseline Frequency -> Minimum Suppression")
% Amygdala
scatter(ones(length(suppressionAmperage_sigAMY),1),[sigSPESstruct(sigAMYLocs & suppressionSigUnits).BLfreqAll], 'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', Muted_Blue)
scatter(2*ones(length(suppressionAmperage_sigAMY),1),[sigSPESstruct(sigAMYLocs & suppressionSigUnits).Minimum_SuppressionAll],'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', Muted_Blue)
line([ones(length(suppressionAmperage_sigAMY),1) 2*ones(length(suppressionAmperage_sigAMY),1)]', [[sigSPESstruct(sigAMYLocs & suppressionSigUnits).BLfreqAll]; [sigSPESstruct(sigAMYLocs & suppressionSigUnits).Minimum_SuppressionAll]], 'Color', Muted_Blue)
% OFC
scatter(3*ones(length(suppressionAmperage_sigOFC),1),[sigSPESstruct(sigOFCLocs & suppressionSigUnits).BLfreqAll], 'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', Muted_Purple)
scatter(4*ones(length(suppressionAmperage_sigOFC),1),[sigSPESstruct(sigOFCLocs & suppressionSigUnits).Minimum_SuppressionAll],'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', Muted_Purple)
line([3*ones(length(suppressionAmperage_sigOFC),1) 4*ones(length(suppressionAmperage_sigOFC),1)]', [[sigSPESstruct(sigOFCLocs & suppressionSigUnits).BLfreqAll]; [sigSPESstruct(sigOFCLocs & suppressionSigUnits).Minimum_SuppressionAll]], 'Color', Muted_Purple)
% Hippocampus
scatter(5*ones(length(suppressionAmperage_sigHIPP),1),[sigSPESstruct(sigHIPPLocs & suppressionSigUnits).BLfreqAll], 'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', Dark_Blue)
scatter(6*ones(length(suppressionAmperage_sigHIPP),1),[sigSPESstruct(sigHIPPLocs & suppressionSigUnits).Minimum_SuppressionAll],'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', Dark_Blue)
line([5*ones(length(suppressionAmperage_sigHIPP),1) 6*ones(length(suppressionAmperage_sigHIPP),1)]', [[sigSPESstruct(sigHIPPLocs & suppressionSigUnits).BLfreqAll]; [sigSPESstruct(sigHIPPLocs & suppressionSigUnits).Minimum_SuppressionAll]], 'Color', Dark_Blue)
% Cingulate
scatter(7*ones(length(suppressionAmperage_sigCING),1),[sigSPESstruct(sigCINGLocs & suppressionSigUnits).BLfreqAll], 'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', Muted_Green)
scatter(8*ones(length(suppressionAmperage_sigCING),1),[sigSPESstruct(sigCINGLocs & suppressionSigUnits).Minimum_SuppressionAll],'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', Muted_Green)
line([7*ones(length(suppressionAmperage_sigCING),1) 8*ones(length(suppressionAmperage_sigCING),1)]', [[sigSPESstruct(sigCINGLocs & suppressionSigUnits).BLfreqAll]; [sigSPESstruct(sigCINGLocs & suppressionSigUnits).Minimum_SuppressionAll]], 'Color', Muted_Green)
% Posterior Cingulate
scatter(9*ones(length(suppressionAmperage_sigpCING),1),[sigSPESstruct(sigpCINGLocs & suppressionSigUnits).BLfreqAll], 'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', Dark_Brown)
scatter(10*ones(length(suppressionAmperage_sigpCING),1),[sigSPESstruct(sigpCINGLocs & suppressionSigUnits).Minimum_SuppressionAll],'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', Dark_Brown)
line([9*ones(length(suppressionAmperage_sigpCING),1) 10*ones(length(suppressionAmperage_sigpCING),1)]', [[sigSPESstruct(sigpCINGLocs & suppressionSigUnits).BLfreqAll]; [sigSPESstruct(sigpCINGLocs & suppressionSigUnits).Minimum_SuppressionAll]], 'Color', Dark_Brown)
xlim([0.5 10.5])
xticks([1.5 3.5 5.5 7.5 9.5])
xticklabels({'Amygdala','OFC & vmPFC','Hippocampus', 'Cingulate', 'posterior Cingulate'})
axis square
hold off
% saving figure 
saveas(22,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\WaveFormClassification\',sprintf('Regional_UnitResponsesFRChange_NEW.pdf')))
close(22)

%% ~~~~ Regional Neuronal firing Metrics~~~~~ %%

% BL
BL_AMY = ([sigSPESstruct(sigAMYLocs & suppressionSigUnits).BLfreqAll])
BL_OFC = ([sigSPESstruct(sigOFCLocs & suppressionSigUnits).BLfreqAll])
BL_HIPP = ([sigSPESstruct(sigHIPPLocs & suppressionSigUnits).BLfreqAll])
BL_CING = ([sigSPESstruct(sigCINGLocs & suppressionSigUnits).BLfreqAll])
BL_pCING = ([sigSPESstruct(sigpCINGLocs & suppressionSigUnits).BLfreqAll])

% ONE-WAY ANOVA FOR BASELINE (removed pCING from this analysis)
y_BL_Region = [BL_AMY, BL_OFC, BL_HIPP, BL_CING];
group_BL_Region = repelem(1:4, 1, [numel(BL_AMY'),numel(BL_OFC'),numel(BL_HIPP'), numel(BL_CING')]);

[p_BL_Region, tbl_BL_Region, stats_BL_Region] = anova1(y_BL_Region,group_BL_Region)
if p_BL_Region <= 0.05
    results_BL_Region = multcompare(stats_AUC_Region)
else
end % No significant differences.

% min suppression
minSuppression_AMY = ([sigSPESstruct(sigAMYLocs & suppressionSigUnits).Minimum_SuppressionAll])
minSuppression_OFC = ([sigSPESstruct(sigOFCLocs & suppressionSigUnits).Minimum_SuppressionAll])
minSuppression_HIPP = ([sigSPESstruct(sigHIPPLocs & suppressionSigUnits).Minimum_SuppressionAll])
minSuppression_CING = ([sigSPESstruct(sigCINGLocs & suppressionSigUnits).Minimum_SuppressionAll])
mminSuppression_pCING = ([sigSPESstruct(sigpCINGLocs & suppressionSigUnits).Minimum_SuppressionAll])

% ONE-WAY ANOVA FOR MINIMUM SUPPRESSION (removed pCING from this analysis)
y_minSupp_Region = [minSuppression_AMY, minSuppression_OFC, minSuppression_HIPP, minSuppression_CING];
group_minSupp_Region = repelem(1:4, 1, [numel(minSuppression_AMY'),numel(minSuppression_OFC'),numel(minSuppression_HIPP'), numel(minSuppression_CING')]);

[p_minSupp_Region, tbl_minSupp_Region, stats_minSupp_Region] = anova1(y_minSupp_Region,group_minSupp_Region)
if p_minSupp_Region <= 0.05
    results_minSupp_Region = multcompare(stats_minSupp_Region)
else
end % not sig

% suppression latency
SuppressionLat_AMY = ([sigSPESstruct(sigAMYLocs & suppressionSigUnits).suppLatencyAll])
SuppressionLat_OFC = ([sigSPESstruct(sigOFCLocs & suppressionSigUnits).suppLatencyAll])
SuppressionLat_HIPP = ([sigSPESstruct(sigHIPPLocs & suppressionSigUnits).suppLatencyAll])
SuppressionLat_CING = ([sigSPESstruct(sigCINGLocs & suppressionSigUnits).suppLatencyAll])
SuppressionLat_pCING = ([sigSPESstruct(sigpCINGLocs & suppressionSigUnits).suppLatencyAll])

% ONE-WAY ANOVA FOR SUPPRESSION LATENCY (removed pCING from this analysis)
y_suppLat_Region = [SuppressionLat_AMY, SuppressionLat_OFC, SuppressionLat_HIPP, SuppressionLat_CING];
group_suppLat_Region = repelem(1:4, 1, [numel(SuppressionLat_AMY'),numel(SuppressionLat_OFC'),numel(SuppressionLat_HIPP'), numel(SuppressionLat_CING')]);

[p_suppLat_Region, tbl_suppLat_Region, stats_suppLat_Region] = anova1(y_suppLat_Region,group_suppLat_Region)
if p_suppLat_Region <= 0.05
    results_suppLat_Region = multcompare(stats_suppLat_Region)
else
end % OFC and CING are significantly different.

% suppression duration
duration_AMY = ([sigSPESstruct(sigAMYLocs & suppressionSigUnits).suppressionDurationTimeAll])
duration_OFC = ([sigSPESstruct(sigOFCLocs & suppressionSigUnits).suppressionDurationTimeAll])
duration_HIPP = ([sigSPESstruct(sigHIPPLocs & suppressionSigUnits).suppressionDurationTimeAll])
duration_CING = ([sigSPESstruct(sigCINGLocs & suppressionSigUnits).suppressionDurationTimeAll])
duration_pCING = ([sigSPESstruct(sigpCINGLocs & suppressionSigUnits).suppressionDurationTimeAll])

% ONE-WAY ANOVA FOR DURATION (removed pCING from this analysis)
y_duration_Region = [duration_AMY, duration_OFC, duration_HIPP, duration_CING];
group_duration_Region = repelem(1:4, 1, [numel(duration_AMY'),numel(duration_OFC'),numel(duration_HIPP'), numel(duration_CING')]);

[p_duration_Region, tbl_duration_Region, stats_duration_Region] = anova1(y_duration_Region,group_duration_Region)
if p_duration_Region <= 0.05
    results_duration_Region = multcompare(stats_duration_Region)
else
end % No significant differences.

% suppression duration threshold
durationthreshold_AMY = ([sigSPESstruct(sigAMYLocs & suppressionSigUnits).durationThresholdAll])
durationthreshold_OFC = ([sigSPESstruct(sigOFCLocs & suppressionSigUnits).durationThresholdAll])
durationthreshold_HIPP = ([sigSPESstruct(sigHIPPLocs & suppressionSigUnits).durationThresholdAll])
durationthreshold_CING = ([sigSPESstruct(sigCINGLocs & suppressionSigUnits).durationThresholdAll])
durationthreshold_pCING = ([sigSPESstruct(sigpCINGLocs & suppressionSigUnits).durationThresholdAll])

% ONE-WAY ANOVA FOR DURATION THRESHOLD (removed pCING from this analysis)
y_durationthreshold_Region = [durationthreshold_AMY, durationthreshold_OFC, durationthreshold_HIPP, durationthreshold_CING];
group_durationthreshold_Region = repelem(1:4, 1, [numel(durationthreshold_AMY'),numel(durationthreshold_OFC'),numel(durationthreshold_HIPP'), numel(durationthreshold_CING')]);

[p_durationthreshold_Region, tbl_durationthreshold_Region, stats_durationthreshold_Region] = anova1(y_durationthreshold_Region,group_durationthreshold_Region)
if p_durationthreshold_Region <= 0.05
    results_durationthreshold_Region = multcompare(stats_durationthreshold_Region)
else
end % No significant differences.

% suppression recoveryRate_Time
recoveryRate_AMY = ([sigSPESstruct(sigAMYLocs & suppressionSigUnits).recoveryRate_Time])
recoveryRate_OFC = ([sigSPESstruct(sigOFCLocs & suppressionSigUnits).recoveryRate_Time])
recoveryRate_HIPP = ([sigSPESstruct(sigHIPPLocs & suppressionSigUnits).recoveryRate_Time])
recoveryRate_CING = ([sigSPESstruct(sigCINGLocs & suppressionSigUnits).recoveryRate_Time])
recoveryRate_pCING = ([sigSPESstruct(sigpCINGLocs & suppressionSigUnits).recoveryRate_Time])

% ONE-WAY ANOVA FOR recoveryRate_Time (removed pCING from this analysis)
y_recoveryRate_Region = [recoveryRate_AMY, recoveryRate_OFC, recoveryRate_HIPP, recoveryRate_CING];
group_recoveryRate_Region = repelem(1:4, 1, [numel(recoveryRate_AMY'),numel(recoveryRate_OFC'),numel(recoveryRate_HIPP'), numel(recoveryRate_CING')]);

[p_recoveryRate_Region, tbl_recoveryRate_Region, stats_recoveryRate_Region] = anova1(y_recoveryRate_Region,group_recoveryRate_Region)
if p_recoveryRate_Region <= 0.05
    results_recoveryRate_Region = multcompare(stats_recoveryRate_Region)
else
end % No significant differences.

% AUC
AUC_AMY = ([sigSPESstruct(sigAMYLocs & suppressionSigUnits).AUC])
AUC_OFC = ([sigSPESstruct(sigOFCLocs & suppressionSigUnits).AUC])
AUC_HIPP = ([sigSPESstruct(sigHIPPLocs & suppressionSigUnits).AUC])
AUC_CING = ([sigSPESstruct(sigCINGLocs & suppressionSigUnits).AUC])
AUC_pCING = ([sigSPESstruct(sigpCINGLocs & suppressionSigUnits).AUC])

% ONE-WAY ANOVA FOR AUC (removed pCING from this analysis)
y_AUC_Region = [AUC_AMY, AUC_OFC, AUC_HIPP, AUC_CING];
group_AUC_Region = repelem(1:4, 1, [numel(AUC_AMY'),numel(AUC_OFC'),numel(AUC_HIPP'), numel(AUC_CING')]);

[p_AUC_Region, tbl_AUC_Region, stats_AUC_Region] = anova1(y_AUC_Region,group_AUC_Region)
if p_AUC_Region <= 0.05
    results_AUC_Region = multcompare(stats_AUC_Region)
else
end % No significant differences.

% suppression Amerage
y_AMP_Region = [suppressionAmperage_sigAMY, suppressionAmperage_sigOFC, suppressionAmperage_sigHIPP, suppressionAmperage_sigCING];
group_AMP_Region = repelem(1:4, 1, [numel(suppressionAmperage_sigAMY'),numel(suppressionAmperage_sigOFC'),numel(suppressionAmperage_sigHIPP'), numel(suppressionAmperage_sigCING')]);

[p_AMP_Region, tbl_AMP_Region, stats_AMP_Region] = anova1(y_AMP_Region,group_AMP_Region)
if p_AMP_Region <= 0.05
    results_AMP_Region = multcompare(stats_AMP_Region)
else
end % No significant differences.

% post firing rate
PSsupp_AMY = [sigSPESstruct(sigAMYLocs & suppressionSigUnits).PSfreqAll];
PSsupp_OFC = [sigSPESstruct(sigOFCLocs & suppressionSigUnits).PSfreqAll];
PSsupp_HIPP = [sigSPESstruct(sigHIPPLocs & suppressionSigUnits).PSfreqAll];
PSsupp_CING = [sigSPESstruct(sigCINGLocs & suppressionSigUnits).PSfreqAll];
PSsupp_pCING = [sigSPESstruct(sigpCINGLocs & suppressionSigUnits).PSfreqAll];

% ONE-WAY ANOVA FOR AUC (removed pCING from this analysis)
y_PS_Region = [PSsupp_AMY, PSsupp_OFC, PSsupp_HIPP, PSsupp_CING];
group_PS_Region = repelem(1:4, 1, [numel(PSsupp_AMY'),numel(PSsupp_OFC'),numel(PSsupp_HIPP'), numel(PSsupp_CING')]);

[p_PS_Region, tbl_PS_Region, stats_PS_Region] = anova1(y_PS_Region,group_PS_Region)
if p_PS_Region <= 0.05
    results_PS_Region = multcompare(stats_PS_Region)
else
end % No significant differences.

%% ~~~ Regional Neuronal Firing Metrics Z-scored) ~~~
% BLz
BL_AMY = ([sigSPESstruct(sigAMYLocs & suppressionSigUnits).BLfreqAllz])
BL_OFC = ([sigSPESstruct(sigOFCLocs & suppressionSigUnits).BLfreqAllz])
BL_HIPP = ([sigSPESstruct(sigHIPPLocs & suppressionSigUnits).BLfreqAllz])
BL_CING = ([sigSPESstruct(sigCINGLocs & suppressionSigUnits).BLfreqAllz])
BL_pCING = ([sigSPESstruct(sigpCINGLocs & suppressionSigUnits).BLfreqAllz])

% ONE-WAY ANOVA FOR BASELINE (removed pCING from this analysis)
y_BL_Region = [BL_AMY, BL_OFC, BL_HIPP, BL_CING];
group_BL_Region = repelem(1:4, 1, [numel(BL_AMY'),numel(BL_OFC'),numel(BL_HIPP'), numel(BL_CING')]);

[p_BL_Region, tbl_BL_Region, stats_BL_Region] = anova1(y_BL_Region,group_BL_Region)
if p_BL_Region <= 0.05
    results_BL_Region = multcompare(stats_AUC_Region)
else
end % Main effect, NS multcomp

% min suppressionz
minSuppression_AMY = ([sigSPESstruct(sigAMYLocs & suppressionSigUnits).Minimum_SuppressionAllz])
minSuppression_OFC = ([sigSPESstruct(sigOFCLocs & suppressionSigUnits).Minimum_SuppressionAllz])
minSuppression_HIPP = ([sigSPESstruct(sigHIPPLocs & suppressionSigUnits).Minimum_SuppressionAllz])
minSuppression_CING = ([sigSPESstruct(sigCINGLocs & suppressionSigUnits).Minimum_SuppressionAllz])
mminSuppression_pCING = ([sigSPESstruct(sigpCINGLocs & suppressionSigUnits).Minimum_SuppressionAllz])

% ONE-WAY ANOVA FOR MINIMUM SUPPRESSION (removed pCING from this analysis)
y_minSupp_Region = [minSuppression_AMY, minSuppression_OFC, minSuppression_HIPP, minSuppression_CING];
group_minSupp_Region = repelem(1:4, 1, [numel(minSuppression_AMY'),numel(minSuppression_OFC'),numel(minSuppression_HIPP'), numel(minSuppression_CING')]);

[p_minSupp_Region, tbl_minSupp_Region, stats_minSupp_Region] = anova1(y_minSupp_Region,group_minSupp_Region)
if p_minSupp_Region <= 0.05
    results_minSupp_Region = multcompare(stats_minSupp_Region)
else
end % Main effect, NS multcomp

% Extract the time data
timeData = [sigSPESstruct(:).suppLatencyAll];
% Calculate the mean and standard deviation
mu = mean(timeData, 'omitnan'); % Omit NaNs if present
sigma = std(timeData, 'omitnan'); % Omit NaNs if present
% Z-score the time data
zScores = (timeData - mu) / sigma;
% If you want to assign it back to the structure
for i = 1:length(sigSPESstruct)
    sigSPESstruct(i).SuppLatencyAllz = zScores(i); % Store in a new field
end
clear timeData

% suppression latencyz
SuppressionLat_AMY = ([sigSPESstruct(sigAMYLocs & suppressionSigUnits).SuppLatencyAllz])
SuppressionLat_OFC = ([sigSPESstruct(sigOFCLocs & suppressionSigUnits).SuppLatencyAllz])
SuppressionLat_HIPP = ([sigSPESstruct(sigHIPPLocs & suppressionSigUnits).SuppLatencyAllz])
SuppressionLat_CING = ([sigSPESstruct(sigCINGLocs & suppressionSigUnits).SuppLatencyAllz])
SuppressionLat_pCING = ([sigSPESstruct(sigpCINGLocs & suppressionSigUnits).SuppLatencyAllz])

% ONE-WAY ANOVA FOR SUPPRESSION LATENCY (removed pCING from this analysis)
y_suppLat_Region = [SuppressionLat_AMY, SuppressionLat_OFC, SuppressionLat_HIPP, SuppressionLat_CING];
group_suppLat_Region = repelem(1:4, 1, [numel(SuppressionLat_AMY'),numel(SuppressionLat_OFC'),numel(SuppressionLat_HIPP'), numel(SuppressionLat_CING')]);

[p_suppLat_Region, tbl_suppLat_Region, stats_suppLat_Region] = anova1(y_suppLat_Region,group_suppLat_Region)
if p_suppLat_Region <= 0.05
    results_suppLat_Region = multcompare(stats_suppLat_Region)
else
end % NS

for i = 1:length(sigSPESstruct)
    if isempty(sigSPESstruct(i).suppressionDurationTimeAll)
        sigSPESstruct(i).suppressionDurationTimeAll = NaN;
    end
end

% Loop through the structure and replace zeros with NaN
for i = 1:length(sigSPESstruct)
    if sigSPESstruct(i).suppressionDurationTimeAll == 0
        sigSPESstruct(i).suppressionDurationTimeAll = NaN;
    end
end

% Extract the time data
timeDurData = ([sigSPESstruct(:).suppressionDurationTimeAll]);
% Calculate the mean and standard deviation
mu = mean(timeDurData, 'omitnan'); % Omit NaNs if present
sigma = std(timeDurData, 'omitnan'); % Omit NaNs if present
% Z-score the time data
zScores = (timeDurData - mu) / sigma;
% If you want to assign it back to the structure
for i = 1:125%length(sigSPESstruct)
    sigSPESstruct(i).suppressionDurationTimeAllz = zScores(i); % Store in a new field
end
clear timeDurData

% suppression duration
duration_AMY = ([sigSPESstruct(sigAMYLocs & suppressionSigUnits).suppressionDurationTimeAllz])
duration_OFC = ([sigSPESstruct(sigOFCLocs & suppressionSigUnits).suppressionDurationTimeAllz])
duration_HIPP = ([sigSPESstruct(sigHIPPLocs & suppressionSigUnits).suppressionDurationTimeAllz])
duration_CING = ([sigSPESstruct(sigCINGLocs & suppressionSigUnits).suppressionDurationTimeAllz])
duration_pCING = ([sigSPESstruct(sigpCINGLocs & suppressionSigUnits).suppressionDurationTimeAllz])

% ONE-WAY ANOVA FOR DURATION (removed pCING from this analysis)
y_duration_Region = [duration_AMY, duration_OFC, duration_HIPP, duration_CING];
group_duration_Region = repelem(1:4, 1, [numel(duration_AMY'),numel(duration_OFC'),numel(duration_HIPP'), numel(duration_CING')]);

[p_duration_Region, tbl_duration_Region, stats_duration_Region] = anova1(y_duration_Region,group_duration_Region)
if p_duration_Region <= 0.05
    results_duration_Region = multcompare(stats_duration_Region)
else
end % No significant differences.

% Extract the time data
timeThresData = [sigSPESstruct(:).durationThresholdAll];
% Calculate the mean and standard deviation
mu = mean(timeThresData, 'omitnan'); % Omit NaNs if present
sigma = std(timeThresData, 'omitnan'); % Omit NaNs if present
% Z-score the time data
zScores = (timeThresData - mu) / sigma;
% If you want to assign it back to the structure
for i = 1:length(sigSPESstruct)
    sigSPESstruct(i).durationThresholdAllz = zScores(i); % Store in a new field
end
clear timeThresData

% suppression duration threshold
durationthreshold_AMY = ([sigSPESstruct(sigAMYLocs & suppressionSigUnits).durationThresholdAllz])
durationthreshold_OFC = ([sigSPESstruct(sigOFCLocs & suppressionSigUnits).durationThresholdAllz])
durationthreshold_HIPP = ([sigSPESstruct(sigHIPPLocs & suppressionSigUnits).durationThresholdAllz])
durationthreshold_CING = ([sigSPESstruct(sigCINGLocs & suppressionSigUnits).durationThresholdAllz])
durationthreshold_pCING = ([sigSPESstruct(sigpCINGLocs & suppressionSigUnits).durationThresholdAllz])

% ONE-WAY ANOVA FOR DURATION THRESHOLD (removed pCING from this analysis)
y_durationthreshold_Region = [durationthreshold_AMY, durationthreshold_OFC, durationthreshold_HIPP, durationthreshold_CING];
group_durationthreshold_Region = repelem(1:4, 1, [numel(durationthreshold_AMY'),numel(durationthreshold_OFC'),numel(durationthreshold_HIPP'), numel(durationthreshold_CING')]);

[p_durationthreshold_Region, tbl_durationthreshold_Region, stats_durationthreshold_Region] = anova1(y_durationthreshold_Region,group_durationthreshold_Region)
if p_durationthreshold_Region <= 0.05
    results_durationthreshold_Region = multcompare(stats_durationthreshold_Region)
else
end % No significant differences.

for i = 1:length(sigSPESstruct)
    if isempty(sigSPESstruct(i).recoveryRate_Time)
        sigSPESstruct(i).recoveryRate_Time = NaN;
    end
end

% Loop through the structure and replace zeros with NaN
for i = 1:length(sigSPESstruct)
    if sigSPESstruct(i).recoveryRate_Time == 0
        sigSPESstruct(i).recoveryRate_Time = NaN;
    end
end

% Extract the time data
timeRecData = [sigSPESstruct(:).recoveryRate_Time];
% Calculate the mean and standard deviation
mu = mean(timeRecData, 'omitnan'); % Omit NaNs if present
sigma = std(timeRecData, 'omitnan'); % Omit NaNs if present
% Z-score the time data
zScores = (timeRecData - mu) / sigma;
% If you want to assign it back to the structure
for i = 1:length(sigSPESstruct)
    sigSPESstruct(i).recoveryRate_Timez = zScores(i); % Store in a new field
end
clear timeThresData

% suppression recoveryRate_Time
recoveryRate_AMY = ([sigSPESstruct(sigAMYLocs & suppressionSigUnits).recoveryRate_Timez])
recoveryRate_OFC = ([sigSPESstruct(sigOFCLocs & suppressionSigUnits).recoveryRate_Timez])
recoveryRate_HIPP = ([sigSPESstruct(sigHIPPLocs & suppressionSigUnits).recoveryRate_Timez])
recoveryRate_CING = ([sigSPESstruct(sigCINGLocs & suppressionSigUnits).recoveryRate_Timez])
recoveryRate_pCING = ([sigSPESstruct(sigpCINGLocs & suppressionSigUnits).recoveryRate_Timez])

% ONE-WAY ANOVA FOR recoveryRate_Time (removed pCING from this analysis)
y_recoveryRate_Region = [recoveryRate_AMY, recoveryRate_OFC, recoveryRate_HIPP, recoveryRate_CING];
group_recoveryRate_Region = repelem(1:4, 1, [numel(recoveryRate_AMY'),numel(recoveryRate_OFC'),numel(recoveryRate_HIPP'), numel(recoveryRate_CING')]);

[p_recoveryRate_Region, tbl_recoveryRate_Region, stats_recoveryRate_Region] = anova1(y_recoveryRate_Region,group_recoveryRate_Region)
if p_recoveryRate_Region <= 0.05
    results_recoveryRate_Region = multcompare(stats_recoveryRate_Region)
else
end % No significant differences.

% AUC
AUC_AMY = sigAUCz(sigAMYLocs & suppressionSigUnits)
AUC_OFC = sigAUCz(sigOFCLocs & suppressionSigUnits)
AUC_HIPP = sigAUCz(sigHIPPLocs & suppressionSigUnits)
AUC_CING = sigAUCz(sigCINGLocs & suppressionSigUnits)
AUC_pCING = sigAUCz(sigpCINGLocs & suppressionSigUnits)

% ONE-WAY ANOVA FOR AUC (removed pCING from this analysis)
y_AUC_Region = [AUC_AMY, AUC_OFC, AUC_HIPP, AUC_CING];
group_AUC_Region = repelem(1:4, 1, [numel(AUC_AMY'),numel(AUC_OFC'),numel(AUC_HIPP'), numel(AUC_CING')]);

[p_AUC_Region, tbl_AUC_Region, stats_AUC_Region] = anova1(y_AUC_Region,group_AUC_Region)
if p_AUC_Region <= 0.05
    results_AUC_Region = multcompare(stats_AUC_Region)
else
end % No significant differences.

% post firing rate
PSsupp_AMY = [sigSPESstruct(sigAMYLocs & suppressionSigUnits).PSfreqAllz];
PSsupp_OFC = [sigSPESstruct(sigOFCLocs & suppressionSigUnits).PSfreqAllz];
PSsupp_HIPP = [sigSPESstruct(sigHIPPLocs & suppressionSigUnits).PSfreqAllz];
PSsupp_CING = [sigSPESstruct(sigCINGLocs & suppressionSigUnits).PSfreqAllz];
PSsupp_pCING = [sigSPESstruct(sigpCINGLocs & suppressionSigUnits).PSfreqAllz];

% ONE-WAY ANOVA FOR AUC (removed pCING from this analysis)
y_PS_Region = [PSsupp_AMY, PSsupp_OFC, PSsupp_HIPP, PSsupp_CING];
group_PS_Region = repelem(1:4, 1, [numel(PSsupp_AMY'),numel(PSsupp_OFC'),numel(PSsupp_HIPP'), numel(PSsupp_CING')]);

[p_PS_Region, tbl_PS_Region, stats_PS_Region] = anova1(y_PS_Region,group_PS_Region)
if p_PS_Region <= 0.05
    results_PS_Region = multcompare(stats_PS_Region)
else
end % No significant differences.


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EUCLIDEAN DISTANCE ACROSS SUBJECTS ANALYSIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%

% load csv with imagesc logical and direction of FR modulation
imagescTable = readtable("\\155.100.91.44\D\Data\Rhiannon\CCEPS\SPES_UNITS\imagesc.csv")

evokedFRLogical = imagescTable.Var1 == 1; % 40
nearFRLogical = imagescTable.Var1 == 1 & strcmp(imagescTable.Var2, 'N'); % 31
farFRLogical = imagescTable.Var1 == 1 & strcmp(imagescTable.Var2, 'F'); % 9
sozFRLogical = imagescTable.Var1 == 1 & strcmp(imagescTable.Var3, 'SOZ'); % 4
nonsozFRlogical = imagescTable.Var1 == 1 & strcmp(imagescTable.Var3, 'nonSOZ'); % 36
sozNearFRLogical = imagescTable.Var1 == 1 & strcmp(imagescTable.Var3, 'SOZ') & strcmp(imagescTable.Var2, 'N'); % 4
nonsozNearFRLogical = imagescTable.Var1 == 1 & strcmp(imagescTable.Var3, 'nonSOZ') & strcmp(imagescTable.Var2, 'N'); % 4

compareLogical = [sigUnits; evokedFRLogical']; % matrix of overlapping logicals.
sum(sum(compareLogical) ==2); % how many overlapping units were there in total: 30.

compareLogicalNear = [sigUnits; evokedFRLogical';nearFRLogical']; % matrix of overlapping logicals.
sum(sum(compareLogicalNear) ==3); % how many overlapping units were there for Nears: 24.

compareLogicalFar = [sigUnits; evokedFRLogical';farFRLogical']; % matrix of overlapping logicals.
sum(sum(compareLogicalFar) ==3); % how many overlapping units were there for Fars: 6.

% Region Evoked Responses
Regions = ["Amygdala"; "OFC & vmPFC"; "Hippocampus"; "Cingulate"; "Posterior Cingulate"];
x_Regions = categorical(Regions);
Muted_Blue =  [0.4, 0.5, 0.6];
Muted_Green = [0.45, 0.55, 0.45];
Dark_Blue = [0.1, 0.15, 0.3];
Muted_Purple = [0.55, 0.406, 0.6];
Dark_Brown = [0.6588,0.5490,0.4902];
AMY = "Amygdala"; % finding locations with Amygdala
AMYLocsAll = contains(unitLocs(1,:), AMY);
sum(AMYLocsAll)
OFC = ("OFC"| "Orbital" | "vmPFC"); % finding locations with OFC
OFCLocsAll = contains(unitLocs(1,:), OFC);
sum(OFCLocsAll)
HIPP = "Hippocampus"; % finding locations with HIPP
HIPPLocsAll = contains(unitLocs(1,:), HIPP);
sum(HIPPLocsAll)
CING = "Cingulate"; % finding locations with Cingulate
CINGLocsAll = contains(unitLocs(1,:), CING);
sum(CINGLocsAll)
pCING = "Posterior"; % finding locations with Cingulate
pCINGLocsAll = contains(unitLocs(1,:), pCING);
sum(pCINGLocsAll)

y_unitCount = [sum(AMYLocsAll(nearFRLogical)); sum(OFCLocsAll(nearFRLogical)); sum(HIPPLocsAll(nearFRLogical)); sum(CINGLocsAll(nearFRLogical)); sum(pCINGLocsAll(nearFRLogical))];
figure(222)
b = bar(x_Regions, y_unitCount,'FaceColor','flat');
% Assign colors to each bar
b.CData(1,:) = Muted_Blue;
b.CData(2,:) = Muted_Green;
b.CData(3,:) = Dark_Blue;
b.CData(4,:) = Muted_Purple;
b.CData(5,:) = Dark_Brown;
title('Evoked FR Unit Counts', 'FontSize', 14) % Increase title font size
xlabel('Regions', 'FontSize', 12) % Increase x-axis label font size
ylabel('Unit Count', 'FontSize', 12) % Increase y-axis label font size
ylim([0 14])
axis square
maximize(222)
% saving figure
saveas(222,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\',sprintf('EvokedFR_GroupedLocations_AcrossSubs.pdf')))
close(222)

% Cell types that had evoked firing generally 
InterneuronsEvokedStim = evokedFRLogical(isolatedCluster); % IN that have evoked FR (n = 8)
PrincipalcellEvokedStim = evokedFRLogical(~isolatedCluster); % IN that have evoked FR (n = 32)

% Cell types that had evoked firing near
InterneuronsEvokedNearStim =  nearFRLogical(isolatedCluster); % IN that have evoked FR (n = 6)
PrincipalcellEvokedNearStim = nearFRLogical(~isolatedCluster); % IN that have evoked FR (n = 22)

% Cell types that had evoked firing far
InterneuronsEvokedFarStim = farFRLogical(isolatedCluster); % IN that have evoked FR (n = 1)
PrincipalcellEvokedFarStim =  farFRLogical(~isolatedCluster); % IN that have evoked FR (n = 8)

% Proportion Test between evoked stim and cell type
% proportion test
X_IN = sum(InterneuronsEvokedNearStim); % 7
X_PC = sum(PrincipalcellEvokedNearStim); % 24
% prop test out of total units
N = [sum(isolatedCluster) sum(~isolatedCluster)]; % 37 191
[h_evokedStim_total, p_evokedStim_total, chi2stat_INevokedStim_total, df_evokedStim_total] = prop_test([X_IN X_PC], N,false);
% prop test out of evoked units.
N_evoked = [sum(InterneuronsEvokedStim) sum(PrincipalcellEvokedStim)]; % 8 32
[h_evokedStim, p_evokedStim, chi2stat_INevokedStim, df_evokedStim] = prop_test([X_IN X_PC], N_evoked,false);

% Where do these near stims come from?: MATTER
matterINstims = {SPESstruct(InterneuronsEvokedStim).electrodeMatter}; % white = 4 gray = 4
matterPCstims = {SPESstruct(PrincipalcellEvokedStim).electrodeMatter}; % white = 15 gray = 17
[uc, ~, idc] = unique(matterPCstims)
counts = accumarray(idc,ones(size(idc))); % white = 15 gray = 17
matterINnearstims = {SPESstruct(InterneuronsEvokedNearStim).electrodeMatter}; % 3/4 
matterPCnearstims = {SPESstruct(PrincipalcellEvokedNearStim).electrodeMatter}; % 10/14
[ucnear, ~, idcnear] = unique(matterPCnearstims)
counts = accumarray(idcnear,ones(size(idcnear))); % grey 14, white 10.
matterINfarstims = {SPESstruct(InterneuronsEvokedFarStim).electrodeMatter}; %white 1
matterPCfarstims = {SPESstruct(PrincipalcellEvokedFarStim).electrodeMatter}; % white 5, gray 3.

% proportion test
X_INmatter = sum(matterINnearstims);
X_PCmatter = sum(matterPCnearstims);
N = [sum(isolatedCluster) sum(~isolatedCluster)];
[h_evokedStimMatter, p_evokedStimMatter, chi2stat_INevokedStimMatter, df_evokedStimMatter] = prop_test([X_INmatter X_PCmatter], N,false);

% Where do these near stims come from?: REGION

regionsEvokedstims = {SPESstruct(nearFRLogical).chanLabel};

nearEvokedOFC = sum(contains(regionsEvokedstims(1,:), OFC));
nearEvokedAMY = sum(contains(regionsEvokedstims(1,:), AMY));
nearEvokedCING = sum(contains(regionsEvokedstims(1,:), CING));
nearEvokedHIPP = sum(contains(regionsEvokedstims(1,:), HIPP));
nearEvokedpCING = sum(contains(regionsEvokedstims(1,:), pCING));

% proportion of near evoked firings by region
propOFC_evoked = (nearEvokedOFC/sum(OFCLocsAll(sigAll)))*100
propHIPP_evoked = (nearEvokedHIPP/sum(HIPPLocsAll(sigAll)))*100
propAMY_evoked = (nearEvokedAMY/sum(AMYLocsAll(sigAll)))*100
propCING_evoked = (nearEvokedCING/sum(CINGLocsAll(sigAll)))*100
proppCING_evoked = (nearEvokedpCING/sum(pCINGLocsAll(sigAll)))*100

% pairwise proportion test [OFC,HIPP; OFC,AMY; OFC,CING; OFC,pCING; HIPP,AMY; HIPP,CING; HIPP,pCING; AMY,CING; AMY, pCING; CING,pCING]
% sum of OFC significant units
OFC_HIPP = [sum(OFCLocsAll(sigAll)) sum(HIPPLocsAll(sigAll))];
OFC_AMY = [sum(OFCLocsAll(sigAll)) sum(AMYLocsAll(sigAll))];
OFC_CING = [sum(OFCLocsAll(sigAll)) sum(CINGLocsAll(sigAll))];
OFC_pCING = [sum(OFCLocsAll(sigAll)) sum(pCINGLocsAll(sigAll))];
HIPP_AMY = [sum(HIPPLocsAll(sigAll)) sum(AMYLocsAll(sigAll))];
HIPP_CING = [sum(HIPPLocsAll(sigAll)) sum(CINGLocsAll(sigAll))];
HIPP_pCING = [sum(HIPPLocsAll(sigAll)) sum(pCINGLocsAll(sigAll))];
AMY_CING = [sum(AMYLocsAll(sigAll)) sum(CINGLocsAll(sigAll))];
AMY_pCING = [sum(AMYLocsAll(sigAll)) sum(pCINGLocsAll(sigAll))];
CING_pCING = [sum(CINGLocsAll(sigAll)) sum(pCINGLocsAll(sigAll))];

% pairwise proportion testing
[h_OFC_HIPP, p_OFC_HIPP, chi2stat_OFC_HIPP, df_OFC_HIPP] = prop_test([nearEvokedOFC nearEvokedHIPP], OFC_HIPP,false);
[h_OFC_AMY, p_OFC_AMY, chi2stat_OFC_AMY, df_OFC_AMY] = prop_test([nearEvokedOFC nearEvokedAMY], OFC_AMY,false);
[h_OFC_CING, p_OFC_CING, chi2stat_OFC_CING, df_OFC_CING] = prop_test([nearEvokedOFC nearEvokedCING], OFC_CING,false);
[h_OFC_pCING, p_OFC_pCING, chi2stat_OFC_pCING, df_OFC_pCING] = prop_test([nearEvokedOFC nearEvokedpCING], OFC_pCING,false);
[h_HIPP_AMY, p_HIPP_AMY, chi2stat_HIPP_AMY, df_HIPP_AMY] = prop_test([nearEvokedHIPP nearEvokedAMY], HIPP_AMY,false);
[h_HIPP_CING, p_HIPP_CING, chi2stat_HIPP_CING, df_HIPP_CING] = prop_test([nearEvokedHIPP nearEvokedCING], HIPP_CING,false);
[h_HIPP_pCING, p_HIPP_pCING, chi2stat_HIPP_pCING, df_HIPP_pCING] = prop_test([nearEvokedHIPP nearEvokedpCING], HIPP_pCING,false);
[h_AMY_CING, p_AMY_CING, chi2stat_AMY_CING, df_AMY_CING] = prop_test([nearEvokedAMY nearEvokedCING], AMY_CING,false);
[h_AMY_pCING, p_AMY_pCING, chi2stat_AMY_pCING, df_AMY_pCING] = prop_test([nearEvokedAMY nearEvokedpCING], AMY_pCING,false);
[h_CING_pCING, p_CING_pCING, chi2stat_CING_pCING, df_CING_pCING] = prop_test([nearEvokedCING nearEvokedpCING], CING_pCING,false);

p_values = [p_OFC_HIPP;p_OFC_AMY;p_OFC_CING;p_OFC_pCING;p_HIPP_AMY;p_HIPP_CING;p_HIPP_pCING;p_AMY_CING;p_AMY_pCING;p_CING_pCING];

% fdr doesnt work, too few comparisons.
[fdr, q] = mafdr(p_values, 'showplot', true, 'BHFDR', false);


regionINstims = {SPESstruct(InterneuronsEvokedStim).chanLabel}; % right OFC 1, left OFC 2, left amygdala 4, 1 left Anterior Hippocampus
regionPCstims = {SPESstruct(PrincipalcellEvokedStim).chanLabel}; 
[uc, ~, idc] = unique(regionPCstims)
counts = accumarray(idc,ones(size(idc))); % 5,2,4,1,1,1,4,4,1,2,3,1,1
regionINnearstims = {SPESstruct(InterneuronsEvokedNearStim).chanLabel}; % 3/3 
regionPCnearstims = {SPESstruct(PrincipalcellEvokedNearStim).chanLabel};
[ucnear, ~, idcnear] = unique(regionPCnearstims)
counts = accumarray(idcnear,ones(size(idcnear))); % grey 12, white 10.
matterINfarstims = {SPESstruct(InterneuronsEvokedFarStim).chanLabel}; %white 1
matterPCfarstims = {SPESstruct(PrincipalcellEvokedFarStim).chanLabel}; % white 5, gray 3

boxplot(idc,counts)

% Histogram of residual error from changept analysis %
residual_changepts = [SPESstruct(:).residual_changepts];
plot(histcounts(residual_changepts(nearFRLogical),20)) % just near evoked stims residual error
plot(residual_changepts(nearFRLogical), 'r')

% Histogram of Rsquared from fit %
rsquare_fit = arrayfun(@(s) s.goodness.rsquare, SPESstruct);
adjrsquare_fit = arrayfun(@(s) s.goodness.adjrsquare, SPESstruct);
histogram(rsquare_fit)
plot(histcounts(rsquare_fit(nearFRLogical),20)) % just near evoked stims residual error

% What is the firing rate within and  outside of the distance??
duration_FRBelow = [SPESstruct(:).duration_FRBelow]
meanFR_belowEvoked = [SPESstruct(nearFRLogical).meanFR_belowEvoked] % averaging the mean firing rates for near stims (or can plot individually?)
meanFR_aboveEvoked = [SPESstruct(nearFRLogical).meanFR_aboveEvoked]  % averaging the mean firing rates for far stims (or can plot individually?)

figure(1234);
% Define the number of points in each group
numBelow = length(meanFR_belowEvoked);
numAbove = length(meanFR_aboveEvoked);
% Scatter plot for meanFR_belowEvoked
hold on;
scatter(ones(numBelow, 1), meanFR_belowEvoked, 70, 'filled', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'yellow');
scatter(2 * ones(numAbove, 1), meanFR_aboveEvoked, 70, 'filled', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'blue');
for i = 1:numBelow
    line([1, 2], [meanFR_belowEvoked(i), meanFR_aboveEvoked(i)], 'Color', 'k'); % Connecting line
end
xlim([0.5 2.5]);
xticks([1, 2]);
xticklabels({'Below Evoked', 'Above Evoked'});
title('Mean FR Below and Above Evoked');
xlabel('Evoked Firing Distance');
ylabel('Mean FR (Hz)');
axis square; 
hold off;
saveas(1234,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\',sprintf('EvokedFR_BelowAboveScatter.pdf')))
close(1234)

% (1) Evoked response vs distance for SOZ  vs nonSOZ

%(A) Amplitude differences
SOZ_FR_belowEvoked = [SPESstruct(sozNearFRLogical).meanFR_belowEvoked]; 
mean(SOZ_FR_belowEvoked) % 1.26 Hz
nonSOZ_FR_belowEvoked = [SPESstruct(nonsozNearFRLogical).meanFR_belowEvoked];
mean(nonSOZ_FR_belowEvoked, 'omitNaN') % 1.90 Hz

%(B) Number of Units
evokedFRLogical
nearFRLogical
farFRLogical
sozFRLogical
nonsozFRlogical
sozNearFRLogical
nonsozNearFRLogical

%(C) Time of evoked firing?
SOZ_FR_below_Duration = [SPESstruct(sozNearFRLogical).duration_FRBelow]; 
mean(SOZ_FR_below_Duration)
nonSOZ_FR_below_Duration = [SPESstruct(nonsozNearFRLogical).duration_FRBelow]; 
mean(nonSOZ_FR_below_Duration)

SOZ_FR_above_Duration = [SPESstruct(sozNearFRLogical).duration_FRAbove]

% Suppression of transient excitation similar to those without?
% ACROSS UNITS
% (1) Baseline Freq
BLFreq_evokedFR = mean([SPESstruct(nearFRLogical).BLfreqAll], "omitnan");
BLFreq_no_evokedFR = mean([SPESstruct(~nearFRLogical).BLfreqAll], "omitnan");
BLFreq_evokedFRz = mean([SPESstruct(nearFRLogical).BLfreqAllz], "omitnan"); % -3.529
BLFreq_no_evokedFRz = mean([SPESstruct(~nearFRLogical).BLfreqAllz], "omitnan"); % -3.591

[p,h,stats] = ranksum([SPESstruct(nearFRLogical).BLfreqAllz],[SPESstruct(~nearFRLogical).BLfreqAllz],'Tail','both'); % not sig

% (2) Post Stimulation Freq
PSFreq_evokedFR = mean([SPESstruct(nearFRLogical).PSfreqAll], "omitnan");
PSFreq_no_evokedFR = mean([SPESstruct(~nearFRLogical).PSfreqAll], "omitnan");
PSFreq_evokedFRz = mean([SPESstruct(nearFRLogical).PSfreqAllz], "omitnan");
PSFreq_no_evokedFRz = mean([SPESstruct(~nearFRLogical).PSfreqAllz], "omitnan");

[p,h,stats] = ranksum([SPESstruct(nearFRLogical).PSfreqAllz],[SPESstruct(~nearFRLogical).PSfreqAllz],'Tail','both');% not sig

% (3) Suppression Amplitude
suppAmp_evokedFR = mean([SPESstruct(nearFRLogical).Minimum_SuppressionAll], "omitnan");
suppAmp_no_evokedFR = mean([SPESstruct(~nearFRLogical).Minimum_SuppressionAll], "omitnan");
suppAmp_evokedFRz = mean([SPESstruct(nearFRLogical).Minimum_SuppressionAllz], "omitnan");
suppAmp_no_evokedFRz = mean([SPESstruct(~nearFRLogical).Minimum_SuppressionAllz], "omitnan");

[p,h,stats] = ranksum([SPESstruct(nearFRLogical).Minimum_SuppressionAllz],[SPESstruct(~nearFRLogical).Minimum_SuppressionAllz],'Tail','both'); % not sig

% (4) Suppression Latency
suppLat_evokedFR = mean([SPESstruct(nearFRLogical).suppLatencyAll], "omitnan");
suppLat_no_evokedFR = mean([SPESstruct(~nearFRLogical).suppLatencyAll], "omitnan");
suppLat_evokedFRz = mean([SPESstruct(nearFRLogical).SuppLatencyAllz], "omitnan");
suppLat_no_evokedFRz = mean([SPESstruct(~nearFRLogical).SuppLatencyAllz], "omitnan");

[p,h,stats] = ranksum([SPESstruct(nearFRLogical).SuppLatencyAllz],[SPESstruct(~nearFRLogical).SuppLatencyAllz],'Tail','both'); % not sig

% (5) Suppression Duration
suppDur_evokedFR = mean([SPESstruct(nearFRLogical).suppressionDurationTimeAll], "omitnan");
suppDur_no_evokedFR = mean([SPESstruct(~nearFRLogical).suppressionDurationTimeAll], "omitnan");
suppDur_evokedFRz = mean([SPESstruct(nearFRLogical).suppressionDurationTimeAllz], "omitnan");
suppDur_no_evokedFRz = mean([SPESstruct(~nearFRLogical).suppressionDurationTimeAllz], "omitnan");

[p,h,stats] = ranksum([SPESstruct(nearFRLogical).suppressionDurationTimeAll],[SPESstruct(~nearFRLogical).suppressionDurationTimeAll],'Tail','both'); % not sig

% (6) Recovery Rate 
recoveryRate_evokedFR = mean([SPESstruct(nearFRLogical).recoveryRate_Time], "omitnan");
recoveryRate_no_evokedFR = mean([SPESstruct(~nearFRLogical).recoveryRate_Time], "omitnan");
recoveryRate_evokedFRz = mean([SPESstruct(nearFRLogical).recoveryRate_Timez], "omitnan");
recoveryRate_no_evokedFRz = mean([SPESstruct(~nearFRLogical).recoveryRate_Timez], "omitnan");

[p,h,stats] = ranksum([SPESstruct(nearFRLogical).recoveryRate_Time],[SPESstruct(~nearFRLogical).recoveryRate_Time],'Tail','both'); % sig....

% WITHIN UNITS.....
% We have 32 units that exhibited near evoked FR increases (within ~40mm of the stimulation site)
% Look at how these units differ, based on suppression withing the evoked
% firing threshold, and out of the evoked firign threshold

% units that have evoked firing below 40mm
evokedFR_below_Duration = [SPESstruct(nearFRLogical).duration_FRBelow]; 
evokedFR_above_Duration = [SPESstruct(nearFRLogical).duration_FRAbove]; 

[p,h,stats] = ranksum([SPESstruct(nearFRLogical).duration_FRBelow],[SPESstruct(nearFRLogical).duration_FRAbove],'Tail','both')

% getting distance indices
for k = 1:length(SPESstruct)
    sorted_EucIdx =  [SPESstruct(k).sortEucIdx];
    sortedEuc_distance =  [SPESstruct(k).sortedEuc_distance];
    closestDistance_Logical = sortedEuc_distance == min(sortedEuc_distance);
    closeIndex = sorted_EucIdx(closestDistance_Logical); % closest macro stims.
    stimTimes = [SPESstruct(k).stimTimes(closeIndex)];
    unitTimes = [SPESstruct(k).unitTimes];
    spikeWaveforms = [SPESstruct(k).spikeWaveforms]; % spikes by wf samples (48)

    timeWindow = .3; % window of interest to see waveforms after stimulation

    for i = 1:length(stimTimes)
        validTimes_tmp{i} = unitTimes(unitTimes >= stimTimes(i) & unitTimes <= stimTimes(i) + timeWindow) - stimTimes(i);
        validWaveforms_tmp{i} = spikeWaveforms(:,unitTimes >= stimTimes(i) & unitTimes <= stimTimes(i) + timeWindow);
    end
    validTimes{k} = cell2mat(validTimes_tmp');
    tmp = cell2mat(validWaveforms_tmp);
    validWaveforms{k} = cell2mat(validWaveforms_tmp);
    waveformHeight{k} = range(tmp);
    nUnits{k} = k;

    clear validTimes_tmp validWaveforms_tmp
end

% cleaning data: drop empty cells and flatten
validMask = ~cellfun(@isempty, validTimes);  % True for non-empty cells

validTimes_cleaned = validTimes(~cellfun(@isempty, validTimes));
validTimes_flattened = vertcat(validTimes_cleaned{:});
waveformHeight_cleaned = waveformHeight(~cellfun(@isempty, waveformHeight));
waveformHeight_flattened = horzcat(waveformHeight_cleaned{:})';
waveformHeight_flattened = double(waveformHeight_flattened);
validWaveforms_cleaned = validWaveforms(~cellfun(@isempty, validWaveforms));
validWaveforms_flattened = horzcat(validWaveforms_cleaned{:});
nUnits_valid  = nUnits(validMask);
isolatedCluster_valid = isolatedCluster(validMask);

% Expand number of units to match waveforms.....
nUnits_expanded = [];
for i = 1:numel(validTimes_cleaned)
    count = numel(validTimes_cleaned{i});
    nUnits_expanded = [nUnits_expanded; repmat(nUnits_valid(i), count, 1)];
end
nUnits_expanded = cell2mat(nUnits_expanded);

% Expand number of cell classification to match waveforms.....
isolatedCluster_expanded = [];
for i = 1:numel(validTimes_cleaned)
    count = numel(validTimes_cleaned{i});
    isolatedCluster_expanded = [isolatedCluster_expanded; repmat(isolatedCluster_valid(i), count, 1)];
end
%isolatedCluster_expanded = cell2mat(isolatedCluster_expanded);
isolatedCluster_Cat_expanded = categorical(isolatedCluster_expanded);
isolatedCluster_expanded_logical = logical(isolatedCluster_expanded);

% flattening/expanding waveforms
allWaveforms = [];
for i = 1:numel(validWaveforms_cleaned)
    waveform = validWaveforms_cleaned{i};
    allWaveforms = [allWaveforms, waveform];
end

% plot of regression and waveforms
figure(1234)
subplot(1,2,1)
cMap = copper(size(allWaveforms,2));
lm = fitlm(validTimes_flattened,waveformHeight_flattened);
h = lm.plot;
h(1).Marker = '.';
h(1).MarkerEdgeColor = 'w';
h(1).MarkerFaceColor = 'w';
h(1).MarkerSize = 30;
hold on 
[a,ind] = sort(validTimes_flattened);
s = scatter(a,waveformHeight_flattened(ind), 70,cMap, 'filled');
xlabel('spike times (s)');
ylabel('waveform height');
subtitle(sprintf('p = %2f', lm.ModelFitVsNullModel.Pvalue));
legend off
axis tight square
% plot waveforms
subplot(1,2,2)
hold on
cMap(:,4) = 0.5;
l = plot(allWaveforms(:,ind),'LineWidth',0.1);
xlabel('wf samples')
ylabel('waveform amplitude (uV)')
set(l, {'color'},num2cell(cMap,2));
axis tight square
% saving plot
saveas(1234,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\',sprintf('waveform_change_stimulation_3.pdf')))

lm2 = fitlm(validTimes_flattened(isolatedCluster_expanded_logical),waveformHeight_flattened(isolatedCluster_expanded_logical));
lm3 = fitlm(validTimes_flattened(~isolatedCluster_expanded_logical),waveformHeight_flattened(~isolatedCluster_expanded_logical));

% mixed effects model:
wf_tbl = table(validTimes_flattened, waveformHeight_flattened,nUnits_expanded,isolatedCluster_Cat_expanded,...
    'VariableNames', {'validTimes', 'WaveformHeight', 'number_Units', 'cell_Type'});

wf_lme = fitlme(wf_tbl, 'validTimes ~ WaveformHeight + (1|number_Units) + (1|cell_Type)');

% plot: suppression onset latency vs distance..

Minimum_SuppressionAll_Near = ([SPESstruct(:).Minimum_SuppressionAll_Near]);
Minimum_SuppressionAll_Far = ([SPESstruct(:).Minimum_SuppressionAll_Far]);

suppLatency_Near = mean([SPESstruct(:).suppLatency_Near])
suppLatency_Near_std = std([SPESstruct(:).suppLatency_Near])

suppLatency_Far = mean([SPESstruct(:).suppLatency_Far])
suppLatency_Far_std = std([SPESstruct(:).suppLatency_Far])

[p,h,stats] = ranksum([SPESstruct(:).suppLatency_Near],[SPESstruct(:).suppLatency_Far],'Tail','both'); % sig....

%% ~~~~~ Euclidean Distance - Sigmodial Model Checks ~~~~ %%
sigmod_rsquare = [SPESstruct(:).rsquare]
sigmod_distance= [SPESstruct(:).sigmod_distance]
sigmod_nearFRLogical = nearFRLogical; % create grouping 
sigmod_evokedFRLogical = evokedFRLogical;

% matrix for scatter plots
sigmod_matrix = [sigmod_rsquare', sigmod_distance', sigmod_nearFRLogical, sigmod_evokedFRLogical]
sigmod_matrix(sigmod_matrix(:,2)>150,:)=[] % drop rows for distances > 150mm
sigmod_matrix(sigmod_matrix(:,2)<0,:)=[] % drop rows for distances < 0mm

sigmod_nearFRLogical = logical(sigmod_matrix(:,3))
sigmod_evokedFRLogical = logical(sigmod_matrix(:,4))

% Extract data based on the logical index
x = sigmod_matrix(:, 1); % X data
y = sigmod_matrix(:, 2); % Y data

% Create a figure
nBins = 10;
figure(456)
subplot(3,2,1)
scatter(x, y, 'k');
title('R^2 vs. Distance')
xlabel('RSquare');
ylabel('Distance');
axis square
subplot(3,2,2)
histogram(sigmod_matrix(:,1), nBins, 'FaceColor', 'k') % rsquare histogram
title('All Units R^2 Histogram')
axis square
subplot(3,2,3)
scatter(x(~sigmod_nearFRLogical), y(~sigmod_nearFRLogical), 'k', 'DisplayName', 'False');
hold on
scatter(x(sigmod_nearFRLogical), y(sigmod_nearFRLogical), 'MarkerFaceColor', [0.9290 0.6940 0.1250], 'DisplayName', 'True');
title('Near Evoked R^2 vs. Distance')
xlabel('RSquare');
ylabel('Distance');
hold off;
axis square
subplot(3,2,4)
histogram(sigmod_matrix(sigmod_nearFRLogical,1), nBins, 'FaceColor', [0.9290 0.6940 0.1250]) % rsquare histogram
title('Near Evoked Units R^2 Histogram')
axis square
subplot(3,2,5)
scatter(x(~sigmod_evokedFRLogical), y(~sigmod_evokedFRLogical), 'k', 'DisplayName', 'False');
hold on
scatter(x(sigmod_evokedFRLogical), y(sigmod_evokedFRLogical), 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'DisplayName', 'True');
title('All Evoked R^2 vs. Distance')
xlabel('RSquare');
ylabel('Distance');
hold off;
axis square
subplot(3,2,6)
histogram(sigmod_matrix(sigmod_evokedFRLogical,1), nBins, 'FaceColor', [0.8500 0.3250 0.0980]) % rsquare histogram
title('All Evoked Units R^2 Histogram')
axis square
% saving figure 
saveas(456,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\',sprintf('Sigmodial_Rsqr_Distance_ScatterHist.pdf')))
close(456)

% MANUSCRIPT FOGURE TO ADD TO EUCLIDEAN DISTANCES PLOT %
figure(457)
subplot(1,3,1)
scatter(x(~sigmod_nearFRLogical), y(~sigmod_nearFRLogical), 'k', 'DisplayName', 'False');
hold on
scatter(x(sigmod_nearFRLogical), y(sigmod_nearFRLogical), 'MarkerFaceColor', [0.9290 0.6940 0.1250], 'DisplayName', 'True');
title('Near Evoked R^2 vs. Distance')
xlabel('RSquare');
ylabel('Distance');
hold off;
axis square
subplot(1,3,2)
histogram(sigmod_matrix(sigmod_nearFRLogical,1), nBins, 'FaceColor', [0.9290 0.6940 0.1250]) % rsquare histogram
title('Near Evoked Units R^2 Histogram')
axis square
subplot(1,3,3)
histogram(sigmod_matrix(sigmod_nearFRLogical,2), nBins, 'FaceColor', [0.8500 0.3250 0.0980]) % rsquare histogram
title('Near Evoked Units Distance Histogram')
axis square
% saving figure 
saveas(457,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\',sprintf('Sigmodial_Rsqr_Distance_ScatterHist.pdf')))
close(457)

% mean distance and rsquare from sigmodal.
mean(sigmod_matrix(sigmod_nearFRLogical, 2)) % calcualting the average distance for near evoked responses. 
std(sigmod_matrix(sigmod_nearFRLogical, 2))
mean(sigmod_matrix(sigmod_nearFRLogical, 1)) % calcualting the average distance for near evoked responses. 
std(sigmod_matrix(sigmod_nearFRLogical, 1))
mean(sigmod_matrix(sigmod_evokedFRLogical, 2)) % calcualting the average distance for near evoked responses. 
std(sigmod_matrix(sigmod_evokedFRLogical, 2))
mean(sigmod_matrix(sigmod_evokedFRLogical, 1)) % calcualting the average distance for near evoked responses. 
std(sigmod_matrix(sigmod_evokedFRLogical, 1))

%% ~~~~~~~~~~~~~~~~~ Create pie chart for stimulation electrodes  ~~~~~~~~~~~~~~~~~ %%

elecdata = readtable('\\155.100.91.44\D\Data\Rhiannon\CCEPS\SPES_UNITS\ElectrodeInfo\AllStimElectrodeInfo.csv'); % load data
elecdataStruct = table2struct(elecdata)

[counts, groupnames] = groupcounts([elecdata.Region]) % get unique counts
sum(counts); % sanity check (see if counts sums to elecdata size)
filter_counts = counts >= 20; % visual adjustment for pie chart.

hexcolors = [1.0000 0.3412 0.2000; 1.0000 0.7647 0.0000; 0.7804 0.0000 0.2235; 0.5647 0.0471 0.2471; 0.3451 0.0941 0.2784;...
 1.0000 0.3882 0.2784; 1.0000 0.6275 0.4784; 1.0000 0.8431 0.0000; 1.0000 0 0; 1.0000 0.0784 0.5765; 1.0000 0.4118 0.7059; 1.0000 0.2706 0;
 1.0000 0.5490 0; 1.0000 0.8431 0.0000; 0.0000 1.0000 0; 0.1961 0.8039 0.1961; 0.0000 0.5020 0; 0 1.0000 1.0000; 0 0.8078 0.8196; 0.2745 0.5098 0.7059; 0 0 1.0000;
 0.5412 0.1686 0.8863; 0.6 0.1961 0.8; 0.5019 0 0.5019; 1 0 1; 1 0.4118 0.7059; 1.0000 0.0784 0.5765; 0.9333 0.5098 0.9333; 0.5019 0 0; 0.5020 0.5020 0.5020; 0 0 0];
% defined from stim colors in plot
pie(counts, groupnames) % remove counts less than 20. 

pie(counts(filter_counts), groupnames(filter_counts)) % remove counts less than 20. 
colormap(hexcolors)

for i = 1:length(counts)
    proportion(i,:) = counts(i,1)/sum(counts)*100; % getting props for each region
end
proportion = round(proportion)

% stacked barchart
bar(1,proportion,'stacked')
ylim([0 100])
colormap(hexcolors)
legend(groupnames,  "Location", "southeast")

proportion = num2cell(proportion); % changing to cell

% TreeMap of stimulation locations
figure(12)
regionStimProps = [groupnames, proportion]; 
colors = hexcolors;
rectangles = treemap([regionStimProps{:,2}]);
labels = regionStimProps(:,1);
cla
plotRectangles(rectangles,labels,colors)
outline(rectangles)
axis([-0.01 1.01 -0.01 1.01])
title('Stimulation Electrodes')
maximize(figure(12))
% saving figure 
saveas(12,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\WaveFormClassification\',sprintf('stimulationSites_TreeMap.pdf')))
close(12)

% fig2plotly(fig)
% % waffle chart
% I = magic(10);
% %generate where each text will go
% [X Y]=meshgrid(1:10,1:10);
% %create the list of text
% string = mat2cell(num2str([1:10*10]'),ones(10*10,1));
% %display the background
% imagesc(I)

 %% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ UNIT MODULATION FIGURES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

 % Define custom colors
sage_green = [119, 171, 86] / 255; % RGB values for sage green
muted_purple = [79, 49, 170] / 255; % RGB values for muted purple

%% MANUSCRIPT FIGURE 2: CELL TYPE FEATURES %%
clear xlim

% CELL TYPE FEATURES
meanFeatures = [mean_suppAmp_IN, mean_suppLatency_IN; mean_suppAmp_PC, mean_suppLatency_PC]; % need to add duration as a third variable here.
sdFeatures = [std_suppAmp_IN,std_suppLatency_IN; std_suppAmp_PC,std_suppLatency_PC]; 
keyboard

% help set axes the same size.
foo1 = randi(255,20,30,'uint8');
foo2 = randi(255,25,30,'uint8');
subplot(1,2,1); image(foo1)
subplot(1,2,2); image(foo2)

% ~~~~~~~~ (3) Significant Cell Type Units ~~~~~~~~~ %
figure(333)
% (1) FR modulation: Suppression (RdiffAllz)
subplot(2,3,1)
hold on 
scatter(ones(sum(isolatedClusterSigSupp),1),[sigSPESstruct(isolatedClusterSigSupp).BLfreqAll], 'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', sage_green)
scatter(2*ones(sum(isolatedClusterSigSupp),1),[sigSPESstruct(isolatedClusterSigSupp).Minimum_SuppressionAll],'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', sage_green)
scatter(3*ones(sum(~isolatedClusterSigSupp),1),[sigSPESstruct(~isolatedClusterSigSupp).BLfreqAll], 'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', muted_purple)
scatter(4*ones(sum(~isolatedClusterSigSupp),1),[sigSPESstruct(~isolatedClusterSigSupp).Minimum_SuppressionAll],'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', muted_purple)
line([ones(sum(isolatedClusterSigSupp),1) 2*ones(sum(isolatedClusterSigSupp),1)]', [[sigSPESstruct(isolatedClusterSigSupp).BLfreqAll]; [sigSPESstruct(isolatedClusterSigSupp).Minimum_SuppressionAll]], 'Color', sage_green)
line([3*ones(sum(~isolatedClusterSigSupp),1) 4*ones(sum(~isolatedClusterSigSupp),1)]', [[sigSPESstruct(~isolatedClusterSigSupp).BLfreqAll]; [sigSPESstruct(~isolatedClusterSigSupp).Minimum_SuppressionAll]], 'Color', muted_purple)
xlim([0.5 4.5])
subtitle("Suppressed",'FontWeight', 'bold')
axis square
hold off
ylabel('pre-post average ^ FR (spks/s)', 'FontSize', 12); % Y-axis label
xticks([1 2 3 4])
xticklabels({'IN Pre-Stim', 'IN Post-Stim', 'PC Pre-Stim', 'PC Post-Stim'}); % X-axis label , 'FontSize', 14, 'FontWeight', 'bold'
axis square
% (2) FR modulation: Excitation (RdiffAll)
subplot(2,3,2)
hold on 
scatter(ones(sum(isolatedClusterSigExci),1),[sigsexci_SPESstruct(isolatedClusterSigExci).BLfreqAll], 'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', sage_green)
scatter(2*ones(sum(isolatedClusterSigExci),1),[sigsexci_SPESstruct(isolatedClusterSigExci).Maximum_FRAll], 'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', sage_green)
scatter(3*ones(sum(~isolatedClusterSigExci),1),[sigsexci_SPESstruct(~isolatedClusterSigExci).BLfreqAll], 'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', muted_purple)
scatter(4*ones(sum(~isolatedClusterSigExci),1),[sigsexci_SPESstruct(~isolatedClusterSigExci).Maximum_FRAll],'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', muted_purple)
line([ones(sum(isolatedClusterSigExci),1) 2*ones(sum(isolatedClusterSigExci),1)]', [[sigsexci_SPESstruct(isolatedClusterSigExci).BLfreqAll]; [sigsexci_SPESstruct(isolatedClusterSigExci).Maximum_FRAll]], 'Color', sage_green)
line([3*ones(sum(~isolatedClusterSigExci),1) 4*ones(sum(~isolatedClusterSigExci),1)]', [[sigsexci_SPESstruct(~isolatedClusterSigExci).BLfreqAll]; [sigsexci_SPESstruct(~isolatedClusterSigExci).Maximum_FRAll]], 'Color', muted_purple)
xlim([0.5 4.5])
subtitle("Enhanced",'FontWeight', 'bold')
axis square
hold off
ylabel('pre-post average ^ FR (spks/s)', 'FontSize', 12); % Y-axis label
xticks([1 2 3 4])
xticklabels({'IN Pre-Stim', 'IN Post-Stim','PC Pre-Stim', 'PC Post-Stim'}); % X-axis label , 'FontSize', 14, 'FontWeight', 'bold'
axis square
% (3) Suppression Amplitude
SuppressionAmplitude_IN = [sigsupp_SPESstruct(isolatedClusterSigSupp).Minimum_SuppressionAll];
SuppressionAmplitude_INbox = SuppressionAmplitude_IN(:); % needed to be column vector for box.
SuppressionAmplitude_PC = [sigsupp_SPESstruct(~isolatedClusterSigSupp).Minimum_SuppressionAll];
SuppressionAmplitude_PCbox = SuppressionAmplitude_PC(:); % needed to be column vector for box.
jitter_amount = 0.1; % Adjust as needed
jitter_suppressIN = jitter_amount * randn(size(SuppressionAmplitude_IN));
jitter_suppressPC = jitter_amount * randn(size(SuppressionAmplitude_PC));
% Create the subplot
subplot(2,3,3)
hold on;
% Scatter plot
scatter(ones(size(SuppressionAmplitude_IN)) + jitter_suppressIN, SuppressionAmplitude_IN, 50, sage_green, 'filled', '^', 'MarkerEdgeColor', 'k');
scatter(2*ones(size(SuppressionAmplitude_PC)) + jitter_suppressPC, SuppressionAmplitude_PC, 50, muted_purple, 'filled', '^', 'MarkerEdgeColor', 'k');
% Boxplot
boxplot(SuppressionAmplitude_INbox, 'colors', sage_green, 'Positions', 1); % 
boxplot(SuppressionAmplitude_PCbox, 'colors', muted_purple, 'Positions', 2);
% Customizing x-axis and y-axis
xticks([1 2]);
xticklabels({'Interneurons','Principal Cells'});
ylabel('Suppression Amplitude (spks/s)', 'FontSize', 12);
title('Suppression Amplitude');
xlim([0.5, 2.5]);
ylim([0, 80]);
axis square;
hold off;
% (4) Suppression Latency
SuppressionLatency_IN = [sigsupp_SPESstruct(isolatedClusterSigSupp).suppLatencyAll];
SuppressionLatency_INbox = SuppressionLatency_IN(:);
SuppressionLatency_PC = [sigsupp_SPESstruct(~isolatedClusterSigSupp).suppLatencyAll];
SuppressionLatency_PCbox = SuppressionLatency_PC(:);
jitter_amount = 0.1; % Adjust as needed
jitter_latencyIN = jitter_amount * randn(size(SuppressionLatency_IN));
jitter_latencyPC = jitter_amount * randn(size(SuppressionLatency_PC));
% Create the subplot
subplot(2,3,4)
hold on;
% Boxplot
boxplot(SuppressionLatency_INbox, 'colors', sage_green, 'orientation','horizontal', 'Colors', sage_green, 'Positions',1);
boxplot(SuppressionLatency_PCbox, 'colors', muted_purple, 'orientation','horizontal', 'Colors', muted_purple,'Positions',2);
% Scatter plot
scatter(SuppressionLatency_IN, ones(size(SuppressionLatency_IN)) + jitter_latencyIN, 50, sage_green, 'filled', '^', 'MarkerEdgeColor', 'k');
scatter(SuppressionLatency_PC, 2*ones(size(SuppressionLatency_PC)) + jitter_latencyPC, 50, muted_purple, 'filled', '^', 'MarkerEdgeColor', 'k');
% Customizing x-axis and y-axis
xlabel('Suppression Latency (s)', 'FontSize', 12);
ylabel('Cell Types', "Rotation",90,'FontSize', 12);
%ytickangle(90)
title('Suppression Latency', 'FontSize', 14);
ylim([0.5, 2.5]);
xlim([0 2.6])
xticks(0:0.2:2.6)
yticks([1 2]);
axis square
hold off;
% (5) Suppression Duration (s)
SuppressionDuration_IN = [sigsupp_SPESstruct(isolatedClusterSigSupp).suppressionDurationTimeAll];
SuppressionDuration_INbox = SuppressionDuration_IN(:); % needed to be column vector for box.
SuppressionDuration_PC = [sigsupp_SPESstruct(~isolatedClusterSigSupp).suppressionDurationTimeAll];
SuppressionDuration_PCbox = SuppressionDuration_PC(:); % needed to be column vector for box.
jitter_durationIN = jitter_amount * randn(size(SuppressionDuration_IN));
jitter_durationPC = jitter_amount * randn(size(SuppressionDuration_PC));
% Create the subplot
subplot(2,3,5)
hold on;
% Boxplot
boxplot(SuppressionDuration_INbox, 'colors', sage_green, 'orientation','horizontal', 'Colors', sage_green, 'Positions',1);
boxplot(SuppressionDuration_PCbox, 'colors', muted_purple, 'orientation','horizontal', 'Colors', muted_purple,'Positions',2);
% Scatter plot
scatter(SuppressionDuration_IN, ones(size(SuppressionDuration_IN)) + jitter_durationIN, 50, sage_green, 'filled', '^', 'MarkerEdgeColor', 'k');
scatter(SuppressionDuration_PC, 2*ones(size(SuppressionDuration_PC)) + jitter_durationPC, 50, muted_purple, 'filled', '^', 'MarkerEdgeColor', 'k');
% Customizing x-axis and y-axis
xlabel('Suppression Duration (s)', 'FontSize', 12);
ylabel('Cell Types', "Rotation",90,'FontSize', 12);
%ytickangle(90)
title('Suppression Duration', 'FontSize', 14);
ylim([0.5, 2.5]);
xlim([0 2.6])
yticks([1 2]);
xticks(0:0.2:2.6)
axis square
hold off
% (5) Area Under The Curve 
AUC_IN = [SPESstruct(isolatedCluster).AUC];
AUC_IN(AUC_IN==0) = [];
AUC_INbox = AUC_IN(:);
AUC_PC = [SPESstruct(~isolatedCluster).AUC];
AUC_PC(AUC_PC==0) = [];
AUC_PCbox = AUC_PC(:);
jitter_amount = 0.1; % Adjust as needed
jitter_AUCIN = jitter_amount * randn(size(AUC_IN));
jitter_AUCPC = jitter_amount * randn(size(AUC_PC));
% Create the subplot
subplot(2,3,6)
hold on;
% Scatter plot
scatter(ones(size(AUC_IN)) + jitter_AUCIN, AUC_IN, 50, sage_green, 'filled', '^', 'MarkerEdgeColor', 'k');
scatter(2*ones(size(AUC_PC)) + jitter_AUCPC, AUC_PC, 50, muted_purple, 'filled', '^', 'MarkerEdgeColor', 'k');
% Boxplot
boxplot(AUC_INbox, 'colors', sage_green, 'Positions', 1);
boxplot(AUC_PCbox, 'colors', muted_purple, 'Positions', 2);
% Customizing x-axis and y-axis
xticks([1 2]);
xticklabels({'Interneurons','Principal Cells'});
ylabel('Area Under the Curve (spks)', 'FontSize', 12);
title('AUC');
xlim([0.5, 2.5]);
yticks(0:2:24);
ylim([0 24])
axis square;
hold off;

maximize(333)

saveas(334,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\Unit_Modulation',sprintf('CellType_Modulation_Characteristics_AcrossSubs.pdf')))
close(333)

% ~~~~~~~~ (4) Modulation of Significant SOZ Cells ~~~~~~~~~ %

figure(444)
% (1) FR modulation: Suppression (RdiffAllz)
subplot(2,3,1)
hold on 
scatter(ones(sum(SOZsigSupp),1),[sigSPESstruct(SOZsigSupp).BLfreqAll], 'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'yellow')
scatter(2*ones(sum(SOZsigSupp),1),[sigSPESstruct(SOZsigSupp).Minimum_SuppressionAll],'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'yellow')
scatter(3*ones(sum(~SOZsigSupp),1),[sigSPESstruct(~SOZsigSupp).BLfreqAll], 'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
scatter(4*ones(sum(~SOZsigSupp),1),[sigSPESstruct(~SOZsigSupp).Minimum_SuppressionAll],'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
line([ones(sum(SOZsigSupp),1) 2*ones(sum(SOZsigSupp),1)]', [[sigSPESstruct(SOZsigSupp).BLfreqAll]; [sigSPESstruct(SOZsigSupp).Minimum_SuppressionAll]], 'Color', 'yellow')
line([3*ones(sum(~SOZsigSupp),1) 4*ones(sum(~SOZsigSupp),1)]', [[sigSPESstruct(~SOZsigSupp).BLfreqAll]; [sigSPESstruct(~SOZsigSupp).Minimum_SuppressionAll]], 'Color', 'k')
xlim([0.5 4.5])
subtitle("Suppressed",'FontWeight', 'bold')
axis square
hold off
ylabel('pre-post average ^ FR (Hz)', 'FontSize', 12); % Y-axis label
xticks([1 2 3 4])
xticklabels({'SOZ Pre-Stim', 'SOZ Post-Stim', 'nonSOZ Pre-Stim', 'nonSOZ Post-Stim'}); % X-axis label , 'FontSize', 14, 'FontWeight', 'bold'
axis square
% (2) FR modulation: Excitation (RdiffAll)
subplot(2,3,2)
hold on 
scatter(ones(sum(SOZsigExci),1),[sigSPESstruct(SOZsigExci).BLfreqAll], 'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'yellow')
scatter(2*ones(sum(SOZsigExci),1),[sigSPESstruct(SOZsigExci).Maximum_FRAll], 'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'yellow')
scatter(3*ones(sum(~SOZsigExci),1),[sigSPESstruct(~SOZsigExci).BLfreqAll], 'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
scatter(4*ones(sum(~SOZsigExci),1),[sigSPESstruct(~SOZsigExci).Maximum_FRAll],'filled', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
line([ones(sum(SOZsigExci),1) 2*ones(sum(SOZsigExci),1)]', [[sigSPESstruct(SOZsigExci).BLfreqAll]; [sigSPESstruct(SOZsigExci).Maximum_FRAll]], 'Color', 'yellow')
line([3*ones(sum(~SOZsigExci),1) 4*ones(sum(~SOZsigExci),1)]', [[sigSPESstruct(~SOZsigExci).BLfreqAll]; [sigSPESstruct(~SOZsigExci).Maximum_FRAll]], 'Color', 'k')
xlim([0.5 4.5])
subtitle("Enhanced",'FontWeight', 'bold')
axis square
hold off
ylabel('pre-post average ^ FR (spks/s)', 'FontSize', 12); % Y-axis label
xticks([1 2 3 4])
xticklabels({'SOZ Pre-Stim', 'SOZ Post-Stim','nonSOZ Pre-Stim', 'nonSOZ Post-Stim'}); % X-axis label , 'FontSize', 14, 'FontWeight', 'bold'
axis square
% (3) Suppression Amplitude
SuppressionAmplitude_SOZ = [sigsupp_SPESstruct(SOZsigSupp).Minimum_SuppressionAll];
SuppressionAmplitude_SOZbox = SuppressionAmplitude_SOZ(:); % needed to be column vector for box.
SuppressionAmplitude_nonSOZ = [sigsupp_SPESstruct(~SOZsigSupp).Minimum_SuppressionAll];
SuppressionAmplitude_nonSOZbox = SuppressionAmplitude_nonSOZ(:); % needed to be column vector for box.
jitter_amount = 0.1; % Adjust as needed
jitter_suppressSOZ = jitter_amount * randn(size(SuppressionAmplitude_SOZ));
jitter_suppressnonSOZ = jitter_amount * randn(size(SuppressionAmplitude_nonSOZ));
% Create the subplot
subplot(2,3,3)
hold on;
% Scatter plot
scatter(ones(size(SuppressionAmplitude_SOZ)) + jitter_suppressSOZ, SuppressionAmplitude_SOZ, 50, 'yellow', 'filled', '^', 'MarkerEdgeColor', 'k');
scatter(2*ones(size(SuppressionAmplitude_nonSOZ)) + jitter_suppressnonSOZ, SuppressionAmplitude_nonSOZ, 50, 'b', 'filled', '^', 'MarkerEdgeColor', 'k');
% Boxplot
boxplot(SuppressionAmplitude_SOZbox, 'Colors', 'k', 'Positions', 1); % 
boxplot(SuppressionAmplitude_nonSOZbox, 'Colors', 'k', 'Positions', 2);
% Customizing x-axis and y-axis
xticks([1 2]);
xticklabels({'SOZ','nonSOZ'});
ylabel('Suppression Amplitude (spks/s)', 'FontSize', 12);
title('Suppression Amplitude');
xlim([0.5, 2.5]);
%ylim([0, 80]);
axis square;
hold off;
% (4) Suppression Latency
SuppressionLatency_SOZ = [sigsupp_SPESstruct(SOZsigSupp).suppLatencyAll];
SuppressionLatency_SOZbox = SuppressionLatency_SOZ(:);
SuppressionLatency_nonSOZ = [sigsupp_SPESstruct(~SOZsigSupp).suppLatencyAll];
SuppressionLatency_nonSOZbox = SuppressionLatency_nonSOZ(:);
jitter_amount = 0.1; % Adjust as needed
jitter_latencySOZ = jitter_amount * randn(size(SuppressionLatency_SOZ));
jitter_latencynonSOZ = jitter_amount * randn(size(SuppressionLatency_nonSOZ));
% Create the subplot
subplot(2,3,4)
hold on;
% Boxplot
boxplot(SuppressionLatency_SOZbox, 'Colors', 'k', 'orientation','horizontal', 'Positions',1);
boxplot(SuppressionLatency_nonSOZbox, 'Colors', 'k', 'orientation','horizontal', 'Positions',2);
% Scatter plot
scatter(SuppressionLatency_SOZ, ones(size(SuppressionLatency_SOZ)) + jitter_latencySOZ, 50, 'yellow', 'filled', '^', 'MarkerEdgeColor', 'k');
scatter(SuppressionLatency_nonSOZ, 2*ones(size(SuppressionLatency_nonSOZ)) + jitter_latencynonSOZ, 50, 'b', 'filled', '^', 'MarkerEdgeColor', 'k');
% Customizing x-axis and y-axis
xlabel('Suppression Latency (s)', 'FontSize', 12);
ylabel('Cell Types', "Rotation",90,'FontSize', 12);
%ytickangle(90)
title('Suppression Latency', 'FontSize', 14);
ylim([0.5, 2.5]);
xlim([0 2.6])
yticks([1 2]);
xticks(0:0.2:2.6);
axis square
hold off
% (5) Suppression Duration (s)
SuppressionDuration_SOZ = [sigsupp_SPESstruct(SOZsigSupp).suppressionDurationTimeAll];
SuppressionDuration_SOZbox = SuppressionDuration_SOZ(:); % needed to be column vector for box.
SuppressionDuration_nonSOZ = [sigsupp_SPESstruct(~SOZsigSupp).suppressionDurationTimeAll];
SuppressionDuration_nonSOZbox = SuppressionDuration_nonSOZ(:); % needed to be column vector for box.
jitter_durationSOZ = jitter_amount * randn(size(SuppressionDuration_SOZ));
jitter_durationnonSOZ = jitter_amount * randn(size(SuppressionDuration_nonSOZ));
% Create the subplot
subplot(2,3,5)
hold on;
% Boxplot
boxplot(SuppressionDuration_SOZbox, 'Colors', 'k', 'orientation','horizontal', 'Positions',1);
boxplot(SuppressionDuration_nonSOZbox, 'Colors', 'k', 'orientation','horizontal', 'Positions',2);
% Scatter plot
scatter(SuppressionDuration_SOZ, ones(size(SuppressionDuration_SOZ)) + jitter_durationSOZ, 50, 'yellow', 'filled', '^', 'MarkerEdgeColor', 'k');
scatter(SuppressionDuration_nonSOZ, 2*ones(size(SuppressionDuration_nonSOZ)) + jitter_durationnonSOZ, 50, 'b', 'filled', '^', 'MarkerEdgeColor', 'k');
% Customizing x-axis and y-axis
xlabel('Suppression Duration (s)', 'FontSize', 12);
ylabel('Cell Types', "Rotation",90,'FontSize', 12);
%ytickangle(90)
title('Suppression Duration', 'FontSize', 14);
ylim([0.5, 2.5]);
xlim([0  2.6]);
yticks([1 2]);
xticks(0:0.2:2.6);
axis square
hold off
% (6) Area Under The Curve 
AUC_SOZ = [SPESstruct(isolatedCluster).AUC];
AUC_SOZ(AUC_SOZ==0) = [];
AUC_SOZbox = AUC_SOZ(:);
AUC_nonSOZ = [SPESstruct(~isolatedCluster).AUC];
AUC_nonSOZ(AUC_nonSOZ==0) = [];
AUC_nonSOZbox = AUC_nonSOZ(:);
jitter_amount = 0.1; % Adjust as needed
jitter_AUCSOZ = jitter_amount * randn(size(AUC_SOZ));
jitter_AUCnonSOZ = jitter_amount * randn(size(AUC_nonSOZ));
% Create the subplot
subplot(2,3,6)
hold on;
% Scatter plot
scatter(ones(size(AUC_SOZ)) + jitter_AUCSOZ, AUC_SOZ, 50, 'yellow', 'filled', '^', 'MarkerEdgeColor', 'k');
scatter(2*ones(size(AUC_nonSOZ)) + jitter_AUCnonSOZ, AUC_nonSOZ, 50, 'b', 'filled', '^', 'MarkerEdgeColor', 'k');
% Boxplot
boxplot(AUC_SOZbox, 'colors', 'k', 'Positions', 1);
boxplot(AUC_nonSOZbox, 'colors', 'k', 'Positions', 2);
% Customizing x-axis and y-axis
xticks([1 2]);
xticklabels({'SOZ','nonSOZ'});
ylabel('Area Under the Curve (Hz)', 'FontSize', 12);
title('AUC');
xlim([0.5, 2.5]);
yticks(0:2:24);
ylim([0 24])
axis square;
hold off;
maximize(444)

saveas(444,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\Unit_Modulation',sprintf('SOZ_Modulation_Characteristics_AcrossSubs.pdf')))
close(444)

%% ~~~~~~~~~~~~~~~~~ ANATOMICAL LOCATIONS OF RESPONSES ~~~~~~~~~~~~~~~~~~ %%

SigAll = [SPESstruct(:).hAllz] == 1;
SigFar = [SPESstruct(:).hFarz] == 1;
SigNear = [SPESstruct(:).hNearz] == 1;

unitLocs = {SPESstruct(:).chanLabel}; % get unit locations
unitLocs = string(unitLocs); % changing to string
SigStimsAllLocs = {unitLocs; SigFar; SigNear}; % adding to matrix with Far Near and Instant stims.
SigStimsFarLocs = {unitLocs; SigFar};
SigStimsNearLocs = {unitLocs; SigNear};
%SigStimsInstantLocs = {unitLocs; SigInstant};

unitLocs(SigNear); % creates logical index of the contact location of where the significant responses happened
[uniqueLocsNear,~,IcNear] = unique(unitLocs(SigNear)); % creates [string of all individual sigcontact locations, ~, # given to each location].
[cnt_uniqueNear, unique_unitLocsNear] = histcounts(IcNear); % [# of XXXX, # of unique locations] % what is count unique?

unitLocs(SigFar);
[uniqueLocsFar,~,IcFar] = unique(unitLocs(SigFar));
[cnt_uniqueFar, unique_unitLocsFar] = histcounts(IcFar);

%% Need to group locations
% ~~~~~~~~~~~~~~~~~~~~ (0) all locations ~~~~~~~~~~~~~~~~~~~~~~~~%
uniqueUnitLocs = unique(unitLocs,'stable'); % finding the total unique locations that were recorded from.
occurencesUniqueUnitLocs = cellfun(@(x) sum(ismember(unitLocs,x)),uniqueUnitLocs,'un',0); % # of units recorded from each unique location.
occurencesNamesUniqueUnitLocs = [uniqueUnitLocs; occurencesUniqueUnitLocs]; % Combining Locations and # of occurences

y = cellfun(@str2num, occurencesNamesUniqueUnitLocs(2,:)); % creating variables for plot.
x = categorical(occurencesNamesUniqueUnitLocs(1,:));

% KEY REGIONS OF INTEREST
AMY = "Amygdala"; % finding locations with Amygdala
AMYLocsAll = contains(unitLocs(1,:), AMY);
sum(AMYLocsAll)
OFC = ("OFC"| "Orbital" | "vmPFC"); % finding locations with OFC
OFCLocsAll = contains(unitLocs(1,:), OFC);
sum(OFCLocsAll)
HIPP = "Hippocampus"; % finding locations with HIPP
HIPPLocsAll = contains(unitLocs(1,:), HIPP);
sum(HIPPLocsAll)
CING = "Cingulate"; % finding locations with Cingulate
CINGLocsAll = contains(unitLocs(1,:), CING);
sum(CINGLocsAll)
pCING = "Posterior"; % finding locations with Cingulate
pCINGLocsAll = contains(unitLocs(1,:), pCING);
sum(pCINGLocsAll)

% ~~~~~~~~ (1) locations with any significant unit responses ~~~~~~~~~ % 
alllocspositions = find(SigStimsAllLocs{2,:}); % finding locations of responses
SigStimsAllLocs{1}(alllocspositions); % All locations for all unit.
uniqueUnitAllLocsSIG = cellstr(unique(SigStimsAllLocs{1}(alllocspositions),'stable')); % finding the total unique locations that had significant far unit responses.
occurencesUniqueUnitAllLocsSIG = cellfun(@(x) sum(ismember(SigStimsAllLocs{1}(alllocspositions),x)),uniqueUnitAllLocsSIG,'un',0); % # of units recorded from each unique location.
occurencesNamesUniqueUnitAllLocsSIG = vertcat(uniqueUnitAllLocsSIG, occurencesUniqueUnitAllLocsSIG); % made it a cell array

%y1 = cellfun(@str2num, occurencesNamesUniqueUnitAllLocsSIG{2,:})
y1 = cell2mat(occurencesNamesUniqueUnitAllLocsSIG(2,:));
x1 = categorical(occurencesNamesUniqueUnitAllLocsSIG(1,:));

% KEY REGIONS OF INTEREST
% Amygdala
AMYLocs1 = occurencesNamesUniqueUnitAllLocsSIG(2,contains(occurencesNamesUniqueUnitAllLocsSIG(1,:), AMY));
AMYLocs1=cellfun(@sum,AMYLocs1);
sum(AMYLocs1)
% OFC
OFCLocs1 = occurencesNamesUniqueUnitAllLocsSIG(2,contains(occurencesNamesUniqueUnitAllLocsSIG(1,:), OFC));
OFCLocs1=cellfun(@sum,OFCLocs1);
sum(OFCLocs1)
% Hippocampus
HIPPLocs1 = occurencesNamesUniqueUnitAllLocsSIG(2,contains(occurencesNamesUniqueUnitAllLocsSIG(1,:), HIPP));
HIPPLocs1=cellfun(@sum,HIPPLocs1);
sum(HIPPLocs1)
% Cingulate
CINGLocs1 = occurencesNamesUniqueUnitAllLocsSIG(2,contains(occurencesNamesUniqueUnitAllLocsSIG(1,:), CING));
CINGLocs1=cellfun(@sum,CINGLocs1);
sum(CINGLocs1)
% Posterior Cingulate
pCINGLocs1 = occurencesNamesUniqueUnitAllLocsSIG(2,contains(occurencesNamesUniqueUnitAllLocsSIG(1,:), pCING));
pCINGLocs1=cellfun(@sum,pCINGLocs1);
sum(pCINGLocs1)

% (1b) proportion of significant unit responses
propAMYLocs1 = (sum(AMYLocs1)/sum(AMYLocsAll))*100
propOFCLocs1 = (sum(OFCLocs1)/sum(OFCLocsAll))*100
propHIPPLocs1 = (sum(HIPPLocs1)/sum(HIPPLocsAll))*100
propCINGLocs1 = (sum(CINGLocs1)/sum(CINGLocsAll))*100
proppCINGLocs1 = (sum(pCINGLocs1)/sum(pCINGLocsAll))*100

% ~~~~~~~~  (2) locations with Far significant unit responses ~~~~~~~~~ %
farlocspositions = find(SigStimsFarLocs{2,:}); % finding locations of responses
SigStimsFarLocs{1}(farlocspositions); % All locations for each far unit.
uniqueUnitFarLocsSIG = cellstr(unique(SigStimsFarLocs{1}(farlocspositions),'stable')); % finding the total unique locations that had significant far unit responses.
occurencesUniqueUnitFarLocsSIG = cellfun(@(x) sum(ismember(SigStimsFarLocs{1}(farlocspositions),x)),uniqueUnitFarLocsSIG,'un',0); % # of units recorded from each unique location.
occurencesNamesUniqueUnitFarLocsSIG = vertcat(uniqueUnitFarLocsSIG, occurencesUniqueUnitFarLocsSIG);

y2 = cell2mat(occurencesNamesUniqueUnitFarLocsSIG(2,:));
x2 = categorical(occurencesNamesUniqueUnitFarLocsSIG(1,:));

% KEY REGIONS OF INTEREST
% Amygdala
AMYLocs2 = occurencesNamesUniqueUnitFarLocsSIG(2,contains(occurencesNamesUniqueUnitFarLocsSIG(1,:), AMY));
AMYLocs2=cellfun(@sum,AMYLocs2);
sum(AMYLocs2)
% OFC
OFCLocs2 = occurencesNamesUniqueUnitFarLocsSIG(2,contains(occurencesNamesUniqueUnitFarLocsSIG(1,:), OFC));
OFCLocs2=cellfun(@sum,OFCLocs2);
sum(OFCLocs2)
% Hippocampus
HIPPLocs2 = occurencesNamesUniqueUnitFarLocsSIG(2,contains(occurencesNamesUniqueUnitFarLocsSIG(1,:), HIPP));
HIPPLocs2=cellfun(@sum,HIPPLocs2);
sum(HIPPLocs2)
% Cingulate
CINGLocs2 = occurencesNamesUniqueUnitFarLocsSIG(2,contains(occurencesNamesUniqueUnitFarLocsSIG(1,:), CING));
CINGLocs2=cellfun(@sum,CINGLocs2);
sum(CINGLocs2)
% Posterior Cingulate
pCINGLocs2 = occurencesNamesUniqueUnitFarLocsSIG(2,contains(occurencesNamesUniqueUnitFarLocsSIG(1,:), pCING));
pCINGLocs2=cellfun(@sum,pCINGLocs2);
sum(pCINGLocs2)

% (2b) proportion of significant Far unit responses
propAMYLocs2 = (sum(AMYLocs2)/sum(AMYLocsAll))*100
propOFCLocs2 = (sum(OFCLocs2)/sum(OFCLocsAll))*100
propHIPPLocs2 = (sum(HIPPLocs2)/sum(HIPPLocsAll))*100
propCINGLocs2 = (sum(CINGLocs2)/sum(CINGLocsAll))*100
proppCINGLocs2 = (sum(pCINGLocs2)/sum(pCINGLocsAll))*100

% ~~~~~~~~  (3) locations with Near significant unit responses ~~~~~~~~~ %
nearlocspositions = find(SigStimsNearLocs{2,:}); % finding locations of responses
SigStimsNearLocs{1}(nearlocspositions);
uniqueUnitNearLocsSIG = cellstr(unique(SigStimsNearLocs{1}(nearlocspositions),'stable')); % finding the total unique locations that had significant unit responses.
occurencesUniqueUnitNearLocsSIG = cellfun(@(x) sum(ismember(SigStimsNearLocs{1}(nearlocspositions),x)),uniqueUnitNearLocsSIG,'un',0); % # of units recorded from each unique location.
occurencesNamesUniqueUnitNearLocsSIG = vertcat(uniqueUnitNearLocsSIG, occurencesUniqueUnitNearLocsSIG);

% setting plot data
y3 = cell2mat(occurencesNamesUniqueUnitNearLocsSIG(2,:));
x3 = categorical(occurencesNamesUniqueUnitNearLocsSIG(1,:));

% KEY REGIONS OF INTEREST
% Amygdala
AMYLocs3 = occurencesNamesUniqueUnitNearLocsSIG(2,contains(occurencesNamesUniqueUnitNearLocsSIG(1,:), AMY));
AMYLocs3=cellfun(@sum,AMYLocs3);
sum(AMYLocs3)
% OFC
OFCLocs3 = occurencesNamesUniqueUnitNearLocsSIG(2,contains(occurencesNamesUniqueUnitNearLocsSIG(1,:), OFC));
OFCLocs3=cellfun(@sum,OFCLocs3);
sum(OFCLocs3)
% Hippocampus
HIPPLocs3 = occurencesNamesUniqueUnitNearLocsSIG(2,contains(occurencesNamesUniqueUnitNearLocsSIG(1,:), HIPP));
HIPPLocs3=cellfun(@sum,HIPPLocs3);
sum(HIPPLocs3)
% Cingulate
CINGLocs3 = occurencesNamesUniqueUnitNearLocsSIG(2,contains(occurencesNamesUniqueUnitNearLocsSIG(1,:), CING));
CINGLocs3=cellfun(@sum,CINGLocs3);
sum(CINGLocs3)
% posterior Cingulate
pCINGLocs3 = occurencesNamesUniqueUnitNearLocsSIG(2,contains(occurencesNamesUniqueUnitNearLocsSIG(1,:), pCING));
pCINGLocs3=cellfun(@sum,pCINGLocs3);
sum(pCINGLocs3)

% (3b) proportion of significant Near unit responses
propAMYLocs3 = (sum(AMYLocs3)/sum(AMYLocsAll))*100
propOFCLocs3 = (sum(OFCLocs3)/sum(OFCLocsAll))*100
propHIPPLocs3 = (sum(HIPPLocs3)/sum(HIPPLocsAll))*100
propCINGLocs3 = (sum(CINGLocs3)/sum(CINGLocsAll))*100
proppCINGLocs3 = (sum(pCINGLocs3)/sum(pCINGLocsAll))*100

% ~~~~~~~~~~ Figure 1 E. BAR PLOTS FOR GROUPED LOCATIONS and RESPONSES ~~~~~~~~~~ %

% plot proportions of sig unit responses
Regions = ["Amygdala"; "OFC & vmPFC"; "Hippocampus"; "Cingulate"; "Posterior Cingulate"];
x_Regions = categorical(Regions);
Muted_Blue =  [0.4, 0.5, 0.6];
Muted_Green = [0.45, 0.55, 0.45];
Dark_Blue = [0.1, 0.15, 0.3];
Muted_Purple = [0.55, 0.406, 0.6];
Dark_Brown = [0.6588,0.5490,0.4902];
% All Unit Counts by grouped region....
y_unitCount = [sum(AMYLocsAll); sum(OFCLocsAll); sum(HIPPLocsAll); sum(CINGLocsAll); sum(pCINGLocsAll)];
figure(2)
b = bar(x_Regions, y_unitCount,'FaceColor','flat');
% Assign colors to each bar
b.CData(1,:) = Muted_Blue;
b.CData(2,:) = Muted_Green;
b.CData(3,:) = Dark_Blue;
b.CData(4,:) = Muted_Purple;
b.CData(5,:) = Dark_Brown;
title('Regional Unit Counts', 'FontSize', 14) % Increase title font size
xlabel('Regions', 'FontSize', 12) % Increase x-axis label font size
ylabel('Unit Count', 'FontSize', 12) % Increase y-axis label font size
ylim([0 100])
% saving figure
saveas(2,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\',sprintf('UnitCounts_GroupedLocations_AcrossSubs.pdf')))
close(2)

%% ~~~~~~~~~~~~~~~~ REGIONAL SUPPRESSION? ~~~~~~~~~~~~~~~~~~~~~~~ %%

% REGIONAL SCATTERS
%(1) AMYGDALA
% scatter (AMY: amplitude vs latency)
scatter([sigsupp_SPESstruct(sigsuppAMYLocs).Minimum_SuppressionAll], [sigsupp_SPESstruct(sigsuppAMYLocs).suppLatencyAll])
mdl = fitlm([sigsupp_SPESstruct(sigsuppAMYLocs).Minimum_SuppressionAll], [sigsupp_SPESstruct(sigsuppAMYLocs).suppLatencyAll])
plot(mdl)
% scatter (AMY: amplitude vs duration)
scatter([sigsupp_SPESstruct(sigsuppAMYLocs).Minimum_SuppressionAll], [sigsupp_SPESstruct(sigsuppAMYLocs).suppressionDurationTimeAll])
mdl = fitlm([sigsupp_SPESstruct(sigsuppAMYLocs).Minimum_SuppressionAll], [sigsupp_SPESstruct(sigsuppAMYLocs).suppressionDurationTimeAll], 'Exclude', [sigsupp_SPESstruct(sigsuppAMYLocs).suppressionDurationTimeAll] == 0)
plot(mdl)
% scatter (AMY: duration vs latency)
scatter([sigsupp_SPESstruct(sigsuppAMYLocs).suppLatencyAll], [sigsupp_SPESstruct(sigsuppAMYLocs).suppressionDurationTimeAll])
mdl = fitlm([sigsupp_SPESstruct(sigsuppAMYLocs).suppLatencyAll], [sigsupp_SPESstruct(sigsuppAMYLocs).suppressionDurationTimeAll], 'Exclude', [sigsupp_SPESstruct(sigsuppAMYLocs).suppressionDurationTimeAll] == 0)
plot(mdl)

%(2) OFC
% scatter (OFC: amplitude vs latency)
scatter([sigsupp_SPESstruct(sigsuppOFCLocs).Minimum_SuppressionAll], [sigsupp_SPESstruct(sigsuppOFCLocs).suppLatencyAll])
mdl = fitlm([sigsupp_SPESstruct(sigsuppOFCLocs).Minimum_SuppressionAll], [sigsupp_SPESstruct(sigsuppOFCLocs).suppLatencyAll])
plot(mdl)
% scatter (OFC: amplitude vs duration) 
scatter([sigsupp_SPESstruct(sigsuppOFCLocs).Minimum_SuppressionAll], [sigsupp_SPESstruct(sigsuppOFCLocs).suppressionDurationTimeAll])
mdl = fitlm([sigsupp_SPESstruct(sigsuppOFCLocs).Minimum_SuppressionAll], [sigsupp_SPESstruct(sigsuppOFCLocs).suppressionDurationTimeAll], 'Exclude', [sigsupp_SPESstruct(sigsuppOFCLocs).suppressionDurationTimeAll] == 0)
plot(mdl)
% scatter (OFC: duration vs latency)
scatter([sigsupp_SPESstruct(sigsuppOFCLocs).suppLatencyAll], [sigsupp_SPESstruct(sigsuppOFCLocs).suppressionDurationTimeAll])
mdl = fitlm([sigsupp_SPESstruct(sigsuppOFCLocs).suppLatencyAll], [sigsupp_SPESstruct(sigsuppOFCLocs).suppressionDurationTimeAll], 'Exclude', [sigsupp_SPESstruct(sigsuppOFCLocs).suppressionDurationTimeAll] == 0)
plot(mdl)

%(3) HIPP
% scatter (HIPP: amplitude vs latency)
scatter([sigsupp_SPESstruct(sigsuppHIPPLocs).Minimum_SuppressionAll], [sigsupp_SPESstruct(sigsuppHIPPLocs).suppLatencyAll])
mdl = fitlm([sigsupp_SPESstruct(sigsuppHIPPLocs).Minimum_SuppressionAll], [sigsupp_SPESstruct(sigsuppHIPPLocs).suppLatencyAll])
plot(mdl)
% scatter (HIPP: amplitude vs duration) (SIG).
scatter([sigsupp_SPESstruct(sigsuppHIPPLocs).Minimum_SuppressionAll], [sigsupp_SPESstruct(sigsuppHIPPLocs).suppressionDurationTimeAll])
mdl = fitlm([sigsupp_SPESstruct(sigsuppHIPPLocs).Minimum_SuppressionAll], [sigsupp_SPESstruct(sigsuppHIPPLocs).suppressionDurationTimeAll], 'Exclude', [sigsupp_SPESstruct(sigsuppHIPPLocs).suppressionDurationTimeAll] == 0)
plot(mdl)
% scatter (HIPP: duration vs latency) (SIG)
scatter([sigsupp_SPESstruct(sigsuppHIPPLocs).suppLatencyAll], [sigsupp_SPESstruct(sigsuppHIPPLocs).suppressionDurationTimeAll])
mdl = fitlm([sigsupp_SPESstruct(sigsuppHIPPLocs).suppLatencyAll], [sigsupp_SPESstruct(sigsuppHIPPLocs).suppressionDurationTimeAll], 'Exclude', [sigsupp_SPESstruct(sigsuppHIPPLocs).suppressionDurationTimeAll] == 0)
plot(mdl)

%(4) CING
% scatter (CING: amplitude vs latency)
scatter([sigsupp_SPESstruct(sigsuppCINGLocs).Minimum_SuppressionAll], [sigsupp_SPESstruct(sigsuppCINGLocs).suppLatencyAll])
mdl = fitlm([sigsupp_SPESstruct(sigsuppCINGLocs).Minimum_SuppressionAll], [sigsupp_SPESstruct(sigsuppCINGLocs).suppLatencyAll])
plot(mdl)
% scatter (CING: amplitude vs duration) (SIG).
scatter([sigsupp_SPESstruct(sigsuppCINGLocs).Minimum_SuppressionAll], [sigsupp_SPESstruct(sigsuppCINGLocs).suppressionDurationTimeAll])
mdl = fitlm([sigsupp_SPESstruct(sigsuppCINGLocs).Minimum_SuppressionAll], [sigsupp_SPESstruct(sigsuppCINGLocs).suppressionDurationTimeAll], 'Exclude', [sigsupp_SPESstruct(sigsuppCINGLocs).suppressionDurationTimeAll] == 0)
plot(mdl)
% scatter (CING: duration vs latency) (SIG)
scatter([sigsupp_SPESstruct(sigsuppCINGLocs).suppLatencyAll], [sigsupp_SPESstruct(sigsuppCINGLocs).suppressionDurationTimeAll])
mdl = fitlm([sigsupp_SPESstruct(sigsuppCINGLocs).suppLatencyAll], [sigsupp_SPESstruct(sigsuppCINGLocs).suppressionDurationTimeAll], 'Exclude', [sigsupp_SPESstruct(sigsuppCINGLocs).suppressionDurationTimeAll] == 0)
plot(mdl)

%% ~~~~~~~~~~~~~ SUPPLEMENTARY FIGURE: WAVEFORM METRICS BY REGION ~~~~~~~~~~~~


AMYLocsIN_Logical = isolatedCluster & AMYLocs';
AMYLocsPC_Logical = ~isolatedCluster & AMYLocs';
OFCLocsIN_Logical = isolatedCluster & OFCLocs';
OFCLocsPC_Logical = ~isolatedCluster & OFCLocs';
HIPPLocsIN_Logical = isolatedCluster & HIPPLocs';
HIPPLocsPC_Logical = ~isolatedCluster & HIPPLocs';
CINGLocsIN_Logical = isolatedCluster & CINGLocs';
CINGLocsPC_Logical = ~isolatedCluster & CINGLocs';
pCINGLocsIN_Logical = isolatedCluster & pCINGLocs';
pCINGLocsPC_Logical = ~isolatedCluster & pCINGLocs';

% segment data (TTP)
principalCellTTP_AMY = metricTable.TTP(AMYLocsPC_Logical);
interneuronTTP_AMY = metricTable.TTP(AMYLocsIN_Logical);
principalCellTTP_OFC = metricTable.TTP(OFCLocsPC_Logical);
interneuronTTP_OFC = metricTable.TTP(OFCLocsIN_Logical);
principalCellTTP_HIPP = metricTable.TTP(HIPPLocsPC_Logical);
interneuronTTP_HIPP = metricTable.TTP(HIPPLocsIN_Logical);
principalCellTTP_CING = metricTable.TTP(CINGLocsPC_Logical);
interneuronTTP_CING = metricTable.TTP(CINGLocsIN_Logical);
principalCellTTP_pCING = metricTable.TTP(pCINGLocsPC_Logical);
interneuronTTP_pCING = metricTable.TTP(pCINGLocsIN_Logical);

% ONE-wAY ANOVA FOR TTP (removed pCING from this analysis)
y_TTP_PC = [principalCellTTP_AMY', principalCellTTP_OFC', principalCellTTP_HIPP', principalCellTTP_CING'];
group_TTP_PC = repelem(1:4, 1, [numel(principalCellTTP_AMY'),numel(principalCellTTP_OFC'),numel(principalCellTTP_HIPP'), numel(principalCellTTP_CING')]);

    [p_TTP_PC_Region, tbl_TTP_PC_Region, stats_TTP_PC_Region] = anova1(y_TTP_PC,group_TTP_PC)
    if p_TTP_PC_Region <= 0.05
        results_TTP_PC_Region = multcompare(stats_TTP_PC_Region)
    else
    end % OFC PCs are significant longer TTPs than other units.

y_TTP_IN = [interneuronTTP_AMY', interneuronTTP_OFC', interneuronTTP_HIPP', interneuronTTP_CING'];
group_TTP_IN = repelem(1:4, 1, [numel(interneuronTTP_AMY'),numel(interneuronTTP_OFC'),numel(interneuronTTP_HIPP'), numel(interneuronTTP_CING')]);

    [p_TTP_IN_Region, tbl_TTP_IN_Region, stats_TTP_IN_Region] = anova1(y_TTP_IN,group_TTP_IN)
    if p_TTP_IN_Region <= 0.05
        results_TTP_IN_Region = multcompare(stats_TTP_IN_Region)
    else
    end % No differences for INs.

% segment data (FWHM)
principalCellFWHM_AMY = metricTable.FWHM(AMYLocsPC_Logical);
interneuronFWHM_AMY = metricTable.FWHM(AMYLocsIN_Logical);
principalCellFWHM_OFC = metricTable.FWHM(OFCLocsPC_Logical);
interneuronFWHM_OFC = metricTable.FWHM(OFCLocsIN_Logical);
principalCellFWHM_HIPP = metricTable.FWHM(HIPPLocsPC_Logical);
interneuronFWHM_HIPP = metricTable.FWHM(HIPPLocsIN_Logical);
principalCellFWHM_CING = metricTable.FWHM(CINGLocsPC_Logical);
interneuronFWHM_CING = metricTable.FWHM(CINGLocsIN_Logical);
principalCellFWHM_pCING = metricTable.FWHM(pCINGLocsPC_Logical);
interneuronFWHM_pCING = metricTable.FWHM(pCINGLocsIN_Logical);

% ONE-wAY ANOVA FOR FWHM (removed pCING from this analysis)
y_FWHM_PC = [principalCellFWHM_AMY', principalCellFWHM_OFC', principalCellFWHM_HIPP', principalCellFWHM_CING'];
group_FWHM_PC = repelem(1:4, 1, [numel(principalCellFWHM_AMY'),numel(principalCellFWHM_OFC'),numel(principalCellFWHM_HIPP'), numel(principalCellFWHM_CING')]);

    [p_FWHM_PC_Region, tbl_FWHM_PC_Region, stats_FWHM_PC_Region] = anova1(y_FWHM_PC,group_FWHM_PC)
    if p_FWHM_PC_Region <= 0.05
        results_FWHM_PC_Region = multcompare(stats_FWHM_PC_Region)
    else
    end % no differences.

y_FWHM_IN = [interneuronFWHM_AMY', interneuronFWHM_OFC', interneuronFWHM_HIPP', interneuronFWHM_CING'];
group_FWHM_IN = repelem(1:4, 1, [numel(interneuronFWHM_AMY'),numel(interneuronFWHM_OFC'),numel(interneuronFWHM_HIPP'), numel(interneuronFWHM_CING')]);

    [p_FWHM_IN_Region, tbl_FWHM_IN_Region, stats_FWHM_IN_Region] = anova1(y_FWHM_IN,group_FWHM_IN)
    if p_FWHM_IN_Region <= 0.05
        results_FWHM_IN_Region = multcompare(stats_FWHM_IN_Region)
    else
    end % No differences for INs.

% segment data (ASYM)
principalCellASYM_AMY = metricTable.ASYM(AMYLocsPC_Logical);
interneuronASYM_AMY = metricTable.ASYM(AMYLocsIN_Logical);
principalCellASYM_OFC = metricTable.ASYM(OFCLocsPC_Logical);
interneuronASYM_OFC = metricTable.ASYM(OFCLocsIN_Logical);
principalCellASYM_HIPP = metricTable.ASYM(HIPPLocsPC_Logical);
interneuronASYM_HIPP = metricTable.ASYM(HIPPLocsIN_Logical);
principalCellASYM_CING = metricTable.ASYM(CINGLocsPC_Logical);
interneuronASYM_CING = metricTable.ASYM(CINGLocsIN_Logical);
principalCellASYM_pCING = metricTable.ASYM(pCINGLocsPC_Logical);
interneuronASYM_pCING = metricTable.ASYM(pCINGLocsIN_Logical);

% ONE-wAY ANOVA FOR FWHM (removed pCING from this analysis)
y_ASYM_PC = [principalCellASYM_AMY', principalCellASYM_OFC', principalCellASYM_HIPP', principalCellASYM_CING'];
group_ASYM_PC = repelem(1:4, 1, [numel(principalCellASYM_AMY'),numel(principalCellASYM_OFC'),numel(principalCellASYM_HIPP'), numel(principalCellASYM_CING')]);

    [p_ASYM_PC_Region, tbl_ASYM_PC_Region, stats_ASYM_PC_Region] = anova1(y_ASYM_PC,group_ASYM_PC)
    if p_ASYM_PC_Region <= 0.05
        results_ASYM_PC_Region = multcompare(stats_ASYM_PC_Region)
    else
    end % OFC and Cingulate Differences.
    
y_ASYM_IN = [interneuronASYM_AMY', interneuronASYM_OFC', interneuronASYM_HIPP', interneuronASYM_CING'];
group_ASYM_IN = repelem(1:4, 1, [numel(interneuronASYM_AMY'),numel(interneuronASYM_OFC'),numel(interneuronASYM_HIPP'), numel(interneuronASYM_CING')]);

    [p_ASYM_IN_Region, tbl_ASYM_IN_Region, stats_ASYM_IN_Region] = anova1(y_ASYM_IN,group_ASYM_IN)
    if p_ASYM_IN_Region <= 0.05
        results_ASYM_IN_Region = multcompare(stats_ASYM_IN_Region)
    else
    end % No differences for INs.

% FIGURE
figure(555)
% Trough To Peak (region and cell type)
subplot(3,1,1)
boxchartData = [repmat(0, length(principalCellTTP_AMY), 1);  
                repmat(1, length(interneuronTTP_AMY), 1);  
                repmat(2, length(principalCellTTP_OFC), 1); 
                repmat(3, length(interneuronTTP_OFC), 1);
                repmat(4, length(principalCellTTP_HIPP), 1);  
                repmat(5, length(interneuronTTP_HIPP), 1);
                repmat(6, length(principalCellTTP_CING), 1); 
                repmat(7, length(interneuronTTP_CING), 1);
                repmat(8, length(principalCellTTP_pCING), 1);  
                repmat(9, length(interneuronTTP_pCING), 1)];  
boxchartValues = [principalCellTTP_AMY; interneuronTTP_AMY; principalCellTTP_OFC; interneuronTTP_OFC; principalCellTTP_HIPP;...
                  interneuronTTP_HIPP; principalCellTTP_CING; interneuronTTP_CING;  principalCellTTP_pCING; interneuronTTP_pCING];
b = boxchart(boxchartData, boxchartValues,'MarkerStyle','none');
hold on;
jitterAmount = 0.3; 
jitteredPrincipalCellX_AMY = repmat(0, length(principalCellTTP_AMY), 1) + jitterAmount * (rand(length(principalCellTTP_AMY), 1) - 0.5);
jitteredInterneuronX_AMY = repmat(1, length(interneuronTTP_AMY), 1) + jitterAmount * (rand(length(interneuronTTP_AMY), 1) - 0.5);
jitteredPrincipalCellX_OFC = repmat(2, length(principalCellTTP_OFC), 1) + jitterAmount * (rand(length(principalCellTTP_OFC), 1) - 0.5);
jitteredInterneuronX_OFC = repmat(3, length(interneuronTTP_OFC), 1) + jitterAmount * (rand(length(interneuronTTP_OFC), 1) - 0.5);
jitteredPrincipalCellX_HIPP = repmat(4, length(principalCellTTP_HIPP), 1) + jitterAmount * (rand(length(principalCellTTP_HIPP), 1) - 0.5);
jitteredInterneuronX_HIPP = repmat(5, length(interneuronTTP_HIPP), 1) + jitterAmount * (rand(length(interneuronTTP_HIPP), 1) - 0.5);
jitteredPrincipalCellX_CING = repmat(6, length(principalCellTTP_CING), 1) + jitterAmount * (rand(length(principalCellTTP_CING), 1) - 0.5);
jitteredInterneuronX_CING = repmat(7, length(interneuronTTP_CING), 1) + jitterAmount * (rand(length(interneuronTTP_CING), 1) - 0.5);
jitteredPrincipalCellX_pCING = repmat(8, length(principalCellTTP_pCING), 1) + jitterAmount * (rand(length(principalCellTTP_pCING), 1) - 0.5);
jitteredInterneuronX_pCING = repmat(9, length(interneuronTTP_pCING), 1) + jitterAmount * (rand(length(interneuronTTP_pCING), 1) - 0.5);
% Create the scatter plot
scatter(jitteredPrincipalCellX_AMY, principalCellTTP_AMY, 40, RegionColors(1,:), 'filled', 'Marker','^');
scatter(jitteredInterneuronX_AMY, interneuronTTP_AMY, 40,RegionColors(1,:), 'filled');
scatter(jitteredPrincipalCellX_OFC, principalCellTTP_OFC, 40, RegionColors(2,:), 'filled', 'Marker','^');
scatter(jitteredInterneuronX_OFC, interneuronTTP_OFC, 40, RegionColors(2,:), 'filled');
scatter(jitteredPrincipalCellX_HIPP, principalCellTTP_HIPP, 40, RegionColors(3,:), 'filled','Marker','^');
scatter(jitteredInterneuronX_HIPP, interneuronTTP_HIPP, 40, RegionColors(3,:), 'filled');
scatter(jitteredPrincipalCellX_CING, principalCellTTP_CING, 40, RegionColors(4,:), 'filled', 'Marker','^');
scatter(jitteredInterneuronX_CING, interneuronTTP_CING, 40, RegionColors(4,:), 'filled');
scatter(jitteredPrincipalCellX_pCING, principalCellTTP_pCING, 40, RegionColors(5,:), 'filled', 'Marker','^');
scatter(jitteredInterneuronX_pCING, interneuronTTP_pCING, 40, RegionColors(5,:), 'filled');
% Customize the plot
xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
xticklabels({'AMY Principal Cells', 'AMY Interneurons', 'OFC Principal Cells', 'OFC Interneurons','HIPP Principal Cells',...
             'HIPP Interneurons', 'CING Principal Cells', 'CING Interneurons', 'pCING Principal Cells', 'pCING Interneurons'});
xlabel('Cell Type');
ylabel('TTP (ms)');
title('TTP by Cell Type and Region');
hold off;
% FWHM (region and cell type)
subplot(3,1,2)
% Combine the data for the boxchart (4 groups: AMY principal cells, AMY interneurons, OFC principal cells, OFC interneurons)
boxchartData = [repmat(0, length(principalCellFWHM_AMY), 1);  
                repmat(1, length(interneuronFWHM_AMY), 1);  
                repmat(2, length(principalCellFWHM_OFC), 1); 
                repmat(3, length(interneuronFWHM_OFC), 1);
                repmat(4, length(principalCellFWHM_HIPP), 1);  
                repmat(5, length(interneuronFWHM_HIPP), 1);
                repmat(6, length(principalCellFWHM_CING), 1); 
                repmat(7, length(interneuronFWHM_CING), 1);
                repmat(8, length(principalCellFWHM_pCING), 1);  
                repmat(9, length(interneuronFWHM_pCING), 1)];  
boxchartValues = [principalCellFWHM_AMY; interneuronFWHM_AMY; principalCellFWHM_OFC; interneuronFWHM_OFC; principalCellFWHM_HIPP;...
                  interneuronFWHM_HIPP; principalCellFWHM_CING; interneuronFWHM_CING;  principalCellFWHM_pCING; interneuronFWHM_pCING];
b = boxchart(boxchartData, boxchartValues,'MarkerStyle','none');
hold on;
jitterAmount = 0.3;
% Create the scatter plot with jitter
scatter(jitteredPrincipalCellX_AMY, principalCellFWHM_AMY, 40, RegionColors(1,:), 'filled', 'Marker','^');
scatter(jitteredInterneuronX_AMY, interneuronFWHM_AMY, 40, RegionColors(1,:), 'filled');
scatter(jitteredPrincipalCellX_OFC, principalCellFWHM_OFC, 40, RegionColors(2,:), 'filled', 'Marker','^');
scatter(jitteredInterneuronX_OFC, interneuronFWHM_OFC, 40, RegionColors(2,:), 'filled');
scatter(jitteredPrincipalCellX_HIPP, principalCellFWHM_HIPP, 40, RegionColors(3,:), 'filled', 'Marker','^');
scatter(jitteredInterneuronX_HIPP, interneuronFWHM_HIPP, 40, RegionColors(3,:), 'filled');
scatter(jitteredPrincipalCellX_CING, principalCellFWHM_CING, 40, RegionColors(4,:), 'filled', 'Marker','^');
scatter(jitteredInterneuronX_CING, interneuronFWHM_CING, 40, RegionColors(4,:), 'filled');
scatter(jitteredPrincipalCellX_pCING, principalCellFWHM_pCING, 40, RegionColors(5,:), 'filled', 'Marker','^');
scatter(jitteredInterneuronX_pCING, interneuronFWHM_pCING, 40, RegionColors(5,:), 'filled');
% Customize the plot
xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
xticklabels({'AMY Principal Cells', 'AMY Interneurons', 'OFC Principal Cells', 'OFC Interneurons','HIPP Principal Cells',...
             'HIPP Interneurons', 'CING Principal Cells', 'CING Interneurons', 'pCING Principal Cells', 'pCING Interneurons'});
xlabel('Cell Type');
ylabel('FWHM (ms)');
title('FWHM by Cell Type and Region');
hold off;
% ASYM (region and cell type)
subplot(3,1,3)
% Combine the data for the boxchart (4 groups: AMY principal cells, AMY interneurons, OFC principal cells, OFC interneurons)
boxchartData = [repmat(0, length(principalCellASYM_AMY), 1);  
                repmat(1, length(interneuronASYM_AMY), 1);  
                repmat(2, length(principalCellASYM_OFC), 1); 
                repmat(3, length(interneuronASYM_OFC), 1);
                repmat(4, length(principalCellASYM_HIPP), 1);  
                repmat(5, length(interneuronASYM_HIPP), 1);
                repmat(6, length(principalCellASYM_CING), 1); 
                repmat(7, length(interneuronASYM_CING), 1);
                repmat(8, length(principalCellASYM_pCING), 1);  
                repmat(9, length(interneuronASYM_pCING), 1)];  
boxchartValues = [principalCellASYM_AMY; interneuronASYM_AMY; principalCellASYM_OFC; interneuronASYM_OFC; principalCellASYM_HIPP;...
                  interneuronASYM_HIPP; principalCellASYM_CING; interneuronASYM_CING;  principalCellASYM_pCING; interneuronASYM_pCING];
b = boxchart(boxchartData, boxchartValues,'MarkerStyle','none');
hold on;
jitterAmount = 0.3;
% Create the scatter plot with jitter
scatter(jitteredPrincipalCellX_AMY, principalCellASYM_AMY, 40, RegionColors(1,:), 'filled', 'Marker','^');
scatter(jitteredInterneuronX_AMY, interneuronASYM_AMY, 40, RegionColors(1,:), 'filled');
scatter(jitteredPrincipalCellX_OFC, principalCellASYM_OFC, 40, RegionColors(2,:),  'filled', 'Marker','^');
scatter(jitteredInterneuronX_OFC, interneuronASYM_OFC, 40, RegionColors(2,:), 'filled');
scatter(jitteredPrincipalCellX_HIPP, principalCellASYM_HIPP, 40, RegionColors(3,:),'filled', 'Marker','^');
scatter(jitteredInterneuronX_HIPP, interneuronASYM_HIPP, 40, RegionColors(3,:),'filled');
scatter(jitteredPrincipalCellX_CING, principalCellASYM_CING, 40, RegionColors(4,:), 'filled', 'Marker','^');
scatter(jitteredInterneuronX_CING, interneuronASYM_CING, 40,  RegionColors(4,:), 'filled');
scatter(jitteredPrincipalCellX_pCING, principalCellASYM_pCING, 40, RegionColors(5,:), 'filled', 'Marker','^');
scatter(jitteredInterneuronX_pCING, interneuronASYM_pCING, 40, RegionColors(5,:), 'filled');
% Customize the plot
xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
xticklabels({'AMY Principal Cells', 'AMY Interneurons', 'OFC Principal Cells', 'OFC Interneurons','HIPP Principal Cells',...
             'HIPP Interneurons', 'CING Principal Cells', 'CING Interneurons', 'pCING Principal Cells', 'pCING Interneurons'});
xlabel('Cell Type');
ylabel('ASYM (ratio)');
title('ASYM by Cell Type and Region');
hold off;

maximize(figure(555))
% save figure
keyboard
saveas(555,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\WaveFormClassification\',sprintf('meanMetrics_CellType_by_Region_AcrossSubs.pdf')))
close(555)


%% suppression latency for each region... because there is a significant difference NOT ANY MORE!

figure(678);
hold on;
jitter_amount = 0.1;
jitter_AMY = jitter_amount * randn(size(SuppressionLat_AMY));
jitter_OFC = jitter_amount * randn(size(SuppressionLat_OFC));
jitter_HIPP = jitter_amount * randn(size(SuppressionLat_HIPP));
jitter_CING = jitter_amount * randn(size(SuppressionLat_CING));
jitter_pCING = jitter_amount * randn(size(SuppressionLat_pCING));
% Boxplot for each region
boxplot(SuppressionLat_AMY, 'colors', Muted_Blue, 'orientation','horizontal', 'Positions',1);
boxplot(SuppressionLat_OFC, 'colors', Muted_Purple, 'orientation','horizontal', 'Positions',2);
boxplot(SuppressionLat_HIPP, 'colors', Dark_Blue, 'orientation','horizontal', 'Positions',3);
boxplot(SuppressionLat_CING, 'colors', Muted_Green, 'orientation','horizontal', 'Positions',4);
boxplot(SuppressionLat_pCING, 'colors', Dark_Brown, 'orientation','horizontal', 'Positions',5);
% Scatter plot for each region
scatter(SuppressionLat_AMY, ones(size(SuppressionLat_AMY)) + jitter_AMY, 70, Muted_Blue, 'filled', '^', 'MarkerEdgeColor', 'k');
scatter(SuppressionLat_OFC, 2*ones(size(SuppressionLat_OFC)) + jitter_OFC, 70, Muted_Purple, 'filled', '^', 'MarkerEdgeColor', 'k');
scatter(SuppressionLat_HIPP, 3*ones(size(SuppressionLat_HIPP)) + jitter_HIPP, 70, Dark_Blue, 'filled', '^', 'MarkerEdgeColor', 'k');
scatter(SuppressionLat_CING, 4*ones(size(SuppressionLat_CING)) + jitter_CING, 70, Muted_Green, 'filled', '^', 'MarkerEdgeColor', 'k');
scatter(SuppressionLat_pCING, 5*ones(size(SuppressionLat_pCING)) + jitter_pCING, 70, Dark_Brown, 'filled', '^', 'MarkerEdgeColor', 'k');
xlabel('Suppression Latency (s)', 'FontSize', 12);
ylabel('Regions', 'FontSize', 12);
title('Suppression Latency by Region', 'FontSize', 14);
ylim([0.5, 5.5]);
xlim([0 2.5])
xticks(0:0.2:1.4)
yticks(1:5);
yticklabels({'AMY', 'OFC', 'HIPP', 'CING', 'pCING'});
axis square;
% add significant markers
%plot([1.3, 2.1], [4, 2], 'k-', 'LineWidth', 2); % Line connecting CING and OFC
%text(1.4, 3.5, '*', 'FontSize', 14, 'HorizontalAlignment', 'center'); % Star above the line indicating significance
hold off;
% sav
saveas(678,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\Unit_Modulation\',sprintf('suppressionLatency_by_Region_AcrossSubs.pdf')))
close(678)

%% latency by TTP....... COM EBACK TO THIS. 

%sigsuppOFCLocs need to use this?

% scatter (all: TTP vs latency)
% mdl = fitlm(metricTable.TTP(OFCLocsPC_Logical), [SPESstruct(OFCLocsPC_Logical).suppLatencyAll]) % NS
% plot(mdl) 
% 
% mdl = fitlm(metricTable.TTP(OFCLocsPC_Logical), [SPESstruct(OFCLocsPC_Logical).Minimum_SuppressionAll]) % SIGNIFICANT
% plot(mdl)
% 
% mdl = fitlm(metricTable.TTP(OFCLocsPC_Logical), [SPESstruct(OFCLocsPC_Logical).suppressionDurationTimeAllz]) % 
% plot(mdl)

%% suppression latency by region and cell type

% suppression latency region and cell type
SuppressionLat_AMY_PC = ([sigSPESstruct(sigAMYLocs & suppressionSigUnits & ~isolatedClusterSig').suppLatencyAll]);
SuppressionLat_AMY_IN = ([sigSPESstruct(sigAMYLocs & suppressionSigUnits & isolatedClusterSig').suppLatencyAll]);
SuppressionLat_OFC_PC = ([sigSPESstruct(sigOFCLocs & suppressionSigUnits & ~isolatedClusterSig').suppLatencyAll]);
SuppressionLat_OFC_IN = ([sigSPESstruct(sigOFCLocs & suppressionSigUnits & isolatedClusterSig').suppLatencyAll]);
SuppressionLat_HIPP_PC = ([sigSPESstruct(sigHIPPLocs & suppressionSigUnits & ~isolatedClusterSig').suppLatencyAll]);
SuppressionLat_HIPP_IN = ([sigSPESstruct(sigHIPPLocs & suppressionSigUnits & isolatedClusterSig').suppLatencyAll]);
SuppressionLat_CING_PC = ([sigSPESstruct(sigCINGLocs & suppressionSigUnits & ~isolatedClusterSig').suppLatencyAll]);
SuppressionLat_CING_IN = ([sigSPESstruct(sigCINGLocs & suppressionSigUnits & isolatedClusterSig').suppLatencyAll]);
SuppressionLat_pCING_PC = ([sigSPESstruct(sigpCINGLocs & suppressionSigUnits & ~isolatedClusterSig').suppLatencyAll]);
SuppressionLat_pCING_IN = ([sigSPESstruct(sigpCINGLocs & suppressionSigUnits & isolatedClusterSig').suppLatencyAll]);

% ONE-WAY ANOVA FOR SUPPRESSION LATENCY (removed pCING from this analysis)
y_suppLatCellType_Region = [SuppressionLat_AMY_PC, SuppressionLat_AMY_IN, SuppressionLat_OFC_PC, SuppressionLat_OFC_IN,SuppressionLat_HIPP_PC, SuppressionLat_HIPP_IN, SuppressionLat_CING_PC, SuppressionLat_CING_IN];
group_suppLatCellType_Region = repelem(1:8, 1, [numel(SuppressionLat_AMY_PC'),numel(SuppressionLat_AMY_IN'),numel(SuppressionLat_OFC_PC'), numel(SuppressionLat_OFC_IN'),numel(SuppressionLat_HIPP_PC'),numel(SuppressionLat_HIPP_IN'),numel(SuppressionLat_CING_PC'),numel(SuppressionLat_CING_IN')]);

[p_suppLatCT_Region, tbl_suppLatCT_Region, stats_suppLatCT_Region] = anova1(y_suppLatCellType_Region,group_suppLatCellType_Region)
if p_suppLatCT_Region <= 0.05
    results_suppLatCT_Region = multcompare(stats_suppLatCT_Region)
else
end % NS.

% plot: supp latency by region and cell type
figure(679);
hold on;
jitter_amount = 0.1;
jitter_AMY_PC = jitter_amount * randn(size(SuppressionLat_AMY_PC));
jitter_AMY_IN = jitter_amount * randn(size(SuppressionLat_AMY_IN));
jitter_OFC_PC = jitter_amount * randn(size(SuppressionLat_OFC_PC));
jitter_OFC_IN = jitter_amount * randn(size(SuppressionLat_OFC_IN));
jitter_HIPP_PC = jitter_amount * randn(size(SuppressionLat_HIPP_PC));
jitter_HIPP_IN = jitter_amount * randn(size(SuppressionLat_HIPP_IN));
jitter_CING_PC = jitter_amount * randn(size(SuppressionLat_CING_PC));
jitter_CING_IN = jitter_amount * randn(size(SuppressionLat_CING_IN));
jitter_pCING_PC = jitter_amount * randn(size(SuppressionLat_pCING_PC));
jitter_pCING_IN = jitter_amount * randn(size(SuppressionLat_pCING_IN));

% Boxplot for each region
boxplot(SuppressionLat_AMY_PC, 'colors', Muted_Blue, 'orientation','horizontal', 'Positions',1);
boxplot(SuppressionLat_AMY_IN, 'colors', Muted_Blue, 'orientation','horizontal', 'Positions',2);
boxplot(SuppressionLat_OFC_PC, 'colors', Muted_Purple, 'orientation','horizontal', 'Positions',3);
boxplot(SuppressionLat_OFC_IN, 'colors', Muted_Purple, 'orientation','horizontal', 'Positions',4);
boxplot(SuppressionLat_HIPP_PC, 'colors', Dark_Blue, 'orientation','horizontal', 'Positions',5);
boxplot(SuppressionLat_HIPP_IN, 'colors', Dark_Blue, 'orientation','horizontal', 'Positions',6);
boxplot(SuppressionLat_CING_PC, 'colors', Muted_Green, 'orientation','horizontal', 'Positions',7);
boxplot(SuppressionLat_CING_IN, 'colors', Muted_Green, 'orientation','horizontal', 'Positions',8);
boxplot(SuppressionLat_pCING_PC, 'colors', Dark_Brown, 'orientation','horizontal', 'Positions',9);
boxplot(SuppressionLat_pCING_IN, 'colors', Dark_Brown, 'orientation','horizontal', 'Positions',10);

% Scatter plot for each region
scatter(SuppressionLat_AMY_PC, ones(size(SuppressionLat_AMY_PC)) + jitter_AMY_PC, 100, Muted_Blue, 'filled', '^', 'MarkerEdgeColor', 'k');
scatter(SuppressionLat_AMY_IN, 2*ones(size(SuppressionLat_AMY_IN)) + jitter_AMY_IN, 100, Muted_Blue, 'filled', 'o', 'MarkerEdgeColor', 'k');
scatter(SuppressionLat_OFC_PC, 3*ones(size(SuppressionLat_OFC_PC)) + jitter_OFC_PC, 100, Muted_Purple, 'filled', '^', 'MarkerEdgeColor', 'k');
scatter(SuppressionLat_OFC_IN, 4*ones(size(SuppressionLat_OFC_IN)) + jitter_OFC_IN, 100, Muted_Purple, 'filled', 'o', 'MarkerEdgeColor', 'k');
scatter(SuppressionLat_HIPP_PC, 5*ones(size(SuppressionLat_HIPP_PC)) + jitter_HIPP_PC, 100, Dark_Blue, 'filled', '^', 'MarkerEdgeColor', 'k');
scatter(SuppressionLat_HIPP_IN, 6*ones(size(SuppressionLat_HIPP_IN)) + jitter_HIPP_IN, 100, Dark_Blue, 'filled', 'o', 'MarkerEdgeColor', 'k');
scatter(SuppressionLat_CING_PC, 7*ones(size(SuppressionLat_CING_PC)) + jitter_CING_PC, 100, Muted_Green, 'filled', '^', 'MarkerEdgeColor', 'k');
scatter(SuppressionLat_CING_IN, 8*ones(size(SuppressionLat_CING_IN)) + jitter_CING_IN, 100, Muted_Green, 'filled', 'o', 'MarkerEdgeColor', 'k');
scatter(SuppressionLat_pCING_PC, 9*ones(size(SuppressionLat_pCING_PC)) + jitter_pCING_PC, 100, Dark_Brown, 'filled', '^', 'MarkerEdgeColor', 'k');
scatter(SuppressionLat_pCING_IN, 10*ones(size(SuppressionLat_pCING_IN)) + jitter_pCING_IN, 100, Dark_Brown, 'filled', 'o', 'MarkerEdgeColor', 'k');

xlabel('Suppression Latency (s)', 'FontSize', 12);
ylabel('Regions and Cell Type', 'FontSize', 12);
title('Suppression Latency by Region and Cell Type', 'FontSize', 14);
ylim([0.5, 10.5]);
xlim([0 2.6])
xticks(0:0.2:2.6)
yticks(1:10);
yticklabels({'AMY PC', 'AMY IN', 'OFC PC', 'OFC IN', 'HIPP PC', 'HIPP IN', 'CING PC', 'CING IN', 'pCING PC', 'pCING IN'});
axis square;
% add significant markers
%plot([1.3, 2.1], [4, 2], 'k-', 'LineWidth', 2); % Line connecting CING and OFC
%text(1.4, 3.5, '*', 'FontSize', 14, 'HorizontalAlignment', 'center'); % Star above the line indicating significance
hold off;
% saving
maximize(679)
saveas(679,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\Unit_Modulation\',sprintf('suppressionLatency_by_RegionCellType_AcrossSubs.pdf')))
close(679)

saveas(679,fullfile('\\155.100.91.44\D\Data\Rhiannon\CCEPS\Figures\AcrossSubjects\Unit_Modulation\',sprintf('test.pdf')))


%% save wide Probs and FSprobs to SPESstruct for further analysis:

save WFprobs.mat FSprobs wideProbs metrics KMeanTwo KMeanThree isolatedCluster;

 %% including Ed's autocorrelation function so as to avoid using his
% beautiful neuron classes.
function [ac,lags] = autocorr(times,bins)

t_subset = [-Inf Inf];
sub_times = times;

d = diff(bins)/2;
edges = [bins(1)-d(1), bins(1:end-1)+d, bins(end)+d(end)];

bigMat = repmat(sub_times,1,length(sub_times));
bigMat = bigMat - diag(bigMat)';
bigMat = bigMat * 1e3; % convert to milliseconds

vals = histcounts(bigMat(:),edges);
[~,wh] = min(abs(edges));
vals(wh) = vals(wh) - length(sub_times); % remove self from AC
ac_data.xc = vals;
ac_data.lags = bins;
ac_data.time_subset = t_subset;
ac = ac_data.xc;
lags = ac_data.lags;

end

% how many macros did each patient have? 
%macro_count = [122 97 98  90 123 118 162 110 184 130 174 111 151 119 133 181 143 155 129 170 135 201 170 150 159 253 145 232];
%sum(macro_count)
%std(macro_count)
%mean(macro_count)
