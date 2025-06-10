function [metrics,outputWaves] = waveformMetrics(meanWaveforms,varargin)
% Calculate the waveshape-based metrics for mean waveforms of single units.
% Usage:
%  [metrics, outputWaves] = waveformMetrics(meanWaveforms [,settings]);
% where meanWaveforms is an [n x m] matrix of n waveforms sampled across m
% data points
%
% Optional inputs:
%   Fs:             sampling rate [default: 30,000]
%   keyPoint:       data point where the trough should be [default: 200]
%   window:         number of data points to keep around trough [default:
%                   -18:45 (-0.6 to 1.5 ms at 30 kHz)]
%   alignWindow:    data points within which to find the updated trough
%                   [default: 10]
%   smoothFactor:   number of data points over which to smooth the mean
%                   waveform [default: 3 (3 at 30 kHz gives 0.1 ms smooth)]
%   idealizedOrder: polynomial order for fitting the "idealized"
%                   repolarization phase [default: 4]
%   oversampling:   multipe to up-sample at [default: ceil(1e5/settings.Fs)
%                   (i.e. go to next highest value above a new sample rate
%                   of 100 kHz)]
%   tolerance:      tolerance for near-zero values in the derivative when
%                   finding the peak in the trough-to-peak value [default:
%                   0 (i.e. isn't used by default)]
%
% Output:
% 1: a metrics struct containing:
%   troughToPeak:       the number of milliseconds from trough until
%                       following peak in the upsampled waveform
%   troughToPeakIdeal:  the number of milliseconds from trough until
%                       following peak in the idealized polynomial fit
%   FWHM:               the number of milliseconds for the "full-width at
%                       half maximum" (i.e. the spike half-width)
%   asymmetry:          the asymmetry between the maximal values pre- and
%                       post-spike (ratio)
% 2: a matrix of the normalized, upsampled waveforms

% Settings: (see above for explanation of each, update them at runtime as
% name, value pairs in input)
settings.Fs = 3e4;
settings.keyPoint = 14; % OG: 200 % WASHU: was 14
settings.window = -11:34; %-10:19; % WASHU: changed from -11:34 for washU units (have less samples)
settings.alignWindow = 10; % WASHU: changed from 10.
settings.smoothFactor = 3;
settings.idealizedOrder = 4;
settings.oversampling = ceil(1e5/settings.Fs);
settings.tolerance = 0;
unitNumber = 0;
unclassifiedUnits = 0;
% Go through and update above settings if user has given any:
allowable = fieldnames(settings);

if mod(length(varargin),2) ~= 0
    error('Settings inputs must be in name, value pairs');
end
for v = 1:2:length(varargin)
    if find(ismember(allowable,varargin{v}))
        settings.(varargin{v}) = varargin{v+1};
    else
        disp([9 'Not assigning ''' varargin{v} ''': not a valid setting in cellTypeSubclassify']);
    end
end

% Calculate timing vectors for original and upsampled data:
timeWaveformIn = settings.window/(settings.Fs/1e3);
timeWaveform = interp1(timeWaveformIn,timeWaveformIn,timeWaveformIn(1):mean(diff(timeWaveformIn))/settings.oversampling:timeWaveformIn(end),'spline');
troughInterval = [find(timeWaveform >= -0.25,1) find(timeWaveform >= 0.25,1)];

% Remove mean offset:
meanWaveforms = meanWaveforms - nanmean(meanWaveforms,2);

metrics.troughToPeak = NaN(1,size(meanWaveforms,1));
metrics.troughToPeakIdeal = NaN(1,size(meanWaveforms,1));
metrics.FWHM = NaN(1,size(meanWaveforms,1));
metrics.asymmetry = NaN(1,size(meanWaveforms,1));
outputWaves = NaN(size(meanWaveforms,1),length(timeWaveform));

% Loop through each waveform and store relevant metrics:
for i = size(meanWaveforms,1):-1:1
    % Smooth the mean waveform, and keep subset around the new trough:
    smoothWave = smooth(meanWaveforms(i,:),settings.smoothFactor);
    [~,w] = min(smoothWave(settings.keyPoint+(-settings.alignWindow:settings.alignWindow)));
    %waveidx = w-1+settings.keyPoint-settings.alignWindow+settings.window; % added in
    %if max(waveidx) > 48; waveidx = waveidx-1; end % added in
    try
        %         w_Check = w-1 + settings.keyPoint - settings.alignWindow + settings.window; % checking to make sure there are no zeros in the data.
        %         % Check if the first element is zero
        %         if w_Check(1) == 0
        %             % Add 1 to all elements in the array
        %             w= w_Check + 1;
        %         else
        %         end

        keptWave = smoothWave(w-1+settings.keyPoint-settings.alignWindow+settings.window); %(waveidx) % doesnt work if there is a zero in calculated.

        % Upsample, z-score and re-find the trough in the smoothed waveform:
        interpWave = interp1(timeWaveformIn,keptWave,timeWaveform,'spline');
        zs = zscore(interpWave);
        [~,trgh] = min(zs(troughInterval(1):troughInterval(2)));
        [maxPost,pk] = max(interpWave(trgh+troughInterval(1):end));
        maxPre = max(interpWave(1:trgh+troughInterval(1)));
        wvTrgh = trgh+troughInterval(1)-1;
        % Calculate the idealized repolarization phase with a polynomial fit:
        t = (0:(length(interpWave)-wvTrgh));
        p = polyfit(t,interpWave(wvTrgh:end),settings.idealizedOrder);
        returnWvIdeal = polyval(p,t);
        [~,retTimeIdeal] = max(returnWvIdeal);
        % Find where the regular waveform hits local maxima (within tolerance):
        dv = diff(zs(wvTrgh:end));
        [~,inds] = find(dv < settings.tolerance & zs(wvTrgh+1:end) > -settings.tolerance);
        if isempty(inds) % It didn't start returning, so take where the peak was instead
            retTime = pk;
        else
            retTime = inds(1);
        end
        % Calculate the half-max point for FWHM:
        normWave = interpWave - interpWave(wvTrgh);
        normWave = normWave/max(normWave(wvTrgh:end)) - 0.5;
        inds = find(normWave >= 0);
        pre_inds = inds(inds < wvTrgh);
        post_inds = inds(inds > wvTrgh);
        if ~isempty(pre_inds) && ~isempty(post_inds)
            pre_ind = pre_inds(end);
            post_ind = post_inds(1);

            n = normWave(pre_ind);
            m = normWave(pre_ind+1);
            addition = n/(n-m); % go to "true" zero rather than the value at the data point closest to zero

            n = normWave(post_ind);
            m = normWave(post_ind-1);
            subtraction = n/(n-m); % go to "true" zero rather than the value at the data point closest to zero

            hw = (post_ind - subtraction) - (pre_ind + addition);
        end
        % Store the waveform-based metrics:
        metrics.troughToPeak(i) = retTime/((settings.Fs/1e3)*settings.oversampling);
        metrics.troughToPeakIdeal(i) = retTimeIdeal/((settings.Fs/1e3)*settings.oversampling);
        metrics.asymmetry(i) = (maxPost-maxPre)/(maxPost+maxPre);
        metrics.FWHM(i) = hw/((settings.Fs/1e3)*settings.oversampling);
        outputWaves(i,:) = interpWave;
        metrics.unitNumber(i) = i;


        % end of try statement for unclassifiable units
    catch
        metrics.unclassifiedUnits = unclassifiedUnits+1;
        % metrics.idxTwo = idxTwo;
        % metrics.CTwo = CTwo;
        % metrics.idxThree = idxThree;
        % metrics.CThree = CThree;
    end
end
