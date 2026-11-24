function [data_filt, peakWidths] = opm_filter_noise_peaks(data, noise, save_name, cfg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [data_filt, peakWidths] = opm_filter_noise_peaks(data, noise, cfg)
%
% Function to calculate the PSD on OPM data and corresponding noise, and
% identify peaks which represent frequencies present in both types data
% and that should be bandstop filtered. The width of the peak is determined
% to be the smallest possible for the peak to no longer be identifiable at 
% half the prominence.
%
% Inputs:
%   data                (structure) MEG data in FieldTrip format
%   noise               (structure) MEG noise data in FieldTrip format
%   save_name           (string)    Save name for pre/post filtering figure
%   cfg                 (structure) Configuration structure containing:
%   (optional)
%   	foi             (array)     Range of frequencies to compute the PSD
%                                   Default: [1, 150]
%   	channel         (logical)   Whether to plot the individual 
%                                   channels, or just the the average
%                                   Default: true
%       trial_length	(float)     Length of segments to compute PSD over
%                                   Default: whole segment
%       average         (string)    'mean' or 'median': how to average over
%                                   trials/channels
%                                   Default: median
%   	minPeakDistance (float)     Minimum distance between peaks
%                                   Default: 3
%       minPeakProm     (float)     Minimum peak prominence
%                                   Default: 3e-15
%
% Outputs:
%   data_filt           (structure) Filtered MEG data in FieldTrip format
%   peakWidths          (array)     Widths of the bandstop filter that's
%                                   necessary to remove noise peaks from 
%                                   the data
%
% Created by: Marlee M Vandewouw
% Created on: January 20th, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DEFAULTS

% Set default trial length
if ~isfield(cfg, 'trial_length')
    cfg.trial_length = 1;
end

% Set default to plot all channels
if ~isfield(cfg, 'channel')
    cfg.channel = true;
end

% Set default frequency of interest
if ~isfield(cfg, 'foi')
    cfg.foi = [1 150];
end

% Set default averaging method
if ~isfield(cfg, 'average')
    cfg.average = 'median';
end

% Set default for minimum peak prominance
if ~isfield(cfg, 'minPeakProm')
    cfg.minPeakProm = 3e-15;
end

% Set the default for minimum peak distanece
if ~isfield(cfg, 'minPeakDistance')
    cfg.minPeakDistance = 3;
end

% Catch if trial length exceeds the data
if cfg.trial_length > range(data.time{1})
    error('You are specifying PSD windows longer than the trial length!')
end

%% GET NOISE PEAKS

% Compute the PSD
cfg_noise = cfg;
cfg_noise.trial_length = min([cfg_noise.trial_length, range(noise.time{1})]);
[pow_noise, freq_noise] = custom_psd(noise, cfg_noise);
    
% Calculate median or mean out of all epochs
if strcmp(cfg.average,'median')
    po_noise = median(pow_noise(:,:,:),3);
    po_chavg_noise = median(po_noise, 2);
else
    po_noise = mean(pow_noise(:,:,:),3);
    po_chavg_noise = mean(po_noise, 2);
end

% Get peaks
[~, peaks_noise] = findpeaks(po_chavg_noise, freq_noise, 'MinPeakDistance', cfg.minPeakDistance, 'MinPeakProminence', cfg.minPeakProm);
if any(peaks_noise == 1)
    peaks_noise(peaks_noise == 1) = [];
end

% Get rid of peaks outside of FOI
peaks_noise(peaks_noise < cfg.foi(1)) = [];
peaks_noise(peaks_noise > cfg.foi(2)) = [];

%% IDENTIFY PEAKS IN RAW DATA

% Get peaks of data from the last iteration
[pow, freq] = custom_psd(data, cfg);

% Calculate median or mean out of all epochs
if strcmp(cfg.average,'median')
    po = median(pow(:,:,:),3);
    po_chavg = median(po, 2);
else
    po = mean(pow(:,:,:),3);
    po_chavg = mean(po, 2);
end
po_chavg_store = po_chavg;

% Get peaks
[~, peaks, ~] = findpeaks(po_chavg, freq, 'MinPeakDistance', cfg.minPeakDistance, 'MinPeakProminence', cfg.minPeakProm);

% Get rid of peaks outside of FOI
peaks(peaks < cfg.foi(1)) = [];
peaks(peaks > cfg.foi(2)) = [];

% Get overlap with noise
peakFreqs = intersect(peaks, peaks_noise);

%% DETERMINE PEAK WIDTH
% Start with a small width for each peak, and iteratively get bigger until
% the peak is undetectable at half the prominence

% Initialize
peakWidths = zeros(length(peakFreqs), 2);

% Iterate over the peaks
for ff = 1:length(peakFreqs)
    
    fprintf('   Peak identified at %d - finding width...', peakFreqs(ff));

    % Start with a narrow, 0.1Hz width peak
    width = 0.1;
    removed = false;
    while ~removed 
        
        % Create the width
        peakWidths(ff, :) = [peakFreqs(ff) - width, peakFreqs(ff) + width];
        
        % Bandstop filter all existing widths
        cfg_preproc = [];
        cfg_preproc.bsfilter = 'yes';
        cfg_preproc.bsfreq = peakWidths(ff, :);

        % Now, apply
        try
            [~, data_filt] = evalc('ft_preprocessing(cfg_preproc, data)');
        catch
            cfg_preproc.bsfiltord = 3;
            [~, data_filt] = evalc('ft_preprocessing(cfg_preproc, data)');
        end
        
        % Recalculat the peaks
        [pow, freq] = custom_psd(data_filt, cfg);
        if strcmp(cfg.average,'median')
            po = median(pow(:,:,:),3);
            po_chavg = median(po, 2);
        else
            po = mean(pow(:,:,:),3);
            po_chavg = mean(po, 2);
        end
        [~, peaks, ~]   = findpeaks(po_chavg, freq, 'MinPeakProminence', cfg.minPeakProm/2);
        peaks(peaks < cfg.foi(1)) = [];
        peaks(peaks > cfg.foi(2)) = [];
        
        % Check if peak is still there        
        removed = ~any(peaks == peakFreqs(ff));
        
        % Update the width
        width = width + 0.05;
        
    end
    
    fprintf('done.\n');

end

% Now, do the combined filtering
cfg_preproc = [];
cfg_preproc.bsfilter = 'yes';
cfg_preproc.bsfreq = peakWidths;
try
    [~, data_filt] = evalc('ft_preprocessing(cfg_preproc, data)');
catch
    cfg_preproc.bsfiltord = 3;
    [~, data_filt] = evalc('ft_preprocessing(cfg_preproc, data)');
end    

%% CALCULATE PSD POST-FILTERING

% Get PSD
[pow_post, ~] = custom_psd(data_filt, cfg);
    
% Calculate median or mean out of all epochs
if strcmp(cfg.average,'median')
    post_po = median(pow_post(:,:,:),3);
    post_po_chavg = median(post_po, 2);
else
    post_po = mean(pow_post(:,:,:),3);
    post_po_chavg = mean(post_po, 2);
end

%% PRE/POST PLOT

% Set the colormap
cmap = turbo(size(pow,2));

% Initialize
fig = figure('units','normalized','outerposition',[0 0 1 1]);
fig.Name = 'PSD';
hold on;

% Plot channels if necessary
if cfg.channel
    
    % Plot
    h = plot(freq, post_po', 'LineWidth', 0.5);
    
    % Change the color, transparency and add data tips
    for ii = 1:size(pow, 2)
        h(ii).Color(1:3) = cmap(ii, :);
        h(ii).Color(4) = 0.25;
        h(ii).DataTipTemplate.DataTipRows(end+1) = data.label{ii};
    end
    
end
    
% Plot the average 
p1 = plot(freq, po_chavg_store, '-k', 'LineWidth', 2);
p2 = plot(freq_noise, po_chavg_noise, '-r', 'LineWidth', 2);
p3 = plot(freq, post_po_chavg, '-b', 'LineWidth', 2);

% Plot the overlap of the peaks
set(gca, 'YScale', 'log');
xlim(cfg.foi);
for ii = 1:length(peakFreqs)
    plot([peakFreqs(ii), peakFreqs(ii)], ylim, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
end

% Axis
set(gca, 'YScale', 'log');
xlim(cfg.foi);
xlabel('Frequency (Hz)');
labY = ['$$PSD (' 'fT' ' \sqrt[-1]{Hz}$$)'];
ylabel(labY,'interpreter','latex')
legend([p1, p2, p3], {'Original', 'Noise', 'Filtered'});

% Save
saveas(gcf, save_name);
    
