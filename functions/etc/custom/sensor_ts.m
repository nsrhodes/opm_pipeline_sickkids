function sensor_ts(data, cfg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = sensor_ts(cfg, rawData)
%
% Function to calculate and visualize sensor-level timeseries for MEG data
%
% Inputs:
%   data                (structure) MEG data in FieldTrip format
%   cfg                 (structure) Configuration structure containing:
%   (optional)
%   	foi             (array)     Range of frequencies to compute the PSD
%                                   Default: [1, 150]
%                                   Default: true
%       trial_length	(float)     Length of segments to compute PSD over
%                                   Default: 1s
%       average         (string)    'mean' or 'median': how to average over
%                                   trials/channels
%                                   Default: median
%
% Outputs:
%   gcf                 (figure)    Figure showing the timeseries
%
% Created by: Marlee M Vandewouw
% Created on: January 20th, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DEFAULTS

% Set default trial length
if ~isfield(cfg, 'trial_length')
    cfg.trial_length = 1;
end

% Set default frequency of interest
if ~isfield(cfg, 'foi')
    cfg.foi = [1 150];
end

% Set default averaging method
if ~isfield(cfg, 'average')
    cfg.average = 'median';
end

%% MAIN

% if cfg.trial_length is an empty array we'll get it to guess the
% length of the trial
if isempty(cfg.trial_length)
    cfg.trial_length = range(data.time{1});
end

% Catch if trial length exceeds the data
if cfg.trial_length > range(data.time{1})
    error('You are specifying PSD windows longer than the trial length!')
end

% Bandpass filter
cfg_bp = [];
cfg_bp.demean = 'yes';
cfg_bp.bpfilter = 'yes';
cfg_bp.bpfreq = cfg.foi;
[~, data] = evalc('ft_preprocessing(cfg_bp, data)');

% Split data into epochs
nsamps = (cfg.trial_length)*data.fsample;
ntrials = numel(data.trial);
for ii = 1:ntrials
    beg = 1:nsamps:size(data.trial{ii},2);
    endsamp =  beg+(nsamps-1);
    inRange = ~(beg>=size(data.trial{ii},2)|endsamp>=size(data.trial{ii},2));
    chans = 1:1:(size(data.trial{ii},1));
    if ii == 1
        eD = zeros(length(chans),nsamps,sum(inRange),ntrials);
    end
    for jj = 1:length(inRange)
        if(inRange(jj))
            eD(:,:,jj,ii)=data.trial{ii}(chans,beg(jj):endsamp(jj),:);
        end
    end
end

% Reshape eD so its 3D
eD = reshape(eD,size(eD,1),size(eD,2),[]);

% Calculate median or mean out of all epochs
if strcmp(cfg.average,'median')
    ts = median(eD(:,:,:),3);
else
    ts = mean(eD(:,:,:),3);
end

%% PLOT

% Set the colormap
cmap = turbo(size(ts,1));

% Initialize
fig = figure('units','normalized','outerposition',[0 0 1 1]);
fig.Name = 'Timeseries';
hold on;

% Plot each sensor
tiledlayout(ceil(size(ts, 1)/10), 10);
for ii = 1:size(ts, 1)
    nexttile;
    plot(ts(ii, :), 'Color', cmap(ii, :));
    ylim([-2e-11, 2e-11]);
    xlim tight
    title(data.label{ii});    
end




