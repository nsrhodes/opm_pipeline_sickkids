function [pow, freq] = custom_psd(data, cfg, fig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [pow, freq] = custom_psd(cfg, data)
%
% Function to calculate and visualize PSD on MEG data
%
% Inputs:
%   data                (structure) MEG data in FieldTrip format
%   cfg                 (structure) Configuration structure containing:
%   (optional)
%   	foi             (array)     Range of frequencies to compute the PSD
%                                   Default: [1, 150]
%   	channel         (logical)   Whether to plot the individual 
%                                   channels, or just the the average
%                                   Default: true
%       trial_length	(float)     Length of segments to compute PSD over
%                                   Default: 1s
%       average         (string)    'mean' or 'median': how to average over
%                                   trials/channels
%                                   Default: median
%       plot            (logical)   Whether or not to plot
%                                   Default: false
%
% Outputs:
%   pow                 (array)     Power spectral density
%   freq                (array)     Frequency vector
%
%   (optional)
%   gcf                 (figure)    Figure showing the PSD
%
% Adapted from ft_opm_psd (Tim Tierney)
%
% Created by: Marlee M Vandewouw
% Created on: January 20th, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DEFAULTS

% Set default trial length
if isempty(cfg.trial_length)
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

% Set the default to plot
if ~isfield(cfg, 'plot')
    cfg.plot = false;
end

% Catch if trial length exceeds the data
if cfg.trial_length > range(data.time{1})
    error('You are specifying PSD windows longer than the trial length!')
end

%% MAIN

% Catch if trial length exceeds the data
if cfg.trial_length > range(data.time{1})
    error('You are specifying PSD windows longer than the trial length!')
end

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

% Get variables for calculating PSD
fs = data.fsample;
N = size(eD,2);
Nf = ceil((N+1)/2);
nepochs = size(eD,3);
pow = zeros(Nf,size(eD,1),nepochs);
wind  = window(@hanning,size(eD,2));
coFac = max(wind)/mean(wind);
wind = repmat(wind,1,size(eD,1));

% Calculate PSD for each epoch
for j = 1:nepochs
    Btemp=eD(:,:,j)';
    Btemp = Btemp.*wind*coFac;
    mu=median(Btemp);
    zf = bsxfun(@minus,Btemp,mu);
    fzf = zf;
    N= length(fzf);
    xdft = fft(fzf);
    xdft=xdft(1:floor(N/2+1),:);
    psdx = abs(xdft)./sqrt(N*fs);
    freq = 0:fs/size(fzf,1):fs/2;
    odd=mod(size(fzf,1),2)==1;
    if(odd)
        psdx(2:end) = psdx(2:end);
    else
        psdx(2:end-1) = psdx(2:end-1);
    end
    pow(:,:,j) =psdx;
end
    
% Calculate median or mean out of all epochs
if strcmp(cfg.average,'median')
    po = median(pow(:,:,:),3);
    po_chavg = median(po, 2);
else
    po = mean(pow(:,:,:),3);
    po_chavg = mean(po, 2);
end

%% PLOT

if cfg.plot

    % Set the colormap
    cmap = turbo(size(pow,2));

    % Initialize
    if ~exist('fig', 'var')
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        fig.Name = 'PSD';
    end
    hold on;

    % Plot channels if necessary
    if cfg.channel

        % Plot
        h = plot(freq, po', 'LineWidth', 0.5);

        % Change the color, transparency and add data tips
        for ii = 1:size(pow, 2)
            h(ii).Color(1:3) = cmap(ii, :);
            h(ii).Color(4) = 0.75;
            h(ii).DataTipTemplate.DataTipRows(end+1) = data.label{ii};
        end

    end

    % Plot the average 
    plot(freq, po_chavg, '-k', 'LineWidth', 2);

    % Axis
    set(gca, 'YScale', 'log');
    xlabel('Frequency (Hz)');
    labY = ['$$PSD (' 'fT' ' \sqrt[-1]{Hz}$$)'];
    ylabel(labY,'interpreter','latex')
    xlim(cfg.foi);
    
end
