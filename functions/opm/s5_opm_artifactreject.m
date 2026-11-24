function [num_trials, HM] = s5_opm_artifactreject(p, sub_num, clobber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s5_opm_artifactreject(p, sub_num, clobber)
%
% Script to perform artifact rejection on OPM data, by identifying trials
% with amplitudes that exceed a specified threahold.
% 
% Inputs:
%   p                   (structure) Study configuration structure
%   sub_num             (ind)       Subject number to process
%   clobber             (logical)   Whether or not to overwrite existing
%                                   files
%
% Outputs:
%   "trial_cfg.mat"     (file)      FieldTrip trial configuration, with bad
%                                   trials removed
%   num_trials          (array)     Number of trials - first value is total
%                                   number, and each value after
%                                   corresponds to the number of trials per
%                                   condition
%                                   (1 x n+1)
%   HM                  (array)     Mead motion - first value is mean
%                                   across all trials, and each value after
%                                   corresponds to the mean condition
%                                   (1 x n+1)
%
% Created by: Marlee M Vandewouw
% Created on: January 20th, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% INITIALIZATION

if ~exist(p.artifactreject.prev_trl(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file') || clobber == 1 

    % Move file back if necessary
    if exist(p.artifactreject.prev_trl(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file') && clobber == 1
        movefile(p.artifactreject.prev_trl(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), p.epoching.trl_cfg(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
    end

    %% LOAD IN THE DATA

    % Now that all the channels have been identified, load in the bad
    % channels and perform MFC - do this for both the data and the noise
    bad_chs = readcell(p.badchannels.bad_channels(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));

    % Read in the header
    hdr = ft_read_header(p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'headerformat', 'opm_fil_mjt');

    % Load in the MEG data, and demean
    cfg_loadds = [];
    cfg_loadds.dataset = p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num));
    cfg_loadds.headerformat = 'opm_fil_mjt';
    cfg_loadds.dataformat = 'opm_fil_mjt';
    cfg_loadds.eventformat = 'opm_fil_mjt';
    cfg_loadds.channel = hdr.grad.label;
    cfg_loadds.demean = 'yes';

    % Load in the raw data
    data = ft_preprocessing(cfg_loadds);

    %% BAD CHANNELS

    % Eliminate the bad channels
    cfg_ch = [];
    cfg_ch.channel = setdiff(data.label, bad_chs);
    data = ft_selectdata(cfg_ch, data);

    % Also apply to grad
    [~, ind] = intersect(data.grad.label, data.label, 'stable');
    grad_orig = data.grad.label;
    f = fieldnames(data.grad);
    for ii = 1:length(f)
        if length(data.grad.(sprintf('%s', f{ii}))) == length(grad_orig)
            data.grad.(sprintf('%s', f{ii})) = data.grad.(sprintf('%s', f{ii}))(ind, :);
        end
    end

    %% INTERFERENCE SUPRESSION
    
    % Interference supression
    if strcmp(p.interference.method, 'HMC')

        % Mean field correct
        N = data.grad.chanori;
        M = eye(length(N)) - N*pinv(N);
        data.trial{1} = M*data.trial{1};

    elseif strcmp(p.interference.method, 'tSSS')

        % Check if we can perform tSSS
        if ((p.interference.tsss.lin+1)^2+(p.interference.tsss.lout+1)^2-2) > size(data.grad.chanpos,1)
            error('Not enough channels remaining for tSSS!');
        end

        % Temporal signal space separation
        sss_outs = tSSS_opm(data.trial{1}, data.grad.chanpos, data.grad.chanori, p.interference.tsss.lin, p.interference.tsss.lout, 1, p.interference.tsss.corr);
        data.trial{1} = sss_outs.tSSS;

    elseif strcmp(p.interference.method, 'tAMM')
    
        % Check if we can perform tAMM
        % if ((p.interference.tamm.lin+1)^2+(p.interference.tamm.lout+1)^2-2) > size(data.grad.chanpos,1)
        %     error('Not enough channels remaining for tamm!');
        % end
        
        % Temporal signal space separation
        amm_outs = tAMM_opm(data.trial{1}, data.grad.chanpos, data.grad.chanori, p.interference.tamm.lin, p.interference.tamm.lout, p.interference.tamm.corr);
        data.trial{1} = amm_outs.tAMM_data;

    else

        % Error out
        error('Specify either HMC or tSSS for the interference supression method!');

    end

    %% FILTERING

    % Prepare preprocessing for bandpass filtering
    f = fieldnames(p.preproc);
    f = f(startsWith(f, 'bp'));
    cfg_preproc = [];
    for kk = 1:length(f)
        cfg_preproc.(sprintf('%s', f{kk})) = p.preproc.(sprintf('%s', f{kk}));
    end

    % check if bandstop, and determine accordingly
    if isfield(p.preproc, 'bsfilter')            
        if strcmp(p.preproc.bsfilter, 'yes')

            % Add to filter
            cfg_preproc.bsfilter = 'yes';

            if isnumeric(p.preproc.bsfreq)

                % Add the bandstop
                cfg_preproc.bsfreq = [p.preproc.bsfreq'-p.preproc.bswidth, p.preproc.bsfreq'+p.preproc.bswidth];

            elseif strcmp(p.preproc.bsfreq, 'noise_peaks')

                % Load in the peak widths
                cfg_preproc.bsfreq = readmatrix(p.filter.bandstop_widths(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));

            end
        end
    end

    % Now filer
    try
        data_filt = ft_preprocessing(cfg_preproc, data);
    catch
        cfg_preproc.bsfiltord = 3;
        data_filt = ft_preprocessing(cfg_preproc, data);
    end   

    %% EPOCH

    % Load in the trial configuration
    cfg_trl = load(p.epoching.trl_cfg(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));

    % Epoch
    cfg_epoch = [];
    cfg_epoch.trl = cfg_trl.trl;
    data_filt = ft_redefinetrial(cfg_epoch, data_filt);        

    %% ICA        

    % Reject bad ICA componentstrig_offset
    % cfg = [];
    % if exist(p.ica.bad_components(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file')
    %     % Load bad components
    %     cfg.component = load(p.ica.bad_components(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
    %     % Load ICA components
    %     components = load(p.ica.ica_components(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
    %     % Reject
    %     data_filt = ft_rejectcomponent(cfg, components, data_filt);
    % else
    %     error('Subject is missing bad components file.');
    % end

    %% ARTIFACT REJECTION

    if ~exist(p.artifactreject.prev_trl(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file') || clobber == 1

        % Initialize
        include = true(size(cfg_trl.trl, 1), 1);

        % Iterate over the epochs
        for ii = 1:length(include)
            
            % Extract
            dat_trial = data_filt.trial{ii};
            dat_time = data_filt.time{ii};

            % Load in the trial configuration json
            if strcmp(p.epoching.type, 'json')
                fid = fopen(p.epoching.json);
                raw = fread(fid, inf);
                str = char(raw');
                fclose(fid);
                trial_info = jsondecode(str);

                % Find closest sample
                ind = zeros(2, 1);
                for tt = 1:2
                    time_diff = abs(dat_time - trial_info.tInterest(tt));
                    [~, ind(tt)] = min(time_diff);
                end
                
                % Only use tInterest
                dat_trial = dat_trial(:, ind(1):ind(2));
                
            end

            % Look if exceeds within the tInterest
            if any(abs(dat_trial(:)) > p.artifactreject.artf_rej_thresh)
                include(ii) = false;
            end

        end

        % Plot sensor data for all trials
        if ~exist([p.artifactreject.trl_fig(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)) '_1.png'], 'file') || clobber == 1
            trial_plot(data_filt, include, p.artifactreject.artf_rej_thresh, p.artifactreject.trl_fig(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
        end

        % Create new configuration
        cfg_trl.trl_orig = cfg_trl.trl;
        cfg_trl.trl = cfg_trl.trl(include, :);

        % Move old configuration
        movefile(p.epoching.trl_cfg(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), p.artifactreject.prev_trl(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));

        % Save new configuration
        parsave_struct(p.epoching.trl_cfg(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), cfg_trl, 'cfg_trl');

    end
    
end

%% GET NUMBER OF TRIALS

% Load in json if specified
if strcmp(p.epoching.type, 'json')

    % Load in the json file
    fid = fopen(p.epoching.json);
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    trial_info = jsondecode(str);

    % Iterate through and get the number of trials
    num_trials = zeros(1, 1 + length(trial_info.details));

    % Load in the trial configuration
    cfg_trl = load(p.epoching.trl_cfg(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));

    % Assign
    num_trials(1) = size(cfg_trl.trl, 1);
    for ii = 1:length(trial_info.details)
        num_trials(ii+1) = sum(cfg_trl.trl(:, 4) == ii);
    end
    
else

    % Load in the trial configuration
    cfg_trl = load(p.epoching.trl_cfg(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));

    % Assign
    num_trials = size(cfg_trl.trl, 1);
    
end
    

%% GET THE HEAD MOTION

% Only run if it exists
if exist(p.raw_data.headmotion(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file')

    % Read in the header
    hdr = ft_read_header(p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'headerformat', 'opm_fil_mjt');

    % Load in the headmotion file
    HeadMotionData = load(p.raw_data.headmotion(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));

    % Get the head positions
    pos = HeadMotionData.rigidbodies.data(:, contains(HeadMotionData.rigidbodies.colheaders, 'Position'));
    
    % If optitrack data is shorter than opm data, pad with NaNs
    opm_length = max(cfg_trl.trl(:, 2));
    optitrack_length = size(pos, 1);
    if opm_length > optitrack_length
        disp('OPM data recording longer than optitrack recording')
        nan_matrix = NaN(opm_length - optitrack_length, size(pos, 2));
        fprintf('Padding size: %d\n', size(nan_matrix))
        pos = [pos; nan_matrix];
    end
    
    % Iterate over the trials
    HM_trial_max = zeros(size(cfg_trl.trl, 1), 1);

    for tt = 1:size(HM_trial_max, 1)
       
        % Get the head positions for this trial
        pos_trial = pos(cfg_trl.trl(tt, 1):cfg_trl.trl(tt, 2), :);

        % Go back to original sampling
        pos_trial = pos_trial(1:round(hdr.Fs/HeadMotionData.exportFrameRate):end, :);

        % Calculate median head position
        med_pos_trial = nanmedian(pos_trial, 1);

        % Calculate the motion from the median head position for each sample
        motion = sqrt(sum((pos_trial - med_pos_trial).^2, 2))'; %pdist2(med_pos_trial, pos_trial); 

        % Calculate the max head motion, in mm
        HM_trial_max(tt) = max(motion)*1000;

    end
    
    % Get rid of crazy high numbers
    HM_trial_max(HM_trial_max > 100) = NaN;

    % Now, calculate the mean head motion across all trials, and across all
    % trials per condition, if specified
    if strcmp(p.epoching.type, 'json')
        if sum(HM_trial_max(:)) == 0
            HM = nan(1, 1 + length(trial_info.details));
        else
            HM = zeros(1, 1 + length(trial_info.details));
            HM(1) = nanmean(HM_trial_max);
            for ii = 1:length(trial_info.details)
               HM(ii+1) = nanmean(HM_trial_max(cfg_trl.trl(:, 4) == ii));
            end
        end
    else
        if sum(HM_trial_max(:)) == 0
            HM = NaN;
        else
            HM = nanmean(HM_trial_max);
        end
    end
    
else
    if strcmp(p.epoching.type, 'json')
        HM = nan(1, 1 + length(trial_info.details));
    else
        HM = NaN;
    end
    
end
    