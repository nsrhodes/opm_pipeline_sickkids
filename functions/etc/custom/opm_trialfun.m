function [trl, behav_table] = opm_trialfun(cfg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [trl] = opm_trialfun(cfg)
%
% Function to create a FieldTrip trial structure from OPM data that's been
% converted into BIDS format
%
% Inputs:
%   cfg.dataset                 (string)    Path of input data (in BIDS format)
%   cfg.trialdef                (string)    Path to trialdef .json file
%
% JSON structure:
% {
%   "t0marker": ["Drop"],       // marker around which trials begin (ie. t = 0 is at this marker)
%   "t0shift": 0.023,           // shift t=0 by 23 msec to compensate for projector delay
%   "tEpoch": [-2.5, 3.5],      // time window around t=0 to include in each trial
%   "tInterest": [-1.5, 2],     // time window that we are interested in
%   "details": [ ... ]
% }
%
% Details structure:
% {
%   "name": "USUW"             // name of specific trial type / category
% }
%
% Outputs:
%   trl         (array)     FieldTrip trl array
%                           (n x 4)
%                               n       Number of trials
%                               col1    Trial begin (samples)
%                               col2    Trial end   (samples)
%                               offset  # samples before t=0
%                               indx    Index of condition
%
% Created by: Marlee M Vandewouw
% Created on: January 20th, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in the data
cfg_load = [];
cfg_load.headerformat = 'opm_fil_mjt';
cfg_load.dataformat = 'opm_fil_mjt';
cfg_load.eventformat = 'opm_fil_mjt';
cfg_load.dataset = cfg.dataset;
data = ft_preprocessing(cfg_load);

% Load in the json
fid = fopen(cfg.trialdef);
raw = fread(fid, inf);
str = char(raw');
fclose(fid);
trial_info = jsondecode(str);

% Get the trigger offset (in samples) - this is the number of samples
% collected prior to trigger onset (t=0)
trig_offset = ceil(trial_info.tEpoch(1)*data.hdr.Fs);

% Create a default behavioural table
behav_table = table('Size', [length(trial_info.t0marker), 13], 'VariableTypes', [{'string'}; repmat({'double'}, [12, 1])], ...
    'VariableNames', {'Condition', 'Initial', 'Final', 'Correct', 'CorrectPercent', 'CorrectRTMean', 'CorrectRTStd', 'CorrectRTMed', 'Incorrect', 'IncorrectPercent', 'IncorrectRTMean', 'IncorrectRTStd', 'IncorrectRTMed'});

% Iterate over the t0markers
trl = zeros(0, 4);
for tt = 1:length(trial_info.t0marker)
    
    % Get the proper trigger
    trigger = data.trial{1}(strcmp(data.label, trial_info.t0marker{tt}), :);
    
    % Scale
    trigger = trigger/max(trigger);
    
    % Get the trial duration (in samples)
    duration = (trial_info.tEpoch(2) - trial_info.tEpoch(1))*data.hdr.Fs;
    
    % Binarize
    trigger(trigger > 0.9) = 5;
    trigger(trigger < 0.9) = 0;
    [~, trig_start] = findpeaks(trigger,'MinPeakHeight',3,'MinPeakDistance', duration./2);
    
    % Now, figure out if we're actually including (meets the "criteria")
    if isfield(trial_info.details(tt), 'require')
        
        % Determine how many correct/incorrect and include/exclude types
        types = {trial_info.details(tt).require(:).type};
        correct_types = sum(strcmp(types, 'correct') | strcmp(types, 'incorrect'));
        include_types = sum(strcmp(types, 'include') | strcmp(types, 'exclude'));
        
        % Initialize two flags - one for include/exclude, and one for
        % correct/incorrect
        found = nan(length(trig_start), include_types);
        include = nan(length(trig_start), include_types);
        correct = nan(length(trig_start), correct_types);
        rt = nan(length(trig_start), correct_types);
        
        % Initialize the counter
        include_count = 0;
        correct_count = 0;
        
        % Iterate over the requirements
        for rr = 1:length(trial_info.details(tt).require)
            
            % Which type are we?
            if strcmp(trial_info.details(tt).require(rr).type, 'include') || strcmp(trial_info.details(tt).require(rr).type, 'exclude')
                include_count = include_count + 1;
            elseif strcmp(trial_info.details(tt).require(rr).type, 'correct') || strcmp(trial_info.details(tt).require(rr).type, 'incorrect')
                correct_count = correct_count + 1;
            end
            
            % Get the proper trigger
            trigger_require = data.trial{1}(strcmp(data.label, trial_info.details(tt).require(rr).trigName), :);
            
            % Scale
            trigger_require = trigger_require/max(trigger_require);
            
            % Get the trial duration (in samples)
            duration_require = (trial_info.tEpoch(2) - trial_info.tEpoch(1))*data.hdr.Fs;
            
            % Binarize
            trigger_require(trigger_require > 0.9) = 5;
            trigger_require(trigger_require < 0.9) = 0;
            [~, trig_start_require] = findpeaks(trigger_require,'MinPeakHeight',3,'MinPeakDistance', duration_require./2);
            
            % Iterate through the trials
            for trial = 1:length(trig_start)
                
                % Get the search window (in samples)
                sampleSearch = trial_info.details(tt).require(rr).tSearch*data.hdr.Fs;
                
                % Check if the required trigger is within the
                % sampleSearch
                if any((trig_start(trial) + sampleSearch(1) <= (trig_start_require)) & (trig_start(trial) + sampleSearch(2) >= (trig_start_require)))
                    
                    % Find the one trial (first occurance)
                    match_ind = find((trig_start(trial) + sampleSearch(1) <= (trig_start_require)) & (trig_start(trial) + sampleSearch(2) >= (trig_start_require)));
                    match_ind = match_ind(1);
                    
                    
                    % We found the trial! Go through the cases:
                    
                    % - type "include" and we found it
                    if strcmp(trial_info.details(tt).require(rr).type, 'include')
                        found(trial, include_count) = 1;
                        include(trial, include_count) = 1;
                        % - type "exclude" and we found it
                    elseif strcmp(trial_info.details(tt).require(rr).type, 'exclude')
                        found(trial, include_count) = 1;
                        include(trial, include_count) = 0;
                        % - type "correct" and we found it
                    elseif strcmp(trial_info.details(tt).require(rr).type, 'correct')
                        correct(trial, correct_count) = 1;
                        rt(trial, correct_count) = (trig_start_require(match_ind) - trig_start(trial))/data.fsample*1000;
                        % - type "incorrect" and we found it
                    elseif strcmp(trial_info.details(tt).require(rr).type, 'incorrect')
                        correct(trial, correct_count) = 0;
                        rt(trial, correct_count) = (trig_start_require(match_ind) - trig_start(trial))/data.fsample*1000;
                    end
                    
                else
                    
                    % We did not find the trial! Go through the cases:
                    
                    % - type "include" and we didn't find it
                    if strcmp(trial_info.details(tt).require(rr).type, 'include')
                        include(trial, include_count) = 0;
                        % - type "exclude" and we didn't find it
                    elseif strcmp(trial_info.details(tt).require(rr).type, 'exclude')
                        include(trial, include_count) = 0;
                    % - type "correct" and we did not find it, which means
                    % this trial is incorrect and should be excluded
                    elseif strcmp(trial_info.details(tt).require(rr).type, 'correct')
                        correct(trial, correct_count) = 0;
                    % - type "incorrect" and we did not find it, which
                    % means this trial is correct and we should include
                    % it
                    elseif strcmp(trial_info.details(tt).require(rr).type, 'incorrect')
                        correct(trial, correct_count) = 1;
                    end
                    
                end
                
            end
            
            
        end
        
        % Get the included indices
        ind_found = sum(found, 2) == include_count;
        ind_include = sum(include, 2) == include_count;
        ind_correct = sum(correct, 2) == correct_count;
        
        % Assign the variables
        behav_table.Condition(tt) = trial_info.details(tt).name;
        behav_table.Initial(tt) = sum(ind_found);
        behav_table.Final(tt) = sum(ind_include & ind_correct);
        behav_table.Correct(tt) = sum(ind_found & ind_correct);
        behav_table.CorrectPercent(tt) = sum(ind_found & ind_correct)/sum(ind_found)*100;
        behav_table.CorrectRTMean(tt) = mean(rt(ind_found & ind_correct));
        behav_table.CorrectRTStd(tt) = std(rt(ind_found & ind_correct));
        behav_table.CorrectRTMed(tt) = median(rt(ind_found & ind_correct));
        behav_table.Incorrect(tt) = sum(ind_found & ~ind_correct);
        behav_table.IncorrectPercent(tt) = sum(ind_found & ~ind_correct)/sum(ind_found)*100;
        behav_table.IncorrectRTMean(tt) = mean(rt(ind_found & ~ind_correct));
        behav_table.IncorrectRTStd(tt) = std(rt(ind_found & ~ind_correct));
        behav_table.IncorrectRTMed(tt) = median(rt(ind_found & ~ind_correct));
        
    else
        
        % Assume we're including everything
        ind_found = true(length(trig_start), 1);
        ind_include = true(length(trig_start), 1);
        ind_correct = true(length(trig_start), 1);
        
        % Assign the variables
        behav_table.Condition(tt) = trial_info.details(tt).name;
        behav_table.Initial(tt) = sum(ind_found);
        behav_table.Final(tt) = sum(ind_include & ind_correct);
        
    end
    
    % Get rid of trials
    trig_start = trig_start(ind_include & ind_correct);
    
    % Adjust for trigger offset
    trig_start = trig_start + trig_offset + ceil(trial_info.t0shift*data.hdr.Fs);
    
    % Compute end time
    trig_end = trig_start + duration - 1;
    
    % Add to temporary table
    trl_temp = [trig_start', trig_end', ones(size(trig_start))'.*trig_offset, repmat(tt, length(trig_start), 1)];
    
    % Add to table
    trl = cat(1, trl, trl_temp);
    
end

