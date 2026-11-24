close all
clear all
clc
%% 
% Run the study configuration
p = opm_study_config();

% Overwrite?
clobber =0;

%% EPOCHING

% Iterate over each row of the subject info 
for ss = 1:size(p.subject_data, 1)

    fprintf('Epoching subject: %s, session: %s, task: %s, run: %s\n', p.subject(ss), p.session(ss), p.task(ss), p.task(ss));
    % Call epoching function
    s1_opm_epoching(p, ss, clobber);

end

% Get counts of original trials
num_trials_orig = zeros(size(p.subject_data, 1), 1);
for ss = 1:size(p.subject_data, 1)
    
    % Load in the trial configuration
    cfg_trl = load(p.epoching.trl_cfg(p.subject(ss), p.session(ss), p.run(ss), p.task(ss)));
    num_trials_orig(ss) = size(cfg_trl.trl, 1);
    
end

%% BAD CHANNELS

for ss = 1:size(p.subject_data, 1)

    fprintf('Identifying bad channels for subject: %s, session: %s, task: %s, run: %s\n', p.subject(ss), p.session(ss), p.task(ss), p.task(ss));
    s2_opm_badchannels(p, ss, clobber);

end

% Get the number of bad channels
num_bad_channels = zeros(size(p.subject_data, 1), 1);
for ss = 1:size(p.subject_data, 1)
    bad_chs = readcell(p.badchannels.bad_channels(p.subject(ss), p.session(ss), p.run(ss), p.task(ss)));
    num_bad_channels(ss) = length(bad_chs);
end


%% FILTERING

% Iterate over each row of the subject info 
for ss = 1:size(p.subject_data, 1)

    fprintf('Filtering subject: %s, session: %s, task: %s, run: %s\n', p.subject(ss), p.session(ss), p.task(ss), p.task(ss));
    s3_opm_filtering(p, ss, 0);

    close all;
end

peaks = cell(size(p.subject_data, 1), 1);
for ss = 1:size(p.subject_data, 1)
    bsfreq = readmatrix(p.filter.bandstop_widths(p.subject(ss), p.session(ss), p.run(ss), p.task(ss)));
    peaks{ss} = strrep(strrep(strrep(mat2str(mean(bsfreq, 2)), ';', ','), ']', ''), '[', '');
    num_peaks(ss) = length(bsfreq);
end


%% ARTIFACT REJECT
% Iterate over each row of the subject info 
% Initialize number of trials and mean head motion

num_trials = zeros(size(p.subject_data, 1), 3);
HM = zeros(size(p.subject_data, 1), 3);
for ss = 1:size(p.subject_data, 1) 

    fprintf('Artifact rejection for subject: %s, session: %s, task: %s, run: %s\n', p.subject(ss), p.session(ss), p.task(ss), p.run(ss));
    [num_trials(ss, :), HM(ss, :)] = s5_opm_artifactreject(p, ss, clobber);

    close all
end

%% ICA
% Iterate over each row of the subject info 
for ss = 1:size(p.subject_data, 1)

    if p.subject_data.run_ica(ss) == 1

        fprintf('ICA for subject: %s, session: %s, task: %s, run: %s\n', p.subject(ss), p.session(ss), p.task(ss), p.task(ss));

        s5a_opm_ica(p,ss,clobber);

    end
end

close all

%% RECONFIGURE

% Re-call the configuration script to create the ICA bad components file
p = opm_study_config();

% Reset epoch/trial config if additional changed were made
%p.epoching.trl_cfg = @(sub,ses,run,task)[p.directories.epoch_pre(sub,ses,run,task),'trials_cfg_updated.mat'];
    
clobber = 0;


%% SOURCE MODELLING
% Iterate over each row of the subject 
% Use parallel processing if available with parfor

for ss = 1:size(p.subject_data, 1)
    
    s6_opm_sourcemodel(p, ss, clobber);
 
end

%% CONNECTIVITY

% Iterate over each row of the subject info
parfor ss = 1:size(p.subject_data, 1)

    fprintf('Computing wPLI connectivity for subject: %s, session: %s, task: %s, run: %s\n', p.subject(ss), p.session(ss), p.task(ss), p.task(ss));
    % Call wPLI connectivity
    s7a_opm_connectivity_over_trials_wpli(p, ss, 1);

end

%Iterate over each row of the subject info
parfor ss = 1:size(p.subject_data, 1)

    fprintf('Computing AEC connectivity for subject: %s, session: %s, task: %s, run: %s\n', p.subject(ss), p.session(ss), p.task(ss), p.task(ss));
    % Call AEC connectivity
    s7b_opm_connectivity_over_time_aec(p, ss, 1);

end

%% TIME-FREQUENCY

% Iterate over each row of the subject info
for ss = 1:size(p.subject_data, 1)

    fprintf('Computing a TFR for subject: %s, task: %s\n', p.subject(ss), p.task(ss));
    % Call source modeling function
    s8_opm_tfr(p, ss, clobber);

end
%% Envelope over frequency

% Iterate over each row of the subject info
for ss = 1:size(p.subject_data, 1)

    fprintf('Computing envelope for subject: %s, task: %s\n', p.subject(ss), p.task(ss));

    % Call source modeling function
    s8a_opm_envelope(p, ss, clobber);

end
