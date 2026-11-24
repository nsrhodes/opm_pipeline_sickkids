function data = load_timeseries_data(p, rois, return_trial_mean, extract_top_trials)

% Optional variable logistics
if nargin ~= 3
    return_trial_mean = true;
    extract_top_trials = false;
end

% Get time vector
time = get_time(p);
num_timepoints = length(time);

num_rois = length(rois);
num_subjects = size(p.subject_data, 1);
n_trials = zeros(1, num_subjects);

if return_trial_mean
    % Preallocate memory
    ts = zeros(num_rois, num_timepoints, num_subjects);
    
    % Loop over subjects
    for ss = 1:num_subjects
        % Get response
         fprintf('Loading data for subject %s %s %s\n', p.subject(ss), p.session(ss), p.run(ss))
        [ts(:, :, ss), n_trials(ss)] = get_timeseries(p, rois, ss, true, extract_top_trials);
    end
else
    % Preallocate memory
    ts = cell(1,num_subjects);
    
    % Loop over subjects
    for ss = 1:num_subjects
        % Get response
         fprintf('Loading data for subject %s %s %s\n', p.subject(ss), p.session(ss), p.run(ss))
        [ts{ss}, n_trials(ss)] = get_timeseries(p, rois, ss, false, extract_top_trials);
    end
end

% Create a structure
data.time = time;
data.ts = ts;
data.num_rois = num_rois;
data.num_timepoints = num_timepoints;
data.num_subjects = num_subjects;
data.num_trials = n_trials;
end