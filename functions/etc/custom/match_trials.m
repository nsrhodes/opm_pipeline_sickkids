function idx_closest = match_trials(trials, matching_group)
% Find closest subject based on trial count
trials = trials*ones(size(matching_group,1), 1);
trial_diff = matching_group.trials - trials;

% Select the one that is just above zero (i.e. the one with slightly more
% trials than the base subject)
idx_positive = find(trial_diff>0);
trial_diff_positive = trial_diff(trial_diff>0);
[~, idx] = min(trial_diff_positive);
idx_closest = idx_positive(idx);

% If we don't find one, select the closest one without constraint
if isempty(idx_closest)
    [~, idx_closest] = min(abs(trial_diff));
end
end