function [ts, no_trials] = get_timeseries(p, roi_indices, sub_num, return_trial_mean, extract_top_trials)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_timeseries(p, roi_indices, sub_num)
%
% Function to extract timeseries data after source modeling
%
% Inputs:
%   p                   (structure) Study configuration structure
%   roi_indices         (1xR array) Array containing the indices of the R 
%                                   regions of interest     
%   sub_num             (ind)       Subject number to process
%
% Outputs:
%   time                (1xT array) Array containing the T time samples in 
%                                   seconds 
%   ts                  (RxT array) Matrix containing the timeseries data 
%                                   for every region of interest
%                                   The timeseries data are averaged across
%                                   all trials for the subject, and
%                                   z-scored with respect to a baseline
%                                   interval of ]-0.1, 0[ seconds
%   no_trials           (integer)   The number of trials available for the 
%                                   subject of interest.                                
%                                 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optional variable logistics
if nargin < 4
    return_trial_mean = true;
    extract_top_trials = false;
end

%% LOAD DATA

time = get_time(p);

% Initialize
num_rois = length(roi_indices);
ts = zeros(num_rois, length(time));
if num_rois == 2
    ffg_l = 1; ffg_r = 2;
else
    ffg_l = 55; ffg_r = 56;
end
    
% Load timeseries data
VE = load(p.sourcemodel.grid.ts.ve(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));

% Keep only the regions of interest
ts_all_trials = VE.ts(roi_indices, :, :);

% Discard trials
if extract_top_trials
    ind_l = load([p.directories.sourcemodel_pre(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)) 'VE_coordinates_' p.sourcemodel.grid.ts.atlas_name '_maxpower' num2str(p.sourcemodel.grid.foi(1)) '_' num2str(p.sourcemodel.grid.foi(2)) 'Hz_' num2str(p.sourcemodel.grid.ts.foi(1)) '_' num2str(p.sourcemodel.grid.ts.foi(2)) 'Hz_top_trials_%0.2f_%0.2f_l.mat'], 'ind');
    ind_r = load([p.directories.sourcemodel_pre(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)) 'VE_coordinates_' p.sourcemodel.grid.ts.atlas_name '_maxpower' num2str(p.sourcemodel.grid.foi(1)) '_' num2str(p.sourcemodel.grid.foi(2)) 'Hz_' num2str(p.sourcemodel.grid.ts.foi(1)) '_' num2str(p.sourcemodel.grid.ts.foi(2)) 'Hz_top_trials_%0.2f_%0.2f_r.mat'], 'ind');
    ind_l = ind_l.ind;
    ind_r = ind_r.ind;
    
    ts_sorted = [ts_all_trials(ffg_l, :, ind_l); ts_all_trials(ffg_r, :, ind_r)];

    % Only keep top 50%
    n_trials = size(ts_all_trials, 3);
    n_trials_to_keep = ceil(0.5*n_trials);
    ts_all_trials = ts_sorted(:, :, 1:n_trials_to_keep);
    disp(['Original number of trials:         ', num2str(n_trials)])
    disp(['Number of trials after extraction: ', num2str(n_trials_to_keep)])
end

if return_trial_mean
    % Average accross trials
    ts_temp = mean(ts_all_trials, 3);
    
    % Z-score
    toi = -0.2 < time & time < 0;
    for roi = 1:num_rois
        toi_mean = mean(ts_temp(roi, toi));
        toi_std = std(ts_temp(roi, toi));
        ts(roi, :) = (ts_temp(roi, :) - toi_mean)/toi_std;
    end
else
    ts = ts_all_trials;
end

% Assign number of trials
no_trials = size(ts_all_trials,3);
end
