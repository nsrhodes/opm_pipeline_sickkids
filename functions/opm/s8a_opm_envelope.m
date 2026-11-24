function s8a_opm_envelope(p, sub_num, clobber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s8_opm_tfr(p, sub_num, clobber)
%
% Script to compute a TFR for a particular region
%
% Inputs:
%   opm_study_config    (structure) Study configuration function
%   clobber             (logical)   Whether or not to overwrite existing
%
% Created by: Marlee M Vandewouw
% Created on: July 18th, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MAKE THE FILTERS

% Load the header to get the sampling rate
hdr = ft_read_header(p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'headerformat', 'opm_fil_mjt');
srate = hdr.Fs;

%% LOAD THE TIMESERIES

% Get the timeseries
if strcmp(p.tfr.timesries, 'coordinates')
    ts = load(p.sourcemodel.coordinates.ve(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
elseif strcmp(p.tfr.timesries, 'grid_ts')
    ts = load(p.sourcemodel.grid.ts.ve(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
else
    error('Please specify either ''coordinates'' or ''grid_ts''');
end

% Define some dimensions
num_samples = size(ts.ts, 2);
num_trials = size(ts.ts, 3);
num_sources = size(ts.ts, 1);
time = ts.time;

% Permute to match previous pipeline
ts.ts = permute(ts.ts, [2, 3, 1]);

% Reshape back into continuous for filtering purposes
catmatrix = reshape(ts.ts, [], num_sources);

%% Hilbert envelope over frequency of interest

    % Make the filter coefficients
    freq = str2double(strsplit(p.envelope.freq{sub_num},','));
    freq_avg = mean(freq, 2);

    % Iterate over the number of rois
    for rr = 1:length(p.tfr.rois)

        % Get the region name
        region = ts.pos_name{p.tfr.rois(rr)};


        % Check what timeseries we're using
        if strcmp(p.tfr.timesries, 'coordinates')
            cur_file = p.envelope.ts_output_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all');
            if p.tfr.stats.flag
                stats_file = p.tfr.tfr_coordinates_stats_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all');
            else
                stats_file = cur_file;
            end
        elseif strcmp(p.tfr.timesries, 'grid_ts')
            cur_file = p.envelope.ts_output_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all');
            if p.tfr.stats.flag
                stats_file = p.tfr.tfr_grid_ts_stats_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all');
            else
                stats_file = cur_file;
            end
        end

        if (~exist(cur_file, 'file') || ~exist(stats_file, 'file')) || clobber == 1

            % Extract the ROI ts - source x time x trial
            catmatrix_cur = squeeze(catmatrix(:, p.tfr.rois(rr)));

            % Get the baseline window
            baseline = [-(abs(diff(p.tfr.stats.active))), 0];

            % Extract trials, filter + Hilbert
            H_avg = zeros(size(freq, 1), num_samples);
            H_avg_trials = zeros(size(freq, 1), num_samples, num_trials);
            for ff_idx = 1:size(freq, 1)

                % Extract and mean center
                bldata = bsxfun(@minus, catmatrix_cur, mean(catmatrix_cur,1));

                % Butterworth filter + Hilbert
                filt_VE = nut_filter3(bldata,'butter','bp',4, freq(ff_idx, 1), freq(ff_idx, 2), srate, 1)';
                H_data_tmp = abs(hilbert(filt_VE));

                % Put back into trials, baselining
                H_data = zeros(num_samples, num_trials);
                for tt = 1:num_trials
                    H_data(:,tt) = H_data_tmp(num_samples*(tt-1)+1:num_samples*tt);
                end

                % Average across trials
                H_avg(ff_idx, :) = mean(H_data, 2);

                % Compute the baseline
                H_avg_baseline = mean(H_avg(ff_idx, time >= baseline(1) & time <= baseline(2)));

                % Baseline
                H_avg(ff_idx, :) = (H_avg(ff_idx, :) - H_avg_baseline)./H_avg_baseline;

                % Now, do this per trial if necessary
                if p.tfr.stats.flag

                    % Iterate over the trials
                    for tt = 1:num_trials

                        % Compute the baseline
                        H_avg_baseline = mean(H_data(time >= baseline(1) & time <= baseline(2), tt));

                        % Baseline
                        H_avg_trials(ff_idx, :, tt) = (H_data(:, tt) - H_avg_baseline)./H_avg_baseline;

                    end

                end

            end


            % Make a plot
            figure;
            plot(time,H_avg);
            xlabel('Time (s)');
            ylabel('Relative Amplitude')
            


            % Save
            if strcmp(p.tfr.timesries, 'coordinates')
                if p.tfr.stats.flag
                    save(p.envelope.ts_output_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'), 'H_avg', 'H_avg_trials');
                else
                    save(p.envelope.ts_output_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'), 'H_avg');
                end
            elseif strcmp(p.tfr.timesries, 'grid_ts')
                if p.tfr.stats.flag
                    save(p.envelope.ts_output_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'), 'H_avg', 'H_avg_trials');
                else
                    save(p.envelope.ts_output_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'), 'H_avg');
                end
            end

        end

    end