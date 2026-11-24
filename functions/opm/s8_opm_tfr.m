function s8_opm_tfr(p, sub_num, clobber)
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

%% TFR: HILBERT

if strcmp(p.tfr.method, 'hilbert')

    % Make the filter coefficients
    freq = [[p.tfr.fois(1):p.tfr.fois_ss:p.tfr.fois(2)]', [p.tfr.fois(1)+p.tfr.fois_width:p.tfr.fois_ss:p.tfr.fois(2)+p.tfr.fois_width]'];
    freq = freq(1:end-1, :);
    freq_avg = mean(freq, 2);

    % Iterate over the number of TFRs
    for rr = 1:length(p.tfr.rois)

        % Get the region name
        region = ts.pos_name{p.tfr.rois(rr)};


        % Check what timeseries we're using
        if strcmp(p.tfr.timesries, 'coordinates')
            cur_file = p.tfr.tfr_coordinates_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all');
            if p.tfr.stats.flag
                stats_file = p.tfr.tfr_coordinates_stats_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all');
            else
                stats_file = cur_file;
            end
        elseif strcmp(p.tfr.timesries, 'grid_ts')
            cur_file = p.tfr.tfr_grid_ts_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all');
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

            % Run stats if necessary
            if p.tfr.stats.flag

                % Mean over the active window
                H_active_trials = squeeze(mean(H_avg_trials(:, time >= p.tfr.active(1) & time <= p.tfr.active(2), :), 2));

                % Mean over the baseline window
                H_baseline_trials = squeeze(mean(H_avg_trials(:, time >= baseline(1) & time <= baseline(2), :), 2));

                % Mean over the frequency window
                H_active_trials = squeeze(mean(H_active_trials(freq_avg >= p.tfr.stats.freq(1) & freq_avg <= p.tfr.stats.freq(2), :), 1));
                H_baseline_trials = squeeze(mean(H_baseline_trials(freq_avg >= p.tfr.stats.freq(1) & freq_avg <= p.tfr.stats.freq(2), :), 1));

                % Run a paired t-test
                [h, pval, ci, stats] = ttest(H_active_trials, H_baseline_trials);

                % Make a table with the results
                results = table();
                results.MeanPercentageChange = mean(H_active_trials);
                results.tstat = stats.tstat;
                results.df = stats.df;
                results.sd = stats.sd;
                results.pval = pval;
                results.cohdensd = (mean(H_active_trials) - mean(H_baseline_trials))/stats.sd;

                % Save
                if strcmp(p.tfr.timesries, 'coordinates')
                    writetable(results, p.tfr.tfr_coordinates_stats_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'));
                elseif strcmp(p.tfr.timesries, 'grid_ts')
                    writetable(results, p.tfr.tfr_grid_ts_stats_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'));
                end

            end

            % Make a plot
            figure;
            pcolor(time, freq_avg, H_avg);
            shading interp
            xlabel('Time (s)');
            ylabel('Frequency (Hz)')
            colorbar;
            caxis([-0.2, 0.2])
            axis fill

            % Save
            if strcmp(p.tfr.timesries, 'coordinates')
                saveas(gcf, p.tfr.tfr_coordinates_fig(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'));
                if p.tfr.stats.flag
                    save(p.tfr.tfr_coordinates_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'), 'H_avg', 'H_avg_trials');
                else
                    save(p.tfr.tfr_coordinates_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'), 'H_avg');
                end
            elseif strcmp(p.tfr.timesries, 'grid_ts')
                saveas(gcf, p.tfr.tfr_grid_ts_fig(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'));
                if p.tfr.stats.flag
                    save(p.tfr.tfr_grid_ts_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'), 'H_avg', 'H_avg_trials');
                else
                    save(p.tfr.tfr_grid_ts_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'), 'H_avg');
                end
            end

        end

    end

elseif strcmp(p.tfr.method, 'multitaper')

    % Iterate over the number of TFRs
    for rr = 1:length(p.tfr.rois)

        % Get the region name
        region = ts.pos_name{p.tfr.rois(rr)};

        % Check what timeseries we're using
        if strcmp(p.tfr.timesries, 'coordinates')
            cur_file = p.tfr.tfr_coordinates_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all');
            if p.tfr.stats.flag
                stats_file = p.tfr.tfr_coordinates_stats_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all');
            else
                stats_file = cur_file;
            end
        elseif strcmp(p.tfr.timesries, 'grid_ts')
            cur_file = p.tfr.tfr_grid_ts_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all');
            if p.tfr.stats.flag
                stats_file = p.tfr.tfr_grid_ts_stats_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all');
            else
                stats_file = cur_file;
            end
        end

        if (~exist(cur_file, 'file') || ~exist(stats_file, 'file')) || clobber == 1

            % Extract the ROI ts - source x time x trial
            catmatrix_cur = squeeze(catmatrix(:, p.tfr.rois(rr)));

            % Put into trials
            catmatrix_cur_trials = cell(1, num_trials);
            time_trials = cell(1, num_trials);
            for tt = 1:num_trials
                catmatrix_cur_trials{tt} = catmatrix_cur(num_samples*(tt-1)+1:num_samples*tt)';
                time_trials{tt} = time;
            end

            % Make data structure
            data = {};
            data.hdr = hdr;
            data.label = ts.pos_name(p.tfr.rois(rr));
            data.trial = catmatrix_cur_trials;
            data.fsample = srate;
            data.time = time_trials;

            % Compute TFR
            cfg = [];
            cfg.output     = 'pow';
            cfg.method     = 'mtmconvol';
            cfg.foi        = 4:1:100;
            cfg.t_ftimwin  = 0.2*ones(size(cfg.foi));
            cfg.tapsmofrq  = 10*ones(size(cfg.foi));
            cfg.toi        = time;
            cfg.pad = 'maxperlen';
            cfg.keeptrials = 'yes';
            TFRmult = ft_freqanalysis(cfg, data);
            TFRmult.powspctrm = squeeze(TFRmult.powspctrm);

            % Average over all trials
            H_avg = squeeze(mean(TFRmult.powspctrm, 1));

            % Baseline each frequency
            baseline = [-(abs(diff(p.tfr.stats.active))), 0];
            pow_baseline = mean(H_avg(:, TFRmult.time >= baseline(1) & TFRmult.time <= baseline(2)), 2);
            for ff = 1:length(pow_baseline)
                H_avg(ff, :) = (H_avg(ff, :) - pow_baseline(ff))./pow_baseline(ff);
            end

            % Each indvidual trial
            H_avg_trials = TFRmult.powspctrm;
            for tt = 1:num_trials
                H_avg_cur = squeeze(H_avg_trials(tt, :, :));
                pow_baseline = mean(H_avg_cur(:, TFRmult.time >= baseline(1) & TFRmult.time <= baseline(2)), 2);
                for ff = 1:length(pow_baseline)
                    H_avg_cur(ff, :) = (H_avg_cur(ff, :) - pow_baseline(ff))./pow_baseline(ff);
                end
                H_avg_trials(tt, :, :) = H_avg_cur;
            end
            H_avg_trials = permute(H_avg_trials, [2, 3, 1]);

            % Run stats if necessary
            if p.tfr.stats.flag

                % Mean over the active window
                H_active_trials = squeeze(mean(H_avg_trials(:, time >= p.tfr.active(1) & time <= p.tfr.active(2), :), 2));

                % Mean over the baseline window
                H_baseline_trials = squeeze(mean(H_avg_trials(:, time >= baseline(1) & time <= baseline(2), :), 2));

                % Mean over the frequency window
                H_active_trials = squeeze(mean(H_active_trials(TFRmult.freq >= p.tfr.stats.freq(1) & TFRmult.freq <= p.tfr.stats.freq(2), :), 1));
                H_baseline_trials = squeeze(mean(H_baseline_trials(TFRmult.freq >= p.tfr.stats.freq(1) & TFRmult.freq <= p.tfr.stats.freq(2), :), 1));

                % Run a paired t-test
                [h, pval, ci, stats] = ttest(H_active_trials, H_baseline_trials);

                % Make a table with the results
                results = table();
                results.MeanPercentageChange = mean(H_active_trials);
                results.tstat = stats.tstat;
                results.df = stats.df;
                results.sd = stats.sd;
                results.pval = pval;
                results.cohdensd = (mean(H_active_trials) - mean(H_baseline_trials))/stats.sd;

                % Sort the trials by the percentage change from baseline
                [~, sort_ind] = sort(H_active_trials, 'descend');

                % Save
                if strcmp(p.tfr.timesries, 'coordinates')
                    writetable(results, p.tfr.tfr_coordinates_stats_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'));
                    save(p.tfr.tfr_coordinates_ind_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'), 'sort_ind');
                elseif strcmp(p.tfr.timesries, 'grid_ts')
                    writetable(results, p.tfr.tfr_grid_ts_stats_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'));
                    save(p.tfr.tfr_grid_ts_ind_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'), 'sort_ind');
                end

            end

            % Make a plot
            figure;
            pcolor(time, TFRmult.freq, H_avg);
            shading interp
            xlabel('Time (s)');
            ylabel('Frequency (Hz)')
            colorbar;
            caxis([-0.5, 0.5])
            axis fill

            % Save
            if strcmp(p.tfr.timesries, 'coordinates')
                saveas(gcf, p.tfr.tfr_coordinates_fig(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'));
                if p.tfr.stats.flag
                    save(p.tfr.tfr_coordinates_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'), 'H_avg', 'H_avg_trials');
                else
                    save(p.tfr.tfr_coordinates_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'), 'H_avg');
                end
            elseif strcmp(p.tfr.timesries, 'grid_ts')
                saveas(gcf, p.tfr.tfr_grid_ts_fig(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'));
                if p.tfr.stats.flag
                    save(p.tfr.tfr_grid_ts_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'), 'H_avg', 'H_avg_trials');
                else
                    save(p.tfr.tfr_grid_ts_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'), 'H_avg');
                end
            end

        end

    end

elseif strcmp(p.tfr.method, 'hanning')

    % Iterate over the number of TFRs
    for rr = 1:length(p.tfr.rois)

        % Get the region name
        region = ts.pos_name{p.tfr.rois(rr)};

        % Check what timeseries we're using
        if strcmp(p.tfr.timesries, 'coordinates')
            cur_file = p.tfr.tfr_coordinates_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all');
            if p.tfr.stats.flag
                stats_file = p.tfr.tfr_coordinates_stats_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all');
            else
                stats_file = cur_file;
            end
        elseif strcmp(p.tfr.timesries, 'grid_ts')
            cur_file = p.tfr.tfr_grid_ts_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all');
            if p.tfr.stats.flag
                stats_file = p.tfr.tfr_grid_ts_stats_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all');
            else
                stats_file = cur_file;
            end
        end

        if (~exist(cur_file, 'file') || ~exist(stats_file, 'file')) || clobber == 1

            % Extract the ROI ts - source x time x trial
            catmatrix_cur = squeeze(catmatrix(:, p.tfr.rois(rr)));

            % Put into trials
            catmatrix_cur_trials = cell(1, num_trials);
            time_trials = cell(1, num_trials);
            for tt = 1:num_trials
                catmatrix_cur_trials{tt} = catmatrix_cur(num_samples*(tt-1)+1:num_samples*tt)';
                time_trials{tt} = time;
            end

            % Make data structure
            data = {};
            data.hdr = hdr;
            data.label = ts.pos_name(p.tfr.rois(rr));
            data.trial = catmatrix_cur_trials;
            data.fsample = srate;
            data.time = time_trials;

            % Compute the biggest timewindow
            twin_max = p.tfr.cycles/p.tfr.fois(1);
            timestep = 1./hdr.Fs;
            t1 = time(1) + twin_max/2 + timestep;
            t2 = time(end) - twin_max/2 - timestep;

            % Compute TFR
            cfg = [];
            cfg.output     = 'pow';
            cfg.taper      = 'hanning';
            cfg.method     = 'mtmconvol';
            cfg.foi        = p.tfr.fois(1):1:p.tfr.fois(2);
            cfg.toi        = t1:timestep:t2;
            cfg.t_ftimwin  = cycles./cfg.foi;
            cfg.keeptrials = 'yes';
            TFRmult = ft_freqanalysis(cfg, data);
            TFRmult.powspctrm = squeeze(TFRmult.powspctrm);

            if length(data.trial) > 1
                H_avg = squeeze(mean(TFRmult.powspctrm, 1));
            else
                H_avg = squeeze(TFRmult.powspctrm);
            end

            % Baseline each frequency
            baseline = [-(abs(diff(p.tfr.stats.active))), 0];
            pow_baseline = mean(H_avg(:, TFRmult.time >= baseline(1) & TFRmult.time <= baseline(2)), 2);
            for ff = 1:length(pow_baseline)
                H_avg(ff, :) = (H_avg(ff, :) - pow_baseline(ff))./pow_baseline(ff);
            end

            % Each indvidual trial
            H_avg_trials = TFRmult.powspctrm;
            for tt = 1:num_trials
                H_avg_cur = squeeze(H_avg_trials(tt, :, :));
                pow_baseline = mean(H_avg_cur(:, TFRmult.time >= baseline(1) & TFRmult.time <= baseline(2)), 2);
                for ff = 1:length(pow_baseline)
                    H_avg_cur(ff, :) = (H_avg_cur(ff, :) - pow_baseline(ff))./pow_baseline(ff);
                end
                H_avg_trials(tt, :, :) = H_avg_cur;
            end
            H_avg_trials = permute(H_avg_trials, [2, 3, 1]);

            % Run stats if necessary
            if p.tfr.stats.flag

                % Mean over the active window
                H_active_trials = squeeze(mean(H_avg_trials(:, time >= p.tfr.active(1) & time <= p.tfr.active(2), :), 2));

                % Mean over the baseline window
                H_baseline_trials = squeeze(mean(H_avg_trials(:, time >= baseline(1) & time <= baseline(2), :), 2));

                % Mean over the frequency window
                H_active_trials = squeeze(mean(H_active_trials(TFRmult.freq >= p.tfr.stats.freq(1) & TFRmult.freq <= p.tfr.stats.freq(2), :), 1));
                H_baseline_trials = squeeze(mean(H_baseline_trials(TFRmult.freq >= p.tfr.stats.freq(1) & TFRmult.freq <= p.tfr.stats.freq(2), :), 1));

                % Run a paired t-test
                [h, pval, ci, stats] = ttest(H_active_trials, H_baseline_trials);

                % Make a table with the results
                results = table();
                results.MeanPercentageChange = mean(H_active_trials);
                results.tstat = stats.tstat;
                results.df = stats.df;
                results.sd = stats.sd;
                results.pval = pval;
                results.cohdensd = (mean(H_active_trials) - mean(H_baseline_trials))/stats.sd;

                % Sort the trials by the percentage change from baseline
                [~, sort_ind] = sort(H_active_trials, 'descend');

                % Save
                if strcmp(p.tfr.timesries, 'coordinates')
                    writetable(results, p.tfr.tfr_coordinates_stats_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'));
                    save(p.tfr.tfr_coordinates_ind_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'), 'sort_ind');
                elseif strcmp(p.tfr.timesries, 'grid_ts')
                    writetable(results, p.tfr.tfr_grid_ts_stats_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'));
                    save(p.tfr.tfr_grid_ts_ind_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'), 'sort_ind');
                end

            end

            % Make a plot
            figure;
            pcolor(TFRmult.time, TFRmult.freq, H_avg);
            shading interp
            xlabel('Time (s)');
            ylabel('Frequency (Hz)')
            colorbar;
            caxis([-0.5, 0.5])
            axis fill

            % Save
            if strcmp(p.tfr.timesries, 'coordinates')
                saveas(gcf, p.tfr.tfr_coordinates_fig(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'));
                if p.tfr.stats.flag
                    save(p.tfr.tfr_coordinates_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'), 'H_avg', 'H_avg_trials');
                else
                    save(p.tfr.tfr_coordinates_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'), 'H_avg');
                end
            elseif strcmp(p.tfr.timesries, 'grid_ts')
                saveas(gcf, p.tfr.tfr_grid_ts_fig(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'));
                save(p.tfr.tfr_grid_ts_file(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), region, 'all'), 'H_avg', 'H_avg_trials');
            end

        end

    end

end
