function s7a_opm_connectivity_over_trials_wpli(p, sub_num, clobber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s7_opm_connectivity_wpli(p, sub_num, clobber)
%
% Script to compute wPLI connectivity matrices from OPM data
%
% Inputs:
%   opm_study_config    (structure) Study configuration function
%   clobber             (logical)   Whether or not to overwrite existing
%
% Created by: Marlee M Vandewouw
% Created on: April 24th, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZATION

% wPLI conenctivity function
connfn = @(H1,H2) abs(mean(sin(angle(H1)-angle(H2)), 2));

%% MAKE FILTERS

% Load the header to get the sampling rate
hdr = ft_read_header(p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'headerformat', 'opm_fil_mjt');
srate = hdr.Fs;

% Make the filter coefficients
maxn = 0;
fir_coef = cell(size(p.connectivity.fois, 1), 1);
for ff_idx = 1:size(p.connectivity.fois, 1)
    
    % Current frequency
    filtband = p.connectivity.fois(ff_idx, :);
    
    % Estimate the filter order
    freqs = [filtband(1) * 0.6, filtband(1), filtband(2), filtband(2) * 1.9];
    fir_n = kaiserord(freqs, [0 1 0], [0.25 0.05 0.25], srate);
    maxn = max(maxn, fir_n);
    
    % Make the filter coefficients
    fir_coef{ff_idx} = fir1(fir_n, filtband*2/srate, 'bandpass');
    
end

%% VE TIMESERIES: ATLAS

% Only run if this was specified
if contains(p.sourcemodel.method, 'coordinates')
    
    % Only run if doesn't exist or if clobbering
    if ~exist(p.conenctivity.wpli_mat_ve_coordinates(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), 'all', num2str(size(p.connectivity.fois, 1))), 'file') || clobber == 1
        
        % Load in the timeseries
        ts = load(p.sourcemodel.coordinates.ve(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
        
        % Define some dimensions
        num_samples = size(ts.ts, 2);
        num_trials = size(ts.ts, 3);
        num_sources = size(ts.ts, 1);
        
        % Permute to match previous pipeline
        ts.ts = permute(ts.ts, [2, 3, 1]);
        
        % Reshape back into continuous for filtering purposes
        catmatrix = reshape(ts.ts, [], num_sources);
        
        % Initialize: always across all conditions, optionally for each
        % condition
        allfq_adjmat = nan(num_sources, num_sources, num_samples, size(p.connectivity.fois, 1));
        if p.connectivity.per_cond
            % Load in the trial configuration
            cfg_trl = load(p.epoching.trl_cfg(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
            cond_ind = unique(cfg_trl.trl(:, 4));
            allfq_adjmat_percond = nan(num_sources, num_sources, num_samples, size(p.connectivity.fois, 1), length(cond_ind));
        end
        
        % Extract trials, filter + Hilbert
        for ff_idx = 1:size(p.connectivity.fois, 1)
            
            % Extract and mean center
            bldata = squeeze(catmatrix);
            bldata = bsxfun(@minus, bldata, mean(bldata,1));
            
            % Take the hilbert - only filter if the frequency doesn't equal
            % the sourmoceling frequency
            if ~isequal(p.connectivity.fois(ff_idx, :), p.sourcemodel.coordinates.foi)
                H_data_tmp = hilbert(filtfilt(fir_coef{ff_idx}, 1, bldata));
            else
                H_data_tmp = hilbert(bldata);
            end
            
            % Put back into trials
            H_data = complex(zeros(num_samples, num_trials, num_sources));
            for tt = 1:num_trials
                H_data(:,tt,:) = H_data_tmp(num_samples*(tt-1)+1:num_samples*tt,:);
            end
        
            % Initialize connectivity matrix
            p_adjmat = nan(num_sources, num_sources, num_samples);
            if p.connectivity.per_cond
                p_adjmat_percond = nan(num_sources, num_sources, num_samples, length(cond_ind));
            end
            
            % Calculate connectivity
            for aa = 1:num_sources
                for bb = aa+1:num_sources
                    p_adjmat(aa,bb,:) = connfn(H_data(:,:,aa), H_data(:,:,bb));
                    p_adjmat(bb,aa,:) = p_adjmat(aa,bb,:);
                    if p.connectivity.per_cond
                        for cc = 1:length(cond_ind)
                            cond_trials = cfg_trl.trl(:, 4) == cond_ind(cc);
                            p_adjmat_percond(aa, bb, :, cc) = connfn(H_data(:,cond_trials,aa), H_data(:,cond_trials,bb));
                            p_adjmat_percond(bb,aa,:,cc) = p_adjmat_percond(aa,bb,:,cc);
                        end
                    end
                end
            end
            
            % Save
            if p.connectivity.per_cond
                for cc = 1:length(cond_ind)
                    adjmat = p_adjmat_percond(:, :, :, cc);
                    parsave(p.conenctivity.wpli_mat_ve_coordinates(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), sprintf('c%d', cc), sprintf('%d', ff_idx)), adjmat, 'adjmat');
                end
            end
            adjmat = p_adjmat;
            parsave(p.conenctivity.wpli_mat_ve_coordinates(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), 'all', sprintf('%d', ff_idx)), adjmat, 'adjmat');
            
        end
        
    end
    
end        

%% VE TIMESERIES: GRID

% Only run if this was specified
if contains(p.sourcemodel.method, 'grid') && p.sourcemodel.grid.ts.flag
    
    % Only run if doesn't exist or if clobbering
    if ~exist(p.conenctivity.wpli_mat_ve_grid(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), 'all', num2str(size(p.connectivity.fois, 1))), 'file') || clobber == 1
        
        % Load in the timeseries
        ts = load(p.sourcemodel.grid.ts.ve(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
        
        % Define some dimensions
        num_samples = size(ts.ts, 2);
        num_trials = size(ts.ts, 3);
        num_sources = size(ts.ts, 1);
        
        % Permute to match previous pipeline
        ts.ts = permute(ts.ts, [2, 3, 1]);
        
        % Reshape back into continuous for filtering purposes
        catmatrix = reshape(ts.ts, [], num_sources);
        
        % Initialize: always across all conditions, optionally for each
        % condition
        allfq_adjmat = nan(num_sources, num_sources, num_samples, size(p.connectivity.fois, 1));
        if p.connectivity.per_cond
            % Load in the trial configuration
            cfg_trl = load(p.epoching.trl_cfg(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
            cond_ind = unique(cfg_trl.trl(:, 4));
            allfq_adjmat_percond = nan(num_sources, num_sources, num_samples, size(p.connectivity.fois, 1), length(cond_ind));
        end
        
        % Extract trials, filter + Hilbert
        for ff_idx = 1:size(p.connectivity.fois, 1)
            
            % Extract and mean center
            bldata = squeeze(catmatrix);
            bldata = bsxfun(@minus, bldata, mean(bldata,1));
            
            % Take the hilbert - only filter if the frequency doesn't equal
            % the sourmoceling frequency
            if ~isequal(p.connectivity.fois(ff_idx, :), p.sourcemodel.grid.ts.foi)
                H_data_tmp = hilbert(filtfilt(fir_coef{ff_idx}, 1, bldata));
            else
                H_data_tmp = hilbert(bldata);
            end
            
            % Put back into trials
            H_data = complex(zeros(num_samples, num_trials, num_sources));
            for tt = 1:num_trials
                H_data(:,tt,:) = H_data_tmp(num_samples*(tt-1)+1:num_samples*tt,:);
            end
        
            % Initialize connectivity matrix
            p_adjmat = nan(num_sources, num_sources, num_samples);
            if p.connectivity.per_cond
                p_adjmat_percond = nan(num_sources, num_sources, num_samples, length(cond_ind));
            end
            
            % Calculate connectivity
            for aa = 1:num_sources
                for bb = aa+1:num_sources
                    p_adjmat(aa,bb,:) = connfn(H_data(:,:,aa), H_data(:,:,bb));
                    p_adjmat(bb,aa,:) = p_adjmat(aa,bb,:);
                    if p.connectivity.per_cond
                        for cc = 1:length(cond_ind)
                            cond_trials = cfg_trl.trl(:, 4) == cond_ind(cc);
                            p_adjmat_percond(aa, bb, :, cc) = connfn(H_data(:,cond_trials,aa), H_data(:,cond_trials,bb));
                            p_adjmat_percond(bb,aa,:,cc) = p_adjmat_percond(aa,bb,:,cc);
                        end
                    end
                end
            end
            
            % Save
            if p.connectivity.per_cond
                for cc = 1:length(cond_ind)
                    adjmat = p_adjmat_percond(:, :, :, cc);
                    parsave(p.conenctivity.wpli_mat_ve_grid(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), sprintf('c%d', cc), sprintf('%d', ff_idx)), adjmat, 'adjmat');
                end
            end
            adjmat = p_adjmat;
            parsave(p.conenctivity.wpli_mat_ve_grid(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), 'all', sprintf('%d', ff_idx)), adjmat, 'adjmat');
            
        end
        
    end
    
end 
