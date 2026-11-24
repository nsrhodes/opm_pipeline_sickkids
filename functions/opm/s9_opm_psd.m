function s9_opm_psd(p, sub_num, clobber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s9_opm_fooof(p, sub_num, clobber)
%
% Script to run FOOOF on OPM data
%
% Inputs:
%   opm_study_config    (structure) Study configuration function
%   clobber             (logical)   Whether or not to overwrite existing
%
% Created by: Marlee M Vandewouw
% Created on: October 30, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% VE TIMESERIES: ATLAS

% Only run if this was specified
if contains(p.sourcemodel.method, 'coordinates')
    
    % Regional
    
    % Only run if doesn't exist or if clobbering
    if ~exist(p.fooof.psd_mat_ve_coordinates(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file') || clobber == 1
        
        % Load in the timeseries
        ts = load(p.sourcemodel.coordinates.ve(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
        srate = ts.srate;
        
        % Define some dimensions
        num_sources = size(ts.ts, 1);
        
        % Permute to match previous pipeline
        ts.ts = permute(ts.ts, [2, 3, 1]);
        
        % Reshape back into continuous for filtering purposes
        ts = reshape(ts.ts, [], num_sources);
        
        % Run to get dimensions
        [~, f] = pwelch(zscore(ts(:, 1)), 1*srate, [], [], srate, 'psd');

        % Iterate over sources
        psd = zeros(num_sources, length(f));
        for rr = 1:num_sources
            
            % use pwelch to get power spectral density estimate
            [pxx, ~] = pwelch(zscore(ts(:, rr)), 1*srate, [], [], srate, 'psd'); % last input parameter is sampling rate in Hz
            psd(rr, :) = pxx;
            
        end
        
        % Save
        parsave_twovar(p.fooof.psd_mat_ve_coordinates(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), psd, 'psd', f, 'f');
        
    end
    
    % Whole brain
    if p.fooof.whole_brain
        
         % Only run if doesn't exist or if clobbering
        if ~exist(p.fooof.psd_mat_ve_coordinates_whole_brain(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file') || clobber == 1

            % Load in the timeseries
            ts = load(p.sourcemodel.coordinates.ve(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
            srate = ts.srate;

            % Define some dimensions
            num_sources = size(ts.ts, 1);

            % Permute to match previous pipeline
            ts.ts = permute(ts.ts, [2, 3, 1]);

            % Reshape back into continuous for filtering purposes
            ts = reshape(ts.ts, [], num_sources);
            
            % Reshape into a single timeseries
            ts = reshape(ts, [], 1);

            % Run to get dimensions
            [psd, f] = pwelch(zscore(ts), 1*srate, [], [], srate, 'psd');
            % [psd, f] = pmtm(zscore(ts), 3, length(ts), srate);

            % Save
            parsave_twovar(p.fooof.psd_mat_ve_coordinates_whole_brain(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), psd, 'psd', f, 'f');

        end

    end
    
end

