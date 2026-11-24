function s3_opm_filtering(p, sub_num, clobber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s3_opm_fil_mjttering(p, sub_num, clobber)
%
% Script to read in OPM MEG data, perform bandpass filtering, and perform
% notch filtering. For notch filtering, if 'noise_peaks' is specified, the
% noise data is used to identify peaks that are present in both the noise 
% and data, and the width of the peak to notch filter is the smallest 
% possible for the peak to be no longer identified at half prominence.
%
% Inputs:
%   p                   (structure) Study configuration structure
%   sub_num             (ind)       Subject number to process
%   clobber             (logical)   Whether or not to overwrite existing
%                                   files
%
% Outputs:
%   "widths"            (file)      Widths of bandstop filters identified
%   "psd_figure"        (file)      Figure showing pre- and post-filter PSD
%
% Created by: Marlee M Vandewouw
% Created on: January 20th, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZATION

if ~exist(p.filter.bandstop_widths(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file') || clobber == 1

    %% LOAD IN THE DATA

    % Now that all the channels have been identified, load in the bad
    % channels and perform MFC - do this for both the data and the noise
    bad_chs = readcell(p.badchannels.bad_channels(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));

    % Read in the header
    hdr = ft_read_header(p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'headerformat', 'opm_fil_mjt');
    hdr_noise = ft_read_header(p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.noise_run(sub_num), 'noise'), 'headerformat', 'opm_fil_mjt');
    
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
    
    % Update the channel
    cfg_loadds.channel = hdr_noise.grad.label;

    % Load in the noise
    cfg_loadds.dataset = p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.noise_run(sub_num), 'noise');
    noise = ft_preprocessing(cfg_loadds);

    %% BAD CHANNELS

    % Eliminate the bad channels
    cfg_ch = [];
    cfg_ch.channel = setdiff(data.label, bad_chs);
    data = ft_selectdata(cfg_ch, data);
    noise = ft_selectdata(cfg_ch, noise); 

    % Also apply to grad
    [~, ind] = intersect(data.grad.label, data.label, 'stable');
    [~, ind_noise] = intersect(noise.grad.label, noise.label, 'stable');
    grad_orig = data.grad.label;
    grad_orig_noise = noise.grad.label;
    f = fieldnames(data.grad);
    for ii = 1:length(f)
        if length(data.grad.(sprintf('%s', f{ii}))) == length(grad_orig)
            data.grad.(sprintf('%s', f{ii})) = data.grad.(sprintf('%s', f{ii}))(ind, :);
        end
    end
    f_noise = fieldnames(noise.grad);
    for ii = 1:length(f_noise)
        if length(noise.grad.(sprintf('%s', f_noise{ii}))) == length(grad_orig_noise)
            noise.grad.(sprintf('%s', f_noise{ii})) = noise.grad.(sprintf('%s', f_noise{ii}))(ind_noise, :);
        end
    end  

    %% INTERFERENCE SUPRESSION
    
     % Interference supression
    if strcmp(p.interference.method, 'HMC')

        % Mean field correct
        N = data.grad.chanori;
        M = eye(length(N)) - N*pinv(N);
        data.trial{1} = M*data.trial{1};
        N_noise = noise.grad.chanori;
        M_noise = eye(length(N_noise)) - N_noise*pinv(N_noise);
        noise.trial{1} = M_noise*noise.trial{1};

    elseif strcmp(p.interference.method, 'tSSS')

        % Check if we can perform tSSS
        if ((p.interference.tsss.lin+1)^2+(p.interference.tsss.lout+1)^2-2) > size(data.grad.chanpos,1)
            error('Not enough channels remaining for tSSS!');
        end
        
        % Some cases where the channel orientation/position for the noise
        % data is all zeros
        if any(sum(noise.grad.chanpos, 2) == 0) && any(sum(noise.grad.chanori, 2) == 0)
            % Make sure that the rest of the orientation/positions are
            % iequal
            if isequal(noise.grad.chanpos((sum(noise.grad.chanpos, 2) ~= 0), :), data.grad.chanpos((sum(noise.grad.chanpos, 2) ~= 0), :))
                noise.grad.chanpos = data.grad.chanpos;
                noise.grad.chanori = data.grad.chanori;
            else
                error('Problem with the noise orientation and positions!');
            end
        end

        % Temporal signal space separation
        sss_outs = tSSS_opm(data.trial{1}, data.grad.chanpos, data.grad.chanori, p.interference.tsss.lin, p.interference.tsss.lout, 1, p.interference.tsss.corr);
        data.trial{1} = sss_outs.tSSS;
        sss_outs_noise = tSSS_opm(noise.trial{1}, noise.grad.chanpos, noise.grad.chanori, p.interference.tsss.lin, p.interference.tsss.lout, 1, p.interference.tsss.corr);
        noise.trial{1} = sss_outs_noise.tSSS;
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

    % Check if bandstop, and determine accordingly
    if isfield(p.preproc, 'bsfilter')
        if strcmp(p.preproc.bsfilter, 'yes')
            if ~isnumeric(p.preproc.bsfreq) && strcmp(p.preproc.bsfreq, 'noise_peaks')

                % Do the bandpass filtering
                data_bp = ft_preprocessing(cfg_preproc, data);
                noise_bp = ft_preprocessing(cfg_preproc, noise);

                % Set some values
                cfg_psd = [];
                cfg_psd.foi = p.preproc.bpfreq;
                cfg_psd.channel = true;
                cfg_psd.trial_length = 1;
                cfg_psd.average = 'median';

                % Do the filtering baed on the peaks
                [~, peakWidths] = opm_filter_noise_peaks(data_bp, noise_bp, p.filtering.fig(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), cfg_psd);

                % Save the peak widths
                writematrix(peakWidths, p.filter.bandstop_widths(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));

            end
        end
    end

end