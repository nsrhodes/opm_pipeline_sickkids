function s4_opm_ica(p, sub_num, clobber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s4_opm_ica(p, sub_num, clobber)
%
% Script to perform ICA on OPM data
%
% Inputs:
%   p                   (structure) Study configuration structure
%   sub_num             (ind)       Subject number to process
%   clobber             (logical)   Whether or not to overwrite existing
%                                   files
%
% Outputs:
%   "ft_ica_comp"       (file)      FieldTrip ICA components from fastica
%   "ica.png"           (files)     Figures of ICA components
%
% Created by: Marlee M Vandewouw
% Created on: January 20th, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZATION

if ~exist(p.ica.ica_components(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file') || clobber == 1

    %% LOAD IN THE DATA

    % Now that all the channels have been identified, load in the bad
    % channels and perform MFC - do this for both the data and the noise
    fid = fopen(p.badchannels.bad_channels(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'r');
    bad_chs = textscan(fid, '%s', 'Delimiter', ',');
    bad_chs = bad_chs{1};
    fclose(fid);

    % Read in the header
    hdr = ft_read_header(p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'headerformat', 'opm_fil_mjt');

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

    %% BAD CHANNELS

    % Eliminate the bad channels
    cfg_ch = [];
    cfg_ch.channel = setdiff(data.label, bad_chs);
    data = ft_selectdata(cfg_ch, data);

    % Also apply to grad
    [~, ind] = intersect(data.grad.label, data.label, 'stable');
    grad_orig = data.grad.label;
    f = fieldnames(data.grad);
    for ii = 1:length(f)
        if length(data.grad.(sprintf('%s', f{ii}))) == length(grad_orig)
            data.grad.(sprintf('%s', f{ii})) = data.grad.(sprintf('%s', f{ii}))(ind, :);
        end
    end

    %% INTERFERENCE SUPRESSION
    
    % Interference supression
    if strcmp(p.interference.method, 'HMC')

        % Mean field correct
        N = data.grad.chanori;
        M = eye(length(N)) - N*pinv(N);
        data.trial{1} = M*data.trial{1};

    elseif strcmp(p.interference.method, 'tSSS')

        % Check if we can perform tSSS
        if ((p.interference.tsss.lin+1)^2+(p.interference.tsss.lout+1)^2-2) > size(data.grad.chanpos,1)
            error('Not enough channels remaining for tSSS!');
        end

        % Temporal signal space separation
        sss_outs = tSSS_opm(data.trial{1}, data.grad.chanpos, data.grad.chanori, p.interference.tsss.lin, p.interference.tsss.lout, 1, p.interference.tsss.corr);
        data.trial{1} = sss_outs.tSSS;

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

            % Add to filter
            cfg_preproc.bsfilter = 'yes';

            if isnumeric(p.preproc.bsfreq)

                % Add the bandstop
                cfg_preproc.bsfreq = [p.preproc.bsfreq'-p.preproc.bswidth, p.preproc.bsfreq'+p.preproc.bswidth];

            elseif strcmp(p.preproc.bsfreq, 'noise_peaks')

                % Load in the peak widths
                cfg_preproc.bsfreq = readmatrix(p.filter.bandstop_widths(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));

            end
        end
    end

    % Now filer
    try
        data_filt = ft_preprocessing(cfg_preproc, data);
    catch
        cfg_preproc.bsfiltord = 3;
        data_filt = ft_preprocessing(cfg_preproc, data);
    end   

    %% EPOCH

    % Load in the trial configuration
    cfg_trl = load(p.epoching.trl_cfg(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));

    % Epoch
    cfg_epoch = [];
    cfg_epoch.trl = cfg_trl.trl;
    data_filt = ft_redefinetrial(cfg_epoch, data_filt);        

    %% ICA

    % Run ICA
    cfg = [];
    cfg.method = 'fastica';
    if isnumeric(p.ica.ica_numcomponents)
        cfg.numcomponent = p.ica.ica_numcomponents;
    elseif strcmp(p.ica.ica_components, 'num_channels')
        cfg.numcomponent = length(data_filt.label);
    end
    comp = ft_componentanalysis(cfg, data_filt);

    % Make plot
    opm_ica_plot(comp, data.grad, @(pp) p.ica.ica_figure(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), pp));    

    % Save
    parsave_struct(p.ica.ica_components(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), comp, 'comp');


end
