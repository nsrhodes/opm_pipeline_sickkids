function s6_opm_sourcemodel_pca(p, sub_num, clobber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% s6_opm_sourcemodel(p, sub_num, clobber)
%
% Script to perform inverse modeling on OPM data to reconstruct sources.
% Depending on the specifications, this script will either (1) perform
% source reconstruction on pre-specified coordinates, (2) perform source
% reconstruction on a grid, or (3) both. If grid reconstruction is
% performed: (a) nifti images will be created containing the average power
% of every voxel in the grid across a time window of interest (computed as
% percentage change from baseline), and/or (b) for a given parcellation,
% the voxel will the largest percentage increase in power compared to
% baseline across a given time window of interest will be identified, and
% its timeseries will be extracted.
%
% Inputs:
%   opm_study_config    (structure) Study configuration function
%   clobber             (logical)   Whether or not to overwrite existing
%                                   files
%
% Outputs:
%   "W"                 (file)      Beamformer weights
%   (optional)
%   "VE_ts"             (file)      Structure containing timeseries from
%                                   pre-specified coordinates, along with
%                                   other pertinent info
%   "nii"               (file)      Nifti images containing the percentage
%                                   change in power across a time window
%   "VE_grid"           (file)      Structure containing timeseries from
%                                   coordinates identified as parcel
%                                   voxcel's with maximal percentage change
%                                   in power, along with other pertinent
%                                   info
%
% Created by: Marlee M Vandewouw
% Created on: January 20th, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZATION

% Load in the trial configuration json
if strcmp(p.epoching.type, 'json')
    fid = fopen(p.epoching.json);
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    trial_info = jsondecode(str);
end

%% CHECK IF FILES EXIST OR OVERWRITING

% Check if all the weight files exist
ex = [false, false];
if contains(p.sourcemodel.method, 'grid')
    if exist(p.sourcemodel.grid.weights(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file') && clobber == 0
        ex(1) = true;
        if p.sourcemodel.grid.avgpow.flag && ~exist([p.sourcemodel.grid.avgpow.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.avgpow.twoi(end, 1), p.sourcemodel.grid.avgpow.twoi(end, 2), 'all') '.nii.gz'], 'file')
            ex(1) = false;
        end
        if p.sourcemodel.grid.ts.flag && ~exist(p.sourcemodel.grid.ts.ve(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file')
            ex(1) = false;
        end
        if p.sourcemodel.grid.seed.flag
            % Load in the coordinate set info
            if ~exist(p.sourcemodel.grid.seed.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.seed.names{end}), 'file')
                ex(1) = false;
            end
            if p.sourcemodel.grid.seed.reg_mni
                if ~exist(p.sourcemodel.grid.seed.nii_mni(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.seed.names{end}), 'file')
                    ex(1) = false;
                end
            end
        end
    end
else
    ex(1) = true;
end
if contains(p.sourcemodel.method, 'coordinates')
    if exist(p.sourcemodel.coordinates.weights(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file') && clobber == 0
        ex(2) = true;
        if ~exist(p.sourcemodel.coordinates.ve(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file')
            ex(2) = false;
        end
    end
else
    ex(2) = true;
end

% Only run if we need to
if ~all(ex)
    
    %% LOAD IN THE DATA
    
    % Now that all the channels have been identified, load in the bad
    % channels and perform MFC - do this for both the data and the noise
    bad_chs = readcell(p.badchannels.bad_channels(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
    
    % Read in the header
    hdr = ft_read_header(p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'headerformat', 'opm_fil_mjt');
    
    % Load in the MEG data, CONTINUOUS, and demean
    cfg_loadds = [];
    cfg_loadds.dataset = p.raw_data.bin(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num));
    cfg_loadds.headerformat = 'opm_fil_mjt';
    cfg_loadds.dataformat = 'opm_fil_mjt';
    cfg_loadds.eventformat = 'opm_fil_mjt';
    cfg_loadds.channel = hdr.grad.label;
    cfg_loadds.demean = 'yes';
    cfg_loadds.continuous = 'yes';
    
    % Load in the raw data
    data_continuous = ft_preprocessing(cfg_loadds);
    
    %% BAD CHANNELS
    
    % Eliminate the bad channels
    cfg_ch = [];
    cfg_ch.channel = setdiff(data_continuous.label, bad_chs);
    data_continuous = ft_selectdata(cfg_ch, data_continuous);
    
    % Also apply to grad
    [~, ind] = intersect(data_continuous.grad.label, data_continuous.label, 'stable');
    grad_orig = data_continuous.grad.label;
    f = fieldnames(data_continuous.grad);
    for ii = 1:length(f)
        if length(data_continuous.grad.(sprintf('%s', f{ii}))) == length(grad_orig)
            data_continuous.grad.(sprintf('%s', f{ii})) = data_continuous.grad.(sprintf('%s', f{ii}))(ind, :);
        end
    end
    
    %% INTERFERENCE SUPRESSION
    
    % Interference supression
    if strcmp(p.interference.method, 'HMC')
        
        % Mean field correct
        N = data_continuous.grad.chanori;
        M = eye(length(N)) - N*pinv(N);
        data_continuous.trial{1} = M*data_continuous.trial{1};
        
    elseif strcmp(p.interference.method, 'tSSS')
        
        % Check if we can perform tSSS
        if ((p.interference.tsss.lin+1)^2+(p.interference.tsss.lout+1)^2-2) > size(data_continuous.grad.chanpos,1)
            error('Not enough channels remaining for tSSS!');
        end
        
        % Temporal signal space separation
        sss_outs = tSSS_opm(data_continuous.trial{1}, data_continuous.grad.chanpos, data_continuous.grad.chanori, p.interference.tsss.lin, p.interference.tsss.lout, 1, p.interference.tsss.corr);
        data_continuous.trial{1} = sss_outs.tSSS;
        
        % Create dummy matrix
        M = eye(length(data_continuous.grad.chanori));
        
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
        data_continuous_filt = ft_preprocessing(cfg_preproc, data_continuous);
    catch
        cfg_preproc.bsfiltord = 3;
        data_continuous_filt = ft_preprocessing(cfg_preproc, data_continuous);
    end
    
    %% ICA
    
    % Reject bad ICA components
    cfg = [];
    if exist(p.ica.bad_components(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file')
        % Load bad components
        cfg.component = load(p.ica.bad_components(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
        % Load ICA components
        components = load(p.ica.ica_components(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
        % Reject
        data_continuous_filt = ft_rejectcomponent(cfg, components, data_continuous_filt);
    else
        %error('Subject is missing bad components file.');
    end
    
    % Store
    data_continuous_filt_orig = data_continuous_filt;
    data_continuous_filt.orig = data_continuous_filt_orig;
    
    %% RE-FILTER
    
    % Re-filter: grid
    if contains(p.sourcemodel.method, 'grid')
        f = fieldnames(p.preproc);
        f = f(startsWith(f, 'bp'));
        cfg_preproc2 = [];
        for kk = 1:length(f)
            cfg_preproc2.(sprintf('%s', f{kk})) = p.preproc.(sprintf('%s', f{kk}));
        end
        cfg_preproc2.bpfreq = p.sourcemodel.grid.foi;
        data_continuous_filt.grid = ft_preprocessing(cfg_preproc2, data_continuous_filt.orig);
    end
    
    % Refilter: coordinates
    if contains(p.sourcemodel.method, 'coordinates')
        f = fieldnames(p.preproc);
        f = f(startsWith(f, 'bp'));
        cfg_preproc2 = [];
        for kk = 1:length(f)
            cfg_preproc2.(sprintf('%s', f{kk})) = p.preproc.(sprintf('%s', f{kk}));
        end
        cfg_preproc2.bpfreq = p.sourcemodel.coordinates.foi;
        data_continuous_filt.coordinates = ft_preprocessing(cfg_preproc2, data_continuous_filt.orig);
    end
    
    % Refilter: timeseries from grid
    if contains(p.sourcemodel.method, 'grid') && p.sourcemodel.grid.ts.flag
        f = fieldnames(p.preproc);
        f = f(startsWith(f, 'bp'));
        cfg_preproc2 = [];
        for kk = 1:length(f)
            cfg_preproc2.(sprintf('%s', f{kk})) = p.preproc.(sprintf('%s', f{kk}));
        end
        cfg_preproc2.bpfreq = p.sourcemodel.grid.ts.foi;
        data_continuous_filt.grid_ts = ft_preprocessing(cfg_preproc2, data_continuous_filt.orig);
    end
    
    %% EPOCH
    
    % Load in the trial configuration
    cfg_trl = load(p.epoching.trl_cfg(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
    
    % Epoch
    cfg_epoch = [];
    cfg_epoch.trl = cfg_trl.trl;
    data_filt.orig = ft_redefinetrial(cfg_epoch, data_continuous_filt.orig);
    
    % Epoch for each frequency band (grid, coordinates, grid.ts) if
    % necessary
    if contains(p.sourcemodel.method, 'grid')
        data_filt.grid = ft_redefinetrial(cfg_epoch, data_continuous_filt.grid);
    end
    if contains(p.sourcemodel.method, 'coordinates')
        data_filt.coordinates = ft_redefinetrial(cfg_epoch, data_continuous_filt.coordinates);
    end
    if contains(p.sourcemodel.method, 'grid') && p.sourcemodel.grid.ts.flag
        data_filt.grid_ts = ft_redefinetrial(cfg_epoch, data_continuous_filt.grid_ts);
    end
    
    %% HEAD MODEL
    
    % Load in MRI
    mri = ft_read_mri([p.headmodel.template_dir(sub_num) '/' p.headmodel.template_nii(sub_num)]);
    
    % Get the coordinate system
    info = niftiinfo([p.headmodel.template_dir(sub_num) '/' p.headmodel.template_nii(sub_num)]);
    coordsys = '';
    for ii = 1:info.Transform.Dimensionality
        col = info.Transform.T(1:3, ii)/info.PixelDimensions(ii);
        col = round(col, 1);
        ind = find(abs(col) == 1);
        if ind == 1
            if col(ind) == 1
                coordsys = strcat(coordsys, 'r');
            else
                coordsys = strcat(coordsys, 'l');
            end
        elseif ind == 2
            if col(ind) == 1
                coordsys = strcat(coordsys, 'a');
            else
                coordsys = strcat(coordsys, 'p');
            end
        else
            if col(ind) == 1
                coordsys = strcat(coordsys, 's');
            else
                coordsys = strcat(coordsys, 'i');
            end
        end
    end
    mri.coordsys = coordsys;
    
    % Make sure MRI segmentation exists
    if ~exist([p.headmodel.template_dir(sub_num) '/meshes.mat'], 'file')
        
        % Get segmented MRI
        if ~exist([p.headmodel.template_dir(sub_num) '/segmentedmri.mat'], 'file')
            disp('Segmenting MRI...')
            cfg_mri = [];
            cfg_mri.output = {'brain','skull','scalp'};
            cfg_mri.scalpthreshold = 0.05;
            cfg_mri.brainthreshold = 0.5;
            cfg_mri.skullthreshold = 0.5;
            segmentedmri  = ft_volumesegment(cfg_mri, mri);
            parsave([p.headmodel.template_dir(sub_num) '/segmentedmri.mat'], segmentedmri, 'segmentedmri')
            disp('Done!')
        end
        
        % Load
        disp('Loading Segmented MRI...')
        segmentedmri = load([p.headmodel.template_dir(sub_num) '/segmentedmri.mat']);
        disp('Done!')
        
        % Get meshes
        cfg_mesh = [];
        cfg_mesh.tissue = {'brain','skull','scalp'};
        cfg_mesh.numvertices = [3000 2500 2500];
        mesh1 = ft_prepare_mesh(cfg_mesh, segmentedmri.segmentedmri);
        meshes = struct.empty(size(mesh1, 2), 0);
        for nn = 1:size(mesh1,2)
            meshes(nn).pnt = mesh1(nn).pos;
            meshes(nn).tri = mesh1(nn).tri;
            meshes(nn).unit = mesh1(nn).unit;
            meshes(nn).name = cfg_mesh.tissue{nn};
            meshes(nn) = ft_convert_units(meshes(nn),'m');
        end
        parsave_twovar([p.headmodel.template_dir(sub_num) '/meshes.mat'], meshes, 'meshes', segmentedmri.segmentedmri, 'segmentedmri');
        
    end
    
    % Load
    meshes = load([p.headmodel.template_dir(sub_num) '/meshes.mat']);
    
    % Extract
    mri_segmented = meshes.segmentedmri;
    
    % Create the headmodel
    cfg_hm = [];
    cfg_hm.method = p.headmodel.method;
    vol = ft_prepare_headmodel(cfg_hm, mri_segmented);
    
    % Convert to the source position units
    vol = ft_convert_units(vol, p.sourcemodel.units);
    
    % Get the MRI -> CTF transformation
    if ~exist(p.headmodel.tform_inv(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file') || clobber == 1
        
        % Load the sensor transformation matrix
        dat = dlmread(p.raw_data.sens2mri(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
        if size(dat, 2) == 5
            if sum(dat(:, 5)) == 0
                dat = dat(1:4, 1:4);
            end
        end
        
        % Convert to mm
        dat(1:3, end) = dat(1:3, end)*1000;
        
        % Save
        tform = dat;
        writematrix(tform, p.headmodel.tform(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
        
        % Get the inverse of the transform - the inverse is the template to
        % OPM data transformation
        tform_inv = pinv(tform);
        writematrix(tform_inv, p.headmodel.tform_inv(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
        
    else
        
        % Load
        tform_inv = readmatrix(p.headmodel.tform_inv(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
        
    end
    
    %% REGISTRATIONS
    
    % Perform the registrations if necessary
    if ~exist([p.headmodel.template_dir(sub_num) '/registrations/T1_template_registeredMNI152_Warped.nii.gz'], 'file') %|| clobber == 1
        opm_mri_mni_registrations([p.headmodel.template_dir(sub_num) '/' p.headmodel.template_nii(sub_num)], p.directories.anat_dir(sub_num));
    end
    
    %% SOURCE MODEL
    
    % Check if sourcing a grid, or sourcing coordinates
    if contains(p.sourcemodel.method, 'grid')
        
        % Check if using a mask or not
        if ~p.sourcemodel.grid.mask
            
            % Prepare the source model - in template space
            cfg_sm = [];
            cfg_sm.grid.resolution = p.sourcemodel.grid.res;
            cfg_sm.grid.unit = p.sourcemodel.units;
            cfg_sm.grid.tight = 'yes';
            cfg_sm.headmodel = vol;
            sourcemodel_grid = ft_prepare_sourcemodel(cfg_sm);
            
            % Convert to the same units
            sourcemodel_grid = ft_convert_units(sourcemodel_grid, vol.unit);
            
        elseif isfield(p.sourcemodel.grid, 'mask_atlas') && ~strcmp(p.sourcemodel.grid.mask_atlas, 'grey')
            
            % Load in the mask
            mask = ft_read_mri(p.sourcemodel.grid.mask_nii(sub_num));
            mask = ft_convert_units(mask, vol.unit);
            
            % Get the coordinate system
            info = niftiinfo(p.sourcemodel.grid.mask_nii(sub_num));
            coordsys = '';
            for ii = 1:info.Transform.Dimensionality
                col = info.Transform.T(1:3, ii)/info.PixelDimensions(ii);
                col = round(col, 1);
                ind = find(abs(col) == 1);
                if ind == 1
                    if col(ind) == 1
                        coordsys = strcat(coordsys, 'r');
                    else
                        coordsys = strcat(coordsys, 'l');
                    end
                elseif ind == 2
                    if col(ind) == 1
                        coordsys = strcat(coordsys, 'a');
                    else
                        coordsys = strcat(coordsys, 'p');
                    end
                else
                    if col(ind) == 1
                        coordsys = strcat(coordsys, 's');
                    else
                        coordsys = strcat(coordsys, 'i');
                    end
                end
            end
            mask.coordsys = coordsys;
            
            % Convert units for later
            mask_mm = ft_convert_units(mask, 'mm');
            
            % Get the voxel coordinates that are not equal to zero
            [xx, yy, zz] = ind2sub(mask.dim, find(mask.anatomy ~= 0));
            
            % Change to mask units
            pos = [xx, yy, zz]*mask.transform(1:3, 1:3) + mask.transform(1:3, end)';
            clearvars xx yy zz
            
            % Create source model
            cfg_sm = [];
            cfg_sm.sourcemodel.pos = pos;
            cfg_sm.headmodel = vol;
            sourcemodel_grid = ft_prepare_sourcemodel(cfg_sm);
            
            % Clear
            clearvars pos
            
        else
            
            
            % Mask from grey
            
        end
        
        % Check whether we're computing the timeseries of the voxel with max
        % power
        if p.sourcemodel.grid.ts.flag
            
            % Load in the parcellation
            parcel = niftiread(p.sourcemodel.grid.ts.atlas(sub_num));
            parcel_hdr = niftiinfo(p.sourcemodel.grid.ts.atlas(sub_num));
            rois = unique(parcel(parcel ~= 0));
            
        end
        
    end
    
    if contains(p.sourcemodel.method, 'coordinates')
        
        % Prepare the source model - in template space
        cfg_sm = [];
        cfg_sm.grid.pos = readmatrix(p.sourcemodel.coordinates.ve_locs(sub_num));
        cfg_sm.grid.unit = p.sourcemodel.units;
        cfg_sm.headmodel = vol;
        sourcemodel_coord = ft_prepare_sourcemodel(cfg_sm);
        
    end
    
    if ~contains(p.sourcemodel.method, 'grid') && ~contains(p.sourcemodel.method, 'coordinates')
        
        error('Must specify whether sourcing "grid" or "coordinates"');
        
    end
    
    % Grid
    if contains(p.sourcemodel.method, 'grid')
        
        % Now warp to CTF space
        sourcemodel_grid_realigned = ft_transform_geometry(tform_inv, sourcemodel_grid);
        vol_realigned = ft_transform_geometry(tform_inv, vol);
        grad = ft_convert_units(data_filt.orig.grad, p.sourcemodel.units);
        
        % Check alignment of source model and head model and save as image
        if ~exist(p.sourcemodel.sourcemodel_fig(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), 'grid'), 'file') || clobber == 1
            figure;
            hold on;
            ft_plot_vol(vol_realigned, 'EdgeColor', [0.7, 0.7, 0.7], 'EdgeAlpha', 0.1, 'FaceAlpha', 0.1);
            ft_plot_mesh(sourcemodel_grid_realigned.pos(sourcemodel_grid_realigned.inside, :), 'VertexColor', 'r', 'VertexSize', 20); %plot source model
            ft_plot_sens(grad,'style','b*'); %plot MEG channels (sensors)
            hold off;
            view(45,10);
            hf = gcf;
            hf.Position(1:2) = [10 10];
            hf.Position(3:4) = (800 / hf.Position(4)) .* hf.Position(3:4);
            saveas(hf, p.sourcemodel.sourcemodel_fig(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), 'grid'));
        end
        
    end
    
    % Coordinates
    if contains(p.sourcemodel.method, 'coordinates')
        
        % Now warp to CTF space
        sourcemodel_coord_realigned = ft_transform_geometry(tform_inv, sourcemodel_coord);
        vol_realigned = ft_transform_geometry(tform_inv, vol);
        grad = ft_convert_units(data_filt.orig.grad, p.sourcemodel.units);
        
        % Check alignment of source model and head model and save as image
        if ~exist(p.sourcemodel.sourcemodel_fig(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), 'coordinates'), 'file') || clobber == 1
            figure;
            hold on;
            ft_plot_vol(vol_realigned, 'EdgeColor', [0.7, 0.7, 0.7], 'EdgeAlpha', 0.9, 'FaceAlpha', 0.1);
            ft_plot_mesh(sourcemodel_coord_realigned.pos(sourcemodel_coord_realigned.inside, :), 'VertexColor', 'r', 'VertexSize', 20); %plot source model
            ft_plot_sens(grad,'style','b*'); %plot MEG channels (sensors)
            hold off;
            view(45,10);
            hf = gcf;
            hf.Position(1:2) = [10 10];
            hf.Position(3:4) = (800 / hf.Position(4)) .* hf.Position(3:4);
            saveas(hf, p.sourcemodel.sourcemodel_fig(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), 'coordinates'));
        end
        
    end
    
    %% COVARIANCE
    
    % Grid
    if contains(p.sourcemodel.method, 'grid')
        
        % Make the epoched data (without bad trials) into continuous
        data_continuous_filt_bf = [data_filt.grid.trial{:}];
        
        % Compute the covariance matrix across the continuous data
        C_store = cov(data_continuous_filt_bf');
        
        % Regularization
        maxEV = max(svd(C_store));
        C = C_store + p.sourcemodel.mu.*maxEV.*eye(size(C_store));
        Cinv.grid = inv(C);
        
    end
    
    % Coordinates
    if contains(p.sourcemodel.method, 'coordinates')
        
        % Make the epoched data (without bad trials) into continuous
        data_continuous_filt_bf = [data_filt.coordinates.trial{:}];
        
        % Compute the covariance matrix across the continuous data
        C_store = cov(data_continuous_filt_bf');
        
        % Regularization
        maxEV = max(svd(C_store));
        C = C_store + p.sourcemodel.mu.*maxEV.*eye(size(C_store));
        Cinv.coordinates = inv(C);
        
    end
    
    % Grid timeseries
    if contains(p.sourcemodel.method, 'grid') && p.sourcemodel.grid.ts.flag
        
        % Make the epoched data (without bad trials) into continuous
        data_continuous_filt_bf = [data_filt.grid_ts.trial{:}];
        
        % Compute the covariance matrix across the continuous data
        C_store = cov(data_continuous_filt_bf');
        
        % Regularization
        maxEV = max(svd(C_store));
        C = C_store + p.sourcemodel.mu.*maxEV.*eye(size(C_store));
        Cinv.grid_ts = inv(C);
        
    end
    
    %% BEAMFORMING
    
    % First, compute the leadfield, and then beamform
    
    % Grid
    if contains(p.sourcemodel.method, 'grid')
        
        if ~exist(p.sourcemodel.grid.weights(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file') || clobber == 1
            
            % Create the leadfield - note, internally, everything is set to m,
            % and then converted to nm
            LF_grid = opm_create_leadfield(grad, vol_realigned, sourcemodel_grid_realigned, M);
            
            % Now, run the beamforming to get the weights
            W_grid = beamform(sourcemodel_grid_realigned, vol_realigned, LF_grid, Cinv.grid, p.sourcemodel.orientation, p.sourcemodel.nai);
            
            % Also run the beamforming to get the weights for the
            % timeseries (if the foi is different)
            if p.sourcemodel.grid.ts.flag
                
                % Check if foi is different
                if p.sourcemodel.grid.ts.foi ~= p.sourcemodel.grid.foi
                    W_grid_ts = beamform(sourcemodel_grid_realigned, vol_realigned, LF_grid, Cinv.grid_ts, p.sourcemodel.orientation, p.sourcemodel.nai);
                else
                    W_grid_ts = W_grid;
                end
                
                % Save
                parsave(p.sourcemodel.grid.ts.weights(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), W_grid_ts, 'W_grid_ts');
                
            end
            
            % Save
            parsave(p.sourcemodel.grid.weights(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), W_grid, 'W_grid');
            
        else
            
            % Load grid beamformer weights
            W_grid = load(p.sourcemodel.grid.weights(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
            W_grid = W_grid.W_grid;
            
            % Load the grid beamformer weights for the timeseries
            if p.sourcemodel.grid.ts.flag
                W_grid_ts = load(p.sourcemodel.grid.ts.weights(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
                W_grid_ts = W_grid_ts.W_grid_ts;
            end
            
        end
        
    end
    
    % Coordinates
    if contains(p.sourcemodel.method, 'coordinates')
        
        if ~exist(p.sourcemodel.coordinates.weights(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file') || clobber == 1
            
            % Create the leadfield - note, internally, everything is set to m,
            % and then converted to nm
            LF_coord = opm_create_leadfield(grad, vol_realigned, sourcemodel_coord_realigned, M);
            
            % Now, run the beamforming to get the weights
            W_coord = beamform(sourcemodel_coord_realigned, vol_realigned, LF_coord, Cinv.coordinates, p.sourcemodel.orientation, p.sourcemodel.nai);
            
            % Save
            parsave(p.sourcemodel.coordinates.weights(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), W_coord, 'W_coord');
            
        else
            
            % Load
            W_coord = load(p.sourcemodel.coordinates.weights(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)));
            W_coord = W_coord.W_coord;
            
        end
        
    end
    
    %% SOURCE MODEL - GRID POWER
    
    % Only run if specified
    if contains(p.sourcemodel.method, 'grid') && p.sourcemodel.grid.avgpow.flag
        
        % All trials
        
        % Iterate over timewindows
        for tt = 1:size(p.sourcemodel.grid.avgpow.twoi, 1)
            
            if ~exist([p.sourcemodel.grid.avgpow.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.avgpow.twoi(tt, 1), p.sourcemodel.grid.avgpow.twoi(tt, 2), 'all') '.nii.gz'], 'file') || clobber == 1
                
                % First, do the baseline
                
                % Extract the time window
                twoi = [-abs(diff(p.sourcemodel.grid.avgpow.twoi(tt, :))), 0];
                
                % Convert to samples
                twoi_samp = twoi*hdr.Fs - cfg_trl.trl(1, 3) + 1;
                
                % Extract the data
                dat = zeros(length(data_filt.grid.label), abs(diff(twoi_samp)) + 1, length(data_filt.grid.trial));
                for ii = 1:length(data_filt.grid.trial)
                    dat(:, :, ii) = data_filt.grid.trial{ii}(:, twoi_samp(1):twoi_samp(2));
                end
                dat = reshape(dat, length(data_filt.grid.label), []);
                
                % Take the covariance
                C_twoi = cov(dat');
                
                % Iterate over the grid, and compute the power for each
                % location using the beamformer weights - do this for each
                % of the tangent vectors, and take the norm
                Q_baseline = zeros(size(W_grid, 2), 1);
                for ii = 1:size(W_grid, 2)
                    Q_temp = zeros(size(W_grid, 3), 1);
                    for jj = 1:size(W_grid, 3)
                        Q_temp(jj) = W_grid(:, ii, jj)'*C_twoi*W_grid(:, ii, jj);
                    end
                    if size(W_grid, 3) > 1
                        Q_baseline(ii) = sqrt(sum(Q_temp(:).^2));
                    else
                        Q_baseline(ii) = Q_temp;
                    end
                end
                Q_baseline(isnan(Q_baseline)) = 0;
                
                % Now: time window
                
                % Extract the time window
                twoi = p.sourcemodel.grid.avgpow.twoi(tt, :);
                
                % Convert to samples
                twoi_samp = twoi*hdr.Fs - cfg_trl.trl(1, 3) + 1;
                
                % Make sure we're within the trial limits (by 1)
                if twoi_samp(2) == size(data_filt.grid.trial{1}, 2)+1
                    twoi_samp(2) = twoi_samp(2)-1;
                end
                
                % Extract the data
                dat = zeros(length(data_filt.grid.label), abs(diff(twoi_samp)) + 1, length(data_filt.grid.trial));
                for ii = 1:length(data_filt.grid.trial)
                    dat(:, :, ii) = data_filt.grid.trial{ii}(:, twoi_samp(1):twoi_samp(2));
                end
                dat = reshape(dat, length(data_filt.grid.label), []);
                
                % Take the covariance
                C_twoi = cov(dat');
                
                % Iterate over the grid, and compute the power for each
                % location using the beamformer weights
                Q_active = zeros(size(W_grid, 2), 1);
                for ii = 1:size(W_grid, 2)
                    Q_temp = zeros(size(W_grid, 3), 1);
                    for jj = 1:size(W_grid, 3)
                        Q_temp(jj) = W_grid(:, ii, jj)'*C_twoi*W_grid(:, ii, jj);
                    end
                    if size(W_grid, 3) > 1
                        Q_active(ii) = sqrt(sum(Q_temp(:).^2));
                    else
                        Q_active(ii) = Q_temp;
                    end
                end
                Q_active(isnan(Q_active)) = 0;
                
                % Convert to percentage change
                Q = (Q_active - Q_baseline)./Q_baseline;
                
                % Change back into mask if specified
                if p.sourcemodel.grid.mask
                    
                    % Assign back to mask voxels
                    mask_pow = mask_mm;
                    t = pinv(mask_mm.transform);
                    for ii = 1:length(Q)
                        pos = sourcemodel_grid.pos(ii, :);
                        pos_vox = round(pos*t(1:3, 1:3) + t(1:3, 4)');
                        mask_pow.anatomy(pos_vox(1), pos_vox(2), pos_vox(3)) = Q(ii);
                    end
                    
                    % Write
                    ft_write_mri([p.sourcemodel.grid.avgpow.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.avgpow.twoi(tt, 1), p.sourcemodel.grid.avgpow.twoi(tt, 2), 'all') '.nii.gz'], mask_pow, 'dataformat', 'nifti_gz');
                    
                else
                    
                    % Put in FieldTrip structure in order to write nifti
                    pow = [];
                    pow.Q = Q;
                    pow.pos = sourcemodel_grid.pos;
                    pow.unit = sourcemodel_grid.unit;
                    pow.dim = sourcemodel_grid.dim;
                    
                    % Write
                    cfg_wn = [];
                    cfg_wn.filename =  p.sourcemodel.grid.avgpow.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.avgpow.twoi(tt, 1), p.sourcemodel.grid.avgpow.twoi(tt, 2), 'all');
                    cfg_wn.parameter = 'Q';
                    cfg_wn.precision = 'double';
                    cfg_wn.filetype = 'nifti';
                    ft_sourcewrite(cfg_wn, pow);
                    cmd = ['gzip ' p.sourcemodel.grid.avgpow.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.avgpow.twoi(tt, 1), p.sourcemodel.grid.avgpow.twoi(tt, 2), 'all') '.nii'];
                    system(cmd);
                    
                end
                
            end
            
            % Register to MNI space if specified
            if p.sourcemodel.grid.reg_mni && (~exist([p.sourcemodel.grid.avgpow.nii_mni(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.avgpow.twoi(tt, 1), p.sourcemodel.grid.avgpow.twoi(tt, 2), 'all') '.nii.gz'], 'file') || clobber == 1)
                
                % Register to MNI space
                cmd = ['antsApplyTransforms -d 3 -i ' p.sourcemodel.grid.avgpow.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.avgpow.twoi(tt, 1), p.sourcemodel.grid.avgpow.twoi(tt, 2), 'all') '.nii.gz -r ' p.headmodel.template_dir(sub_num) '/registrations/T1_template_registeredMNI152_Warped.nii.gz -t ' p.headmodel.template_dir(sub_num) '/registrations/T1_template_registeredMNI152_1Warp.nii.gz -t ' p.headmodel.template_dir(sub_num) '/registrations/T1_template_registeredMNI152_0GenericAffine.mat -o ' p.sourcemodel.grid.avgpow.nii_mni(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.avgpow.twoi(tt, 1), p.sourcemodel.grid.avgpow.twoi(tt, 2), 'all') '.nii.gz'];
                system(cmd);
                
            end
            
        end
        
        % Determine whether we're doing this per condition as well
        if p.sourcemodel.grid.avgpow.per_cond
            
            % Per condition
            
            % Iterate over the conditions
            cond_ind = unique(cfg_trl.trl(:, 4));
            for cc = 1:length(cond_ind)
                
                % Iterate over timewindows
                for tt = 1:size(p.sourcemodel.grid.avgpow.twoi, 1)
                    
                    if ~exist([p.sourcemodel.grid.avgpow.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.avgpow.twoi(tt, 1), p.sourcemodel.grid.avgpow.twoi(tt, 2), trial_info.details(cc).name) '.nii.gz'], 'file') || clobber == 1
                        
                        % First, do the baseline
                        
                        % Extract the time window
                        twoi = [-abs(diff(p.sourcemodel.grid.avgpow.twoi(tt, :))), 0];
                        
                        % Convert to samples
                        twoi_samp = twoi*hdr.Fs - cfg_trl.trl(1, 3) + 1;
                        
                        % Get the data for the current condition
                        data_filt_cond = data_filt.grid.trial(cfg_trl.trl(:, 4) == cond_ind(cc));
                        
                        % Extract the data
                        dat = zeros(length(data_filt.grid.label), abs(diff(twoi_samp)) + 1, length(data_filt_cond));
                        for ii = 1:length(data_filt_cond)
                            dat(:, :, ii) = data_filt_cond{ii}(:, twoi_samp(1):twoi_samp(2));
                        end
                        dat = reshape(dat, length(data_filt.grid.label), []);
                        
                        % Take the covariance
                        C_twoi = cov(dat');
                        
                        % Iterate over the grid, and compute the power for each
                        % location using the beamformer weights
                        Q_baseline = zeros(size(W_grid, 2), 1);
                        for ii = 1:size(W_grid, 2)
                            Q_temp = zeros(size(W_grid, 3), 1);
                            for jj = 1:size(W_grid, 3)
                                Q_temp(jj) = W_grid(:, ii, jj)'*C_twoi*W_grid(:, ii, jj);
                            end
                            if size(W_grid, 3) > 1
                                Q_baseline(ii) = sqrt(sum(Q_temp(:).^2));
                            else
                                Q_baseline(ii) = Q_temp;
                            end
                        end
                        Q_baseline(isnan(Q_baseline)) = 0;
                        
                        % Now: time window
                        
                        % Extract the time window
                        twoi = p.sourcemodel.grid.avgpow.twoi(tt, :);
                        
                        % Convert to samples
                        twoi_samp = twoi*hdr.Fs - cfg_trl.trl(1, 3) + 1;
                        
                        % Make sure we're within the trial limits (by 1)
                        if twoi_samp(2) == size(data_filt.grid.trial{1}, 2)+1
                            twoi_samp(2) = twoi_samp(2)-1;
                        end
                        
                        % Get the data for the current condition
                        data_filt_cond = data_filt.grid.trial(cfg_trl.trl(:, 4) == cond_ind(cc));
                        
                        % Extract the data
                        dat = zeros(length(data_filt.grid.label), abs(diff(twoi_samp)) + 1, length(data_filt_cond));
                        for ii = 1:length(data_filt_cond)
                            dat(:, :, ii) = data_filt_cond{ii}(:, twoi_samp(1):twoi_samp(2));
                        end
                        dat = reshape(dat, length(data_filt.grid.label), []);
                        
                        % Take the covariance
                        C_twoi = cov(dat');
                        
                        % Iterate over the grid, and compute the power for each
                        % location using the beamformer weights
                        Q_active = zeros(size(W_grid, 2), 1);
                        for ii = 1:size(W_grid, 2)
                            Q_temp = zeros(size(W_grid, 3), 1);
                            for jj = 1:size(W_grid, 3)
                                Q_temp(jj) = W_grid(:, ii, jj)'*C_twoi*W_grid(:, ii, jj);
                            end
                            if size(W_grid, 3) > 1
                                Q_active(ii) = sqrt(sum(Q_temp(:).^2));
                            else
                                Q_active(ii) = Q_temp;
                            end
                        end
                        Q_active(isnan(Q_active)) = 0;
                        
                        % Convert to percentage change
                        Q = (Q_active - Q_baseline)./Q_baseline;
                        
                        % Change back into mask if specified
                        if p.sourcemodel.grid.mask
                            
                            % Assign back to mask voxels
                            mask_pow = mask_mm;
                            t = pinv(mask_mm.transform);
                            for ii = 1:length(Q)
                                pos = sourcemodel_grid.pos(ii, :);
                                pos_vox = round(pos*t(1:3, 1:3) + t(1:3, 4)');
                                mask_pow.anatomy(pos_vox(1), pos_vox(2), pos_vox(3)) = Q(ii);
                            end
                            
                            % Write
                            ft_write_mri([p.sourcemodel.grid.avgpow.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.avgpow.twoi(tt, 1), p.sourcemodel.grid.avgpow.twoi(tt, 2), trial_info.details(cc).name) '.nii.gz'], mask_pow, 'dataformat', 'nifti_gz');
                            
                        else
                            
                            % Put in FieldTrip structure in order to write nifti
                            pow = [];
                            pow.Q = Q;
                            pow.pos = sourcemodel_grid.pos;
                            pow.unit = sourcemodel_grid.unit;
                            pow.dim = sourcemodel_grid.dim;
                            
                            % Write
                            cfg_wn = [];
                            cfg_wn.filename =  p.sourcemodel.grid.avgpow.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.avgpow.twoi(tt, 1), p.sourcemodel.grid.avgpow.twoi(tt, 2), trial_info.details(cc).name);
                            cfg_wn.parameter = 'Q';
                            cfg_wn.precision = 'double';
                            cfg_wn.filetype = 'nifti';
                            ft_sourcewrite(cfg_wn, pow);
                            cmd = ['gzip ' p.sourcemodel.grid.avgpow.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.avgpow.twoi(tt, 1), p.sourcemodel.grid.avgpow.twoi(tt, 2), trial_info.details(cc).name) '.nii'];
                            system(cmd);
                            
                        end
                        
                    end
                    
                    % Register to MNI space if specified
                    if p.sourcemodel.grid.reg_mni && (~exist([p.sourcemodel.grid.avgpow.nii_mni(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.avgpow.twoi(tt, 1), p.sourcemodel.grid.avgpow.twoi(tt, 2), trial_info.details(cc).name) '.nii.gz'], 'file') || clobber == 1)
                        
                        % Register to MNI space
                        cmd = ['antsApplyTransforms -d 3 -i ' p.sourcemodel.grid.avgpow.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.avgpow.twoi(tt, 1), p.sourcemodel.grid.avgpow.twoi(tt, 2), trial_info.details(cc).name) '.nii.gz -r ' p.headmodel.template_dir(sub_num) '/registrations/T1_template_registeredMNI152_Warped.nii.gz -t ' p.headmodel.template_dir(sub_num) '/registrations/T1_template_registeredMNI152_1Warp.nii.gz -t ' p.headmodel.template_dir(sub_num) '/registrations/T1_template_registeredMNI152_0GenericAffine.mat -o ' p.sourcemodel.grid.avgpow.nii_mni(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.avgpow.twoi(tt, 1), p.sourcemodel.grid.avgpow.twoi(tt, 2), trial_info.details(cc).name) '.nii.gz'];
                        system(cmd);
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
    %% SOURCEMODEL - GRID TIMESERIES
    
    % Only run if specified
    if contains(p.sourcemodel.method, 'grid') && p.sourcemodel.grid.ts.flag
        
        if ~exist(p.sourcemodel.grid.ts.ve(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file') || clobber == 1
            
            % Compute the power over the whole baseline period
            
            % Get the time window
            twoi = [-abs(diff( p.sourcemodel.grid.ts.twoi)), 0];
            
            % Convert to samples
            twoi_samp = twoi*hdr.Fs - cfg_trl.trl(1, 3) + 1;
            
            % Extract the data
            dat = zeros(length(data_filt.grid.label), abs(diff(twoi_samp)) + 1, length(data_filt.grid.trial));
            for ii = 1:length(data_filt.grid.trial)
                dat(:, :, ii) = data_filt.grid.trial{ii}(:, twoi_samp(1):twoi_samp(2));
            end
            dat = reshape(dat, length(data_filt.grid.label), []);
            
            % Take the covariance
            C_twoi = cov(dat');
            
            % Iterate over the grid, and compute the power for each
            % location using the beamformer weights
            Q_baseline = zeros(size(W_grid, 2), 1);
            for ii = 1:size(W_grid, 2)
                Q_temp = zeros(size(W_grid, 3), 1);
                for jj = 1:size(W_grid, 3)
                    Q_temp(jj) = W_grid(:, ii, jj)'*C_twoi*W_grid(:, ii, jj);
                end
                if size(W_grid, 3) > 1
                    Q_baseline(ii) = sqrt(sum(Q_temp(:).^2));
                else
                    Q_baseline(ii) = Q_temp;
                end
            end
            Q_baseline(isnan(Q_baseline)) = 0;
            
            % Compute the power over the whole active period
            
            % Get the time window
            twoi = p.sourcemodel.grid.ts.twoi;
            
            % Convert to samples
            twoi_samp = twoi*hdr.Fs - cfg_trl.trl(1, 3) + 1;
            
            % Make sure we're within the trial limits (by 1)
            if twoi_samp(2) == size(data_filt.grid.trial{1}, 2)+1
                twoi_samp(2) = twoi_samp(2)-1;
            end
            
            % Extract the data
            dat = zeros(length(data_filt.grid.label), abs(diff(twoi_samp)) + 1, length(data_filt.grid.trial));
            for ii = 1:length(data_filt.grid.trial)
                dat(:, :, ii) = data_filt.grid.trial{ii}(:, twoi_samp(1):twoi_samp(2));
            end
            dat = reshape(dat, length(data_filt.grid.label), []);
            
            % Take the covariance
            C_twoi = cov(dat');
            
            % Iterate over the grid, and compute the power for each
            % location using the beamformer weights
            Q_active = zeros(size(W_grid, 2), 1);
            for ii = 1:size(W_grid, 2)
                Q_temp = zeros(size(W_grid, 3), 1);
                for jj = 1:size(W_grid, 3)
                    Q_temp(jj) = W_grid(:, ii, jj)'*C_twoi*W_grid(:, ii, jj);
                end
                if size(W_grid, 3) > 1
                    Q_active(ii) = sqrt(sum(Q_temp(:).^2));
                else
                    Q_active(ii) = Q_temp;
                end
            end
            Q_active(isnan(Q_active)) = 0;
            
            % Convert to percentage change
            Q = (Q_active - Q_baseline)./Q_baseline;
            
            % Change back into mask if specified
            if p.sourcemodel.grid.mask
                
                % Assign back to mask voxels
                Q_reshape = zeros(size(mask_mm.anatomy));
                W_grid_ts_cell = cell(size(W_grid_ts, 3), 1);
                t = pinv(mask_mm.transform);
                for jj = 1:size(W_grid_ts, 3)
                    W_reshape = zeros(size(W_grid_ts, 1), size(mask_mm.anatomy, 1), size(mask_mm.anatomy, 2), size(mask_mm.anatomy, 3));
                    for ii = 1:length(Q)
                        pos = sourcemodel_grid.pos(ii, :);
                        pos_vox = round(pos*t(1:3, 1:3) + t(1:3, 4)');
                        W_reshape(:, pos_vox(1), pos_vox(2), pos_vox(3)) = W_grid_ts(:, ii, jj);
                    end
                    W_grid_ts_cell{jj} = reshape(W_reshape, size(W_grid_ts, 1), []);
                end
                for ii = 1:length(Q)
                    pos = sourcemodel_grid.pos(ii, :);
                    pos_vox = round(pos*t(1:3, 1:3) + t(1:3, 4)');
                    Q_reshape(pos_vox(1), pos_vox(2), pos_vox(3)) = Q(ii);
                end
                Q = reshape(Q_reshape, [], 1);
                
            else
                
                W_grid_ts_cell = cell(size(W_grid_ts, 3), 1);
                for jj = 1:size(W_grid_ts, 3)
                    W_grid_ts_cell{jj} = W_grid_ts(:, :, jj);
                end
                
            end
            
            % Load in parcel names
            temp = readtable(p.sourcemodel.grid.ts.atlas_regions(sub_num), 'ReadVariableNames', false);
            
            % Load in the parcellation
            parcel = niftiread(p.sourcemodel.grid.ts.atlas(sub_num));
            parcel_hdr = niftiinfo(p.sourcemodel.grid.ts.atlas(sub_num));
            if p.sourcemodel.grid.ts.all
                if p.sourcemodel.grid.ts.all_rois
                parcel(~((parcel >= p.sourcemodel.grid.ts.all_rois(1)) & (parcel <= p.sourcemodel.grid.ts.all_rois(2)))) = 0;
                else
                end
                rois = parcel(parcel ~= 0);
            else

                rois = unique(parcel(parcel ~= 0));
            end 
            % Create the save structure
            VE = [];
            VE.ts = zeros(length(rois), length(data_filt.grid_ts.time{1}), length(data_filt.grid_ts.trial));
            VE.pos_coord = zeros(length(rois), 3);
            VE.pos_name = temp.Var2;
            VE.srate = data_filt.grid_ts.fsample;
            VE.time = data_filt.grid_ts.time{1};
            inds = find(parcel);
            % Iterate over the regions and compute the timeseries - sources x time x trials
            for rr = 1:length(rois)
                if p.sourcemodel.grid.ts.all
                    
                    q_val = Q(inds(rr));
                    [xx,yy,zz] = ind2sub(size(parcel),inds(rr));
                    coord = [xx,yy,zz]-1;
                    % Convert to mm
                    coord_mm = parcel_hdr.Transform.T(1:3, 1:3)*coord' + parcel_hdr.Transform.T(4, 1:3)';
                    VE.pos_coord(rr, :) = coord_mm;

                    % Compute the timeseries for each trial
                    for tt = 1:length(data_filt.grid_ts.trial)
                        ts = zeros(length(data_filt.grid_ts.time{1}), length(W_grid_ts_cell));
                        for jj = 1:length(W_grid_ts_cell)
                            ts(:, jj) = W_grid_ts_cell{jj}(:, inds(rr))'*data_filt.grid_ts.trial{tt};
                        end
                        if length(W_grid_ts_cell) > 1
                            VE.ts(rr, :, tt) = sqrt(sum(ts.^2, 2));
                        else
                            VE.ts(rr, :, tt) = ts;
                        end
                    end

                else

                    % Get the indices of this parcel
                    ind = find(parcel == rois(rr));

                    % Get the Q values
                    q_val = Q(ind);

                    % Get the maximum
                    [~, max_ind] = max(q_val);
                    if isfield(p.sourcemodel,'min_flag')
                        [~, max_ind] = min(q_val);
                    end
                    % Get the coordinates in xyz
                    [xx, yy, zz] = ind2sub(size(parcel), ind(max_ind));
                    coord = [xx, yy, zz]-1;

                    % Convert to mm
                    coord_mm = parcel_hdr.Transform.T(1:3, 1:3)*coord' + parcel_hdr.Transform.T(4, 1:3)';
                    VE.pos_coord(rr, :) = coord_mm;

                    % Compute the timeseries for each trial
                    for tt = 1:length(data_filt.grid_ts.trial)
                        ts = zeros(length(data_filt.grid_ts.time{1}), length(W_grid_ts_cell));
                        for jj = 1:length(W_grid_ts_cell)
                            ts(:, jj) = W_grid_ts_cell{jj}(:, ind(max_ind))'*data_filt.grid_ts.trial{tt};
                        end
                        if length(W_grid_ts_cell) > 1
                            VE.ts(rr, :, tt) = sqrt(sum(ts.^2, 2));
                        else
                            VE.ts(rr, :, tt) = ts;
                        end
                    end
                end
            end

            % Flip x and y to convert nifti to ITK
            VE.pos_coord(:, 1) = -VE.pos_coord(:, 1);
            VE.pos_coord(:, 2) = -VE.pos_coord(:, 2);

            % Transform to MNI
            writetable(array2table(VE.pos_coord, 'VariableNames', {'x', 'y', 'z'}), [p.directories.sub_dir(p.subject(sub_num), p.session(sub_num)) '/temp.csv']);
            cmd = ['antsApplyTransformsToPoints -d 3 -i ' p.directories.sub_dir(p.subject(sub_num), p.session(sub_num)) '/temp.csv -o ' p.directories.sub_dir(p.subject(sub_num), p.session(sub_num)) '/temp2.csv -t [' p.headmodel.template_dir(sub_num) '/registrations/T1_template_registeredMNI152_0GenericAffine.mat, 1] -t ' p.headmodel.template_dir(sub_num) '/registrations/T1_template_registeredMNI152_1InverseWarp.nii.gz'];
            system(cmd);
            tab = readtable([p.directories.sub_dir(p.subject(sub_num), p.session(sub_num)) '/temp2.csv']);

            % Flip back (ITK -> nifti)
            tab.x = -tab.x;
            tab.y = -tab.y;

            % Assign
            VE.pos_coord_MNI = table2array(tab);
            delete([p.directories.sub_dir(p.subject(sub_num), p.session(sub_num)) '/temp.csv']);
            
            % Save
            parsave_struct(p.sourcemodel.grid.ts.ve(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), VE, 'VE');
            
        end
        
    end
    
    %% SOURCEMODEL - EXTRACT TIMESERIES FOR EACH VOXEL
    
    if contains(p.sourcemodel.method, 'grid') && (p.sourcemodel.grid.seed.flag || p.sourcemodel.grid.psd.flag)
        
        % Check if we're missing anything
        ex = [false, false];
        if p.sourcemodel.grid.seed.flag
            if exist(p.sourcemodel.grid.seed.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.seed.names{end}), 'file') && clobber == 0
                ex(1) = true;
            end
        else
            ex(1) = true;
        end
        if p.sourcemodel.grid.psd.flag
            if exist(p.sourcemodel.grid.psd.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.psd.foi(end, 1), p.sourcemodel.grid.psd.foi(end, 2)), 'file') && clobber == 0
                ex(2) = true;
            end
        else
            ex(2) = true;
        end
        
        % Only proceed if last one doesn't exist
        if ~all(ex)
            
            % Construct filter coefficients if necessary
            if p.sourcemodel.grid.seed.flag && p.sourcemodel.grid.seed.filter
                
                % Current frequency
                filtband = p.sourcemodel.grid.seed.foi;
                
                % Estimate the filter order
                freqs = [filtband(1) * 0.6, filtband(1), filtband(2), filtband(2) * 1.9];
                fir_n = kaiserord(freqs, [0 1 0], [0.25 0.05 0.25], hdr.Fs);
                
                % Make the filter coefficients
                fir_coef = fir1(fir_n, filtband*2/hdr.Fs, 'bandpass');
                
            end
            
            % Compute the timeseries
            inside = find(sourcemodel_grid.inside);
            bldata_store = zeros(length(data_continuous_filt.grid.time{1}), length(inside));
            bldata_filt_store = zeros(length(data_continuous_filt.grid.time{1}), length(inside));
            for ii = 1:length(inside)
                
                % Now, compute the timeseries - time x trials
                ts_temp = zeros(length(data_continuous_filt.grid.time{1}), size(W_grid, 3));
                for jj = 1:size(W_grid, 3)
                    ts_temp(:, jj) = W_grid(:, inside(ii), jj)'*data_continuous_filt.grid.trial{1};
                end
                if size(W_grid, 3) > 1
                    ts = sqrt(sum(ts_temp.^2, 3));
                else
                    ts = ts_temp;
                end
                
                % Mean center
                bldata_store(:, ii) = ts' - mean(ts);
                
                % Filter if necessary
                if p.sourcemodel.grid.seed.filter
                    bldata_filt_store(:, ii) = filtfilt(fir_coef, 1, bldata_store(:, ii));
                end
                
            end
            
        end
        
    end
    
    
    %% SOURCEMODEL - GRID SEED
    
    % Only run if specified
    if contains(p.sourcemodel.method, 'grid') && p.sourcemodel.grid.seed.flag
        
        % Load in the coordinates
        coords = readmatrix(p.sourcemodel.grid.seed.coords(sub_num));
        coords = coords(p.sourcemodel.grid.seed.coord_ids, :);
        
        % Only proceed if last one doesn't exist
        if ~exist(p.sourcemodel.grid.seed.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.seed.names{end}), 'file') || clobber == 1
            
            % Use filtered or not
            if p.sourcemodel.grid.seed.filter
                bldata = bldata_filt_store;
            else
                bldata = bldata_store;
            end
            
            % Iterate over the seeds
            for ss = 1:size(coords, 1)
                
                if ~exist(p.sourcemodel.grid.seed.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.seed.names{ss}), 'file') || clobber == 1
                    
                    % Find closest index
                    dist = sqrt((sourcemodel_grid.pos(:, 1) - coords(ss, 1)).^2 + (sourcemodel_grid.pos(:, 2) - coords(ss, 2)).^2 + (sourcemodel_grid.pos(:, 3) - coords(ss, 3)).^2);
                    dist(~sourcemodel_grid.inside) = [];
                    [~, ind] = min(dist);
                    
                    % Compute the orthogonalization factor
                    fac = (conj(hilbert(bldata(:, ind)))./abs(hilbert(bldata(:, ind))));
                    
                    % Take the absolute value of the analytic signal
                    H_data_o_seed = abs(hilbert(bldata(:, ind)));
                    
                    % Brookes does an extra de-mean before correlating - will also perform here
                    H_data_o_seed = H_data_o_seed - mean(H_data_o_seed);
                    
                    % Redefine into trials (only those we're keeping)
                    H_data_o_t_seed = zeros(length(data_filt.grid.time{1}), length(data_filt.grid.trial));
                    for tt = 1:length(data_filt.grid.trial)
                        H_data_o_t_seed(:, tt) = H_data_o_seed(cfg_trl.trl(tt, 1):cfg_trl.trl(tt, 2));
                    end
                    
                    % Downsample if specified
                    if p.sourcemodel.grid.seed.downsample.flag
                        timet = timetable(H_data_o_t_seed, 'SampleRate', data_filt.orig.fsample);
                        tt_avg = retime(timet, 'regular', 'mean', 'SampleRate', p.sourcemodel.grid.seed.downsample.freq);
                        H_data_o_t_seed = [tt_avg.H_data_o_t_seed];
                    end
                    
                    % Reshape
                    H_data_o_seed = reshape(H_data_o_t_seed, [], 1);
                    
                    % Iterate over the voxels: slow, but otherwise there
                    % are memory issues
                    R = zeros(size(bldata, 2), 1);
                    for ii = 1:size(bldata, 2)
                        
                        % Orthogonalize
                        H_data_o = abs(imag(hilbert(bldata(:, ii)).*fac));
                        
                        % Brookes does an extra de-mean before correlating - will also perform here
                        H_data_o = H_data_o - mean(H_data_o);
                        
                        % Redefine into trials (only those we're keeping)
                        H_data_o_t = zeros(length(data_filt.grid.time{1}), length(data_filt.grid.trial));
                        for tt = 1:length(data_filt.grid.trial)
                            H_data_o_t(:, tt) = H_data_o(cfg_trl.trl(tt, 1):cfg_trl.trl(tt, 2));
                        end
                        
                        % Downsample if specified
                        if p.sourcemodel.grid.seed.downsample.flag
                            timet = timetable(H_data_o_t, 'SampleRate', data_filt.orig.fsample);
                            tt_avg = retime(timet, 'regular', 'mean', 'SampleRate', 5);
                            H_data_o_t = [tt_avg.H_data_o_t];
                        end
                        
                        % Reshape
                        H_data_o = reshape(H_data_o_t, [], 1);
                        
                        % Correlate
                        R(ii) = corr(H_data_o_seed, H_data_o);
                        
                    end
                    
                    % Assign back the inside and average across trials
                    R_avg = zeros(size(W_grid, 2), 1);
                    R_avg(sourcemodel_grid.inside) = mean(R, 2);
                    
                    % Change back into mask if specified
                    if p.sourcemodel.grid.mask
                        
                        % Assign back to mask voxels
                        mask_pow = mask_mm;
                        t = pinv(mask_mm.transform);
                        for ii = 1:length(R_avg)
                            pos = sourcemodel_grid.pos(ii, :);
                            pos_vox = round(pos*t(1:3, 1:3) + t(1:3, 4)');
                            mask_pow.anatomy(pos_vox(1), pos_vox(2), pos_vox(3)) = R_avg(ii);
                        end
                        
                        % Write
                        ft_write_mri(p.sourcemodel.grid.seed.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.seed.names{ss}), mask_pow, 'dataformat', 'nifti_gz');
                        
                    else
                        
                        % Put in FieldTrip structure in order to write nifti
                        pow = [];
                        pow.R_avg = R_avg;
                        pow.pos = sourcemodel_grid.pos;
                        pow.unit = sourcemodel_grid.unit;
                        pow.dim = sourcemodel_grid.dim;
                        
                        % Write
                        fname = p.sourcemodel.grid.seed.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.seed.names{ss});
                        cfg_wn = [];
                        cfg_wn.filename =  fname(1:end-3);
                        cfg_wn.parameter = 'R_avg';
                        cfg_wn.precision = 'double';
                        cfg_wn.filetype = 'nifti';
                        ft_sourcewrite(cfg_wn, pow);
                        cmd = ['gzip ' fname];
                        system(cmd);
                        
                    end
                    
                end
            end
            
        end
        
        % Iterate over the seeds
        for ss = 1:size(coords, 1)
            
            % Register to MNI space if specified
            if p.sourcemodel.grid.seed.reg_mni && (~exist(p.sourcemodel.grid.seed.nii_mni(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.seed.names{ss}), 'file') || clobber == 1)
                
                % Register to MNI space
                cmd = ['antsApplyTransforms -d 3 -i ' p.sourcemodel.grid.seed.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.seed.names{ss}) ' -r ' p.headmodel.template_dir(sub_num) '/registrations/T1_template_registeredMNI152_Warped.nii.gz -t ' p.headmodel.template_dir(sub_num) '/registrations/T1_template_registeredMNI152_1Warp.nii.gz -t ' p.headmodel.template_dir(sub_num) '/registrations/T1_template_registeredMNI152_0GenericAffine.mat -o ' p.sourcemodel.grid.seed.nii_mni(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.seed.names{ss})];
                system(cmd);
                
            end
            
        end
        
    end
    
    %% SOURCEMODEL - GRID PSD
    
    % Only run if specified
    if contains(p.sourcemodel.method, 'grid') && p.sourcemodel.grid.psd.flag
        
        % Only proceed if last one doesn't exist
        if ~exist(p.sourcemodel.grid.psd.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.psd.foi(end, 1), p.sourcemodel.grid.psd.foi(end, 2)), 'file') || clobber == 1
            
            % Compute the pwelch one trial for sizing
            ts = bldata_store(cfg_trl.trl(1, 1):cfg_trl.trl(1, 2), 1);
            [~, psd_freq] = pwelch(zscore(ts), [], [], [], hdr.Fs, 'psd');
            
            % Compute the pwelch for each voxel, averaged over trials
            pxx = zeros(length(psd_freq), size(bldata_store, 2));
            for ii = 1:size(bldata_store, 2)
                
                % Compute the pwelch for each trial
                pxx_t = zeros(size(pxx, 1), length(data_filt.grid.trial));
                for tt = 1:length(data_filt.grid.trial)
                    ts = bldata_store(cfg_trl.trl(tt, 1):cfg_trl.trl(tt, 2), ii);
                    [pxx_t(:, tt), ~] = pwelch(zscore(ts), [], [], [], hdr.Fs, 'psd');
                end
                
                % Average over trials
                pxx(:, ii) = mean(pxx_t, 2);
                
            end
            
            % Iterate over the frequencies
            for ff_idx = 1:size(p.sourcemodel.grid.psd.foi, 1)
                
                if ~exist(p.sourcemodel.grid.psd.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.psd.foi(ff_idx, 1), p.sourcemodel.grid.psd.foi(ff_idx, 2)), 'file') || clobber == 1
                    
                    % Assign back the inside and average across the
                    % frequency bins
                    ind = psd_freq >= p.sourcemodel.grid.psd.foi(ff_idx, 1) & psd_freq <= p.sourcemodel.grid.psd.foi(ff_idx, 2);
                    pxx_avg = zeros(size(W_grid, 2), 1);
                    pxx_avg(sourcemodel_grid.inside) = mean(pxx(ind, :), 1);
                    
                    % Change back into mask if specified
                    if p.sourcemodel.grid.mask
                        
                        % Assign back to mask voxels
                        mask_pow = mask_mm;
                        t = pinv(mask_mm.transform);
                        for ii = 1:length(pxx_avg)
                            pos = sourcemodel_grid.pos(ii, :);
                            pos_vox = round(pos*t(1:3, 1:3) + t(1:3, 4)');
                            mask_pow.anatomy(pos_vox(1), pos_vox(2), pos_vox(3)) = pxx_avg(ii);
                        end
                        
                        % Write
                        ft_write_mri(p.sourcemodel.grid.psd.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.psd.foi(ff_idx, 1), p.sourcemodel.grid.psd.foi(ff_idx, 2)), mask_pow, 'dataformat', 'nifti_gz');
                        
                    else
                        
                        % Put in FieldTrip structure in order to write nifti
                        pow = [];
                        pow.pxx_avg = pxx_avg;
                        pow.pos = sourcemodel_grid.pos;
                        pow.unit = sourcemodel_grid.unit;
                        pow.dim = sourcemodel_grid.dim;
                        
                        % Write
                        fname = p.sourcemodel.grid.psd.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.psd.foi(ff_idx, 1), p.sourcemodel.grid.psd.foi(ff_idx, 2));
                        cfg_wn = [];
                        cfg_wn.filename =  fname(1:end-3);
                        cfg_wn.parameter = 'pxx_avg';
                        cfg_wn.precision = 'double';
                        cfg_wn.filetype = 'nifti';
                        ft_sourcewrite(cfg_wn, pow);
                        cmd = ['gzip ' fname];
                        system(cmd);
                        
                    end
                    
                end
            end
            
        end
        
        % Iterate over the frequency bands
        for ss = 1:size(p.sourcemodel.grid.psd.foi, 1)
            
            % Register to MNI space if specified
            if p.sourcemodel.grid.psd.reg_mni && (~exist(p.sourcemodel.grid.psd.nii_mni(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.psd.foi(ff_idx, 1), p.sourcemodel.grid.psd.foi(ff_idx, 2)), 'file') || clobber == 1)
                
                % Register to MNI space
                cmd = ['antsApplyTransforms -d 3 -i ' p.sourcemodel.grid.psd.nii(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.psd.foi(ff_idx, 1), p.sourcemodel.grid.psd.foi(ff_idx, 2)) ' -r ' p.headmodel.template_dir(sub_num) '/registrations/T1_template_registeredMNI152_Warped.nii.gz -t ' p.headmodel.template_dir(sub_num) '/registrations/T1_template_registeredMNI152_1Warp.nii.gz -t ' p.headmodel.template_dir(sub_num) '/registrations/T1_template_registeredMNI152_0GenericAffine.mat -o ' p.sourcemodel.grid.psd.nii_mni(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num), p.sourcemodel.grid.psd.foi(ff_idx, 1), p.sourcemodel.grid.psd.foi(ff_idx, 2))];
                system(cmd);
                
            end
            
        end
        
    end
    
    %% SOURCEMODEL - COORDINATES
    
    if contains(p.sourcemodel.method, 'coordinates')
        
        if ~exist(p.sourcemodel.coordinates.ve(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), 'file') || clobber == 1
            
            % Load the names
            temp = readtable(p.sourcemodel.coordinates.ve_names(sub_num), 'ReadVariableNames', false);
            
            % Now, compute the timeseries - sources x time x trials
            ts = zeros(size(W_coord, 2), length(data_filt.coordinates.time{1}), length(data_filt.coordinates.trial));
            for tt = 1:length(data_filt.coordinates.trial)
                ts_temp = zeros(size(temp, 1), length(data_filt.coordinates.time{1}), size(W_coord, 3));
                for jj = 1:size(W_coord, 3)
                    ts_temp(:, :, jj) = W_coord(:, :, jj)'*data_filt.coordinates.trial{tt};
                end
                if size(W_coord, 3) > 1
                    ts(:, :, tt) = sqrt(sum(ts_temp.^2, 3));
                else
                    ts(:, :, tt) = ts_temp;
                end
            end
            
            % Create save structure
            VE = [];
            VE.ts = ts;
            VE.pos_coord = sourcemodel_coord.pos;
            VE.pos_name = temp.Var2;
            VE.srate = data_filt.coordinates.fsample;
            VE.time = data_filt.coordinates.time{1};
            
            % Save
            parsave_struct(p.sourcemodel.coordinates.ve(p.subject(sub_num), p.session(sub_num), p.run(sub_num), p.task(sub_num)), VE, 'VE');
            
        end
        
    end
    
end
