function p = configure_study(p, modality)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p = opm_configure(p)
%
% Script to configure output directories and files for an OPM study
%
% Inputs:
%   p	(structure) Study configuration structure
%   
% Outputs:
%   p  (structure)  Study configuration structure with additional directory
%                   and file names
%
% Created by: Marlee M Vandewouw
% Created on: December 16, 2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ASSIGN DEFAULT

if nargin < 2
    modality = 'opm';
end

%% PIPELINE VERSION

p.pipeline.version = '8.0.0';

%% ADD DIRECTORIES TO PATH

% Add fieldtrip defaults
if isempty(which('ft_defaults'))
    addpath(p.fieldtrip_path);
    ft_defaults
end

% Add to the current path
addpath(genpath(p.pipeline_path));

% Remove ROINets+
rmpath([p.pipeline_path '/functions/etc/external/MEG-ROI-nets/+ROInets']);

%% EASY VARIABLES

% Get subject
p.subject = @(ind) ['sub-' p.subject_data.subject{ind}];

% Get session
p.session = @(ind) ['ses-' num2str(p.subject_data.session(ind), '%03d')];

% Get run
p.run = @(ind) ['run-' num2str(p.subject_data.run(ind), '%03d')];

% Get task
p.task = @(ind) p.subject_data.task{ind};

% Get the noise run
p.noise_run = @(ind) ['run-' num2str(p.subject_data.noise_run(ind), '%03d')];

% Get pipeline code
p.pipeline_code = [p.pipeline_name '_v' p.pipeline.version];

%% DIRECTORIES

% Raw data directory (BIDS structure)
p.directories.raw_dir = @(sub, ses) [p.directories.data_dir '/' sub '/' ses '/' modality];

% Anatomical directory
p.directories.anat_dir = @(sub) [p.directories.data_dir '/' p.subject(sub) '/' p.session(sub) '/anat'];

% Subject derivatives directory
p.directories.sub_dir = @(sub, ses) [p.directories.data_dir '/derivatives/' p.pipeline_code '/' sub '/' ses '/' modality];

%% RAW DATA

% Prefix
p.directories.raw_pre = @(sub, ses, run, task) [p.directories.raw_dir(sub, ses) '/' sub '_' ses '_task-' task '_' run '_'];

% Binary file
if strcmp(modality, 'opm')
    p.raw_data.bin = @(sub, ses, run, task) [p.directories.raw_pre(sub, ses, run, task) 'meg.bin'];
else
    p.raw_data.bin = @(sub, ses, run, task) [p.directories.raw_pre(sub, ses, run, task) 'meg.ds'];
end

% .json
p.raw_data.json = @(sub, ses, run, task) [p.directories.raw_pre(sub, ses, run, task) 'meg.json'];

% Positions
p.raw_data.positions = @(sub, ses, run, task) [p.directories.raw_pre(sub, ses, run, task) 'positions.tsv'];

% Channels
p.raw_data.positions = @(sub, ses, run, task) [p.directories.raw_pre(sub, ses, run, task) 'channels.tsv'];

% Sensor-to-mri
p.raw_data.sens2mri = @(sub, ses, run, task) [p.directories.raw_pre(sub, ses, run, task) 'sens2mri_transform.txt'];

% Head motion file
p.raw_data.headmotion = @(sub, ses, run, task) [p.directories.raw_pre(sub, ses, run, task) 'headmotion.mat'];

%% EPOCHING

% Prefix
p.directories.epoch_pre = @(sub, ses, run, task) [p.directories.sub_dir(sub, ses) '/' sub '_' ses '_task-' task '_' run '_1-epoching_desc-'];

% Trial configuration file
p.epoching.trl_cfg = @(sub, ses, run, task) [p.directories.epoch_pre(sub, ses, run, task) 'trials_cfg.mat'];

%% BAD CHANNELS

% Prefix
p.directories.channels_pre = @(sub, ses, run, task) [p.directories.sub_dir(sub, ses) '/' sub '_' ses '_task-' task '_' run '_2-channels_desc-'];

% Bad channels file
p.badchannels.bad_channels = @(sub, ses, run, task) [p.directories.channels_pre(sub, ses, run, task) 'bad_channels.txt']; 
p.badchannels.bad_channels_last = @(sub, ses, run, task) [p.directories.channels_pre(sub, ses, run, task) 'bad_channels_last.txt']; 

% PSD after removing bad channels
p.badchannels.bad_channels_psd = @(sub, ses, run, task) [p.directories.channels_pre(sub, ses, run, task) 'psd.png'];

%% FILTERING

% Prefix
p.directories.filter_pre = @(sub, ses, run, task) [p.directories.sub_dir(sub, ses) '/' sub '_' ses '_task-' task '_' run '_3-filtering_desc-'];

% PSD pre- and post- filtering
p.filtering.fig =  @(sub, ses, run, task) [p.directories.filter_pre(sub, ses, run, task) 'psd.png'];

% Bandstop widths
p.filter.bandstop_widths = @(sub, ses, run, task) [p.directories.filter_pre(sub, ses, run, task) 'widths.txt'];

%% ICA

% Prefix
p.directories.ica_pre = @(sub, ses, run, task) [p.directories.sub_dir(sub, ses) '/' sub '_' ses '_task-' task '_' run '_4-ica_desc-'];

% ICA file
p.ica.ica_components = @(sub, ses, run, task) [p.directories.ica_pre(sub, ses, run, task) 'components.mat']; 

% ICA figures
p.ica.ica_figure = @(sub, ses, run, task, num) [p.directories.ica_pre(sub, ses, run, task) 'ica_' num2str(num) '.png'];

% Bad components
p.ica.bad_components = @(sub, ses, run, task) [p.directories.ica_pre(sub, ses, run, task) 'bad_components.txt'];

% Make ICA files if column exists in the masterlist
if length(intersect('ica', p.subject_data.Properties.VariableNames)) == 1
    try
        if ~all(cellfun(@isempty, p.subject_data.ica))
            for ss = 1:size(p.subject_data, 1)
                if ~isempty(p.subject_data.ica{ss})
                    if ~exist(p.ica.ica_figure(p.subject(ss), p.session(ss), p.run(ss), p.task(ss), 1), 'file')
                        warning('Warning: ICAs are specified but files do not exist, skipping...');
                        continue;
                    end
                    comp = p.subject_data.ica{ss};
                    if isequal(comp, 'NaN')
                        comp = '';
                    end
                    fid = fopen(p.ica.bad_components(p.subject(ss), p.session(ss), p.run(ss), p.task(ss)), 'w');
                    fprintf(fid, "%s\n", comp);
                    fclose(fid);
                end
            end
        end
    end
end    

%% ARTIFACT REJECT

% Prefix
p.directories.artrej_pre = @(sub, ses, run, task) [p.directories.sub_dir(sub, ses) '/' sub '_' ses '_task-' task '_' run '_5-artrej_desc-'];

% Name of trial configuration prior to artifact reject
p.artifactreject.prev_trl = @(sub, ses, run, task) [p.directories.artrej_pre(sub, ses, run, task) 'trials_cfg_noartrej.mat'];

% Trial figures
p.artifactreject.trl_fig = @(sub, ses, run, task) [p.directories.artrej_pre(sub, ses, run, task) 'trials']; 

%% HEADMODEL

% Prefix
p.directories.headmod_pre = @(sub, ses, run, task) [p.directories.sub_dir(sub, ses) '/' sub '_' ses '_task-' task '_' run '_6-headmodel_desc-'];

% Template
if p.headmodel.custom
    p.headmodel.template_dir = @(sub) p.directories.anat_dir(sub);
    p.headmodel.template_nii = @(sub) [p.subject(sub) '_' p.session(sub) '_T1w.nii'];
else
    p.headmodel.template_dir = @(sub) [p.template_path '/' p.subject_data.template{sub}];
    p.headmodel.template_nii = @(sub) 'T1_template.nii';
end

% Tform (CTF -> MRI)
p.headmodel.tform = @(sub, ses, run, task) [p.directories.headmod_pre(sub, ses, run, task) 'tform.txt'];

% Inverse
p.headmodel.tform_inv = @(sub, ses, run, task) [p.directories.headmod_pre(sub, ses, run, task) 'tform_inv.txt'];

%% SOURCEMODEL

% Prefix
p.directories.sourcemodel_pre = @(sub, ses, run, task) [p.directories.sub_dir(sub, ses) '/' sub '_' ses '_task-' task '_' run '_6-sourcemodel_desc-'];

% Source model figure
p.sourcemodel.sourcemodel_fig = @(sub, ses, run, task, type) [p.directories.sourcemodel_pre(sub, ses, run, task) 'sourcefigure_' type '.png'];

% Coordinates

% Virtual electrode coordinates 
p.sourcemodel.coordinates.ve_locs = @(sub) [p.headmodel.template_dir(sub) '/' p.sourcemodel.coordinates.ve_atlas '_COG.txt'];

% Name for each coordinate
p.sourcemodel.coordinates.ve_names = @(sub) [p.headmodel.template_dir(sub) '/' p.sourcemodel.coordinates.ve_atlas '.txt'];

% Coordinate beamformer weights
p.sourcemodel.coordinates.weights = @(sub, ses, run, task) [p.directories.sourcemodel_pre(sub, ses, run, task) 'weights_coordinates_' p.sourcemodel.coordinates.ve_atlas '_' num2str(p.sourcemodel.coordinates.foi(1)) '_' num2str(p.sourcemodel.coordinates.foi(2)) 'Hz.mat'];

% Source model coordinates output (timeseries)
p.sourcemodel.coordinates.ve = @(sub, ses, run, task) [p.directories.sourcemodel_pre(sub, ses, run, task) 'VE_coordinates_' p.sourcemodel.coordinates.ve_atlas '_' num2str(p.sourcemodel.coordinates.foi(1)) '_' num2str(p.sourcemodel.coordinates.foi(2)) 'Hz.mat'];

% Grid

% Grid mask
p.sourcemodel.grid.mask_nii = @(sub) [p.headmodel.template_dir(sub) '/'  p.sourcemodel.grid.mask_atlas '_' num2str(p.sourcemodel.grid.res, '%d') 'mm.nii.gz'];

% Grid beamformer weights
p.sourcemodel.grid.weights = @(sub, ses, run, task) [p.directories.sourcemodel_pre(sub, ses, run, task) 'weights_grid_' num2str(p.sourcemodel.grid.foi(1)) '_' num2str(p.sourcemodel.grid.foi(2)) 'Hz.mat'];

% Source model grid output nifti
p.sourcemodel.grid.avgpow.nii = @(sub, ses, run, task, tw1, tw2, cond) [p.directories.sourcemodel_pre(sub, ses, run, task) 'gridpower_' num2str(p.sourcemodel.grid.foi(1)) '_' num2str(p.sourcemodel.grid.foi(2)) 'Hz_' num2str(tw1) '_' num2str(tw2) '_' cond];

% Source model grid output nifti in MNI space
p.sourcemodel.grid.avgpow.nii_mni = @(sub, ses, run, task, tw1, tw2, cond) [p.directories.sourcemodel_pre(sub, ses, run, task) 'gridpower_' num2str(p.sourcemodel.grid.foi(1)) '_' num2str(p.sourcemodel.grid.foi(2)) 'Hz_' num2str(tw1) '_' num2str(tw2) '_' cond '_MNI'];

% Grid timeseries

% Parcel to use to extract the timeseries for the max power voxel
p.sourcemodel.grid.ts.atlas = @(sub) [p.headmodel.template_dir(sub) '/'  p.sourcemodel.grid.ts.atlas_name '_' num2str(p.sourcemodel.grid.res, '%d') 'mm.nii.gz'];

% Names for parcels
p.sourcemodel.grid.ts.atlas_regions = @(sub) [p.headmodel.template_dir(sub) '/' p.sourcemodel.grid.ts.atlas_name '.txt'];

% Grid timeseries beamformer weights
p.sourcemodel.grid.ts.weights = @(sub, ses, run, task) [p.directories.sourcemodel_pre(sub, ses, run, task) 'weights_grid_ts_' num2str(p.sourcemodel.grid.ts.foi(1)) '_' num2str(p.sourcemodel.grid.ts.foi(2)) 'Hz.mat'];

% Timeseries of voxels with maximum power change
p.sourcemodel.grid.ts.ve = @(sub, ses, run, task) [p.directories.sourcemodel_pre(sub, ses, run, task) 'VE_coordinates_' p.sourcemodel.grid.ts.atlas_name '_maxpower' num2str(p.sourcemodel.grid.foi(1)) '_' num2str(p.sourcemodel.grid.foi(2)) 'Hz_' num2str(p.sourcemodel.grid.ts.foi(1)) '_' num2str(p.sourcemodel.grid.ts.foi(2)) 'Hz.mat'];

% Grid seeds

% File containing the coordinates
p.sourcemodel.grid.seed.coords = @(sub) [p.headmodel.template_dir(sub) '/' p.sourcemodel.grid.seed.coord_atlas '_COG.txt'];

% Correlation map
p.sourcemodel.grid.seed.nii = @(sub, ses, run, task, seed) [p.directories.sourcemodel_pre(sub, ses, run, task) 'seedcorr_' seed '_grid' num2str(p.sourcemodel.grid.res) '_' num2str(p.sourcemodel.grid.foi(1)) '-' num2str(p.sourcemodel.grid.foi(2)) 'Hz.nii.gz'];

% Correlation map: MNI
p.sourcemodel.grid.seed.nii_mni = @(sub, ses, run, task, seed) [p.directories.sourcemodel_pre(sub, ses, run, task) 'seedcorr_' seed '_grid' num2str(p.sourcemodel.grid.res) '_' num2str(p.sourcemodel.grid.foi(1)) '-' num2str(p.sourcemodel.grid.foi(2)) 'Hz_MNI.nii.gz'];

% Grid PSD

% Correlation map
p.sourcemodel.grid.psd.nii = @(sub, ses, run, task, freq1, freq2) [p.directories.sourcemodel_pre(sub, ses, run, task) 'psd_' num2str(freq1) '_' num2str(freq2) 'Hz_grid' num2str(p.sourcemodel.grid.res) '_' num2str(p.sourcemodel.grid.foi(1)) '-' num2str(p.sourcemodel.grid.foi(2)) 'Hz.nii.gz'];

% Correlation map: MNI
p.sourcemodel.grid.psd.nii_mni = @(sub, ses, run, task, freq1, freq2) [p.directories.sourcemodel_pre(sub, ses, run, task) 'psd_' num2str(freq1) '_' num2str(freq2) 'Hz_grid' num2str(p.sourcemodel.grid.res) '_' num2str(p.sourcemodel.grid.foi(1)) '-' num2str(p.sourcemodel.grid.foi(2)) 'Hz_MNI.nii.gz'];

%% CONNECTIVITY

% Prefix
p.directories.connectivity_pre = @(sub, ses, run, task) [p.directories.sub_dir(sub, ses) '/' sub '_' ses '_task-' task '_' run '_7-connectivity_desc-'];

% wPLI for the VE atlas
p.conenctivity.wpli_mat_ve_coordinates = @(sub, ses, run, task, cond, freq) [p.directories.connectivity_pre(sub, ses, run, task) 'VE_coordinates_' p.sourcemodel.coordinates.ve_atlas '_' num2str(p.sourcemodel.coordinates.foi(1)) '_' num2str(p.sourcemodel.coordinates.foi(2)) 'Hz_wpli_connectivity_f' freq '_' cond '.mat'];

% wPLI for the VE grid
p.conenctivity.wpli_mat_ve_grid = @(sub, ses, run, task, cond, freq) [p.directories.connectivity_pre(sub, ses, run, task) 'VE_coordinates_' p.sourcemodel.grid.ts.atlas_name '_maxpower' num2str(p.sourcemodel.grid.foi(1)) '_' num2str(p.sourcemodel.grid.foi(2)) 'Hz_' num2str(p.sourcemodel.grid.ts.foi(1)) '_' num2str(p.sourcemodel.grid.ts.foi(2)) 'Hz_maxpower_wpli_connectivity_f' freq '_' cond '.mat'];

% AEC for the VE atlas
p.conenctivity.aec_mat_ve_coordinates = @(sub, ses, run, task, cond, freq) [p.directories.connectivity_pre(sub, ses, run, task) 'VE_coordinates_' p.sourcemodel.coordinates.ve_atlas '_' num2str(p.sourcemodel.coordinates.foi(1)) '_' num2str(p.sourcemodel.coordinates.foi(2)) 'Hz_aec_connectivity_f' freq '_' cond '.mat'];

% AEC for the VE grid
p.conenctivity.aec_mat_ve_grid = @(sub, ses, run, task, cond, freq) [p.directories.connectivity_pre(sub, ses, run, task) 'VE_coordinates_' p.sourcemodel.grid.ts.atlas_name '_maxpower' num2str(p.sourcemodel.grid.foi(1)) '_' num2str(p.sourcemodel.grid.foi(2)) 'Hz_' num2str(p.sourcemodel.grid.ts.foi(1)) '_' num2str(p.sourcemodel.grid.ts.foi(2)) 'Hz_maxpower_aec_connectivity_f' freq '_' cond '.mat'];

%% TFR

% Prefix
p.directories.tfr_pre = @(sub, ses, run, task) [p.directories.sub_dir(sub, ses) '/' sub '_' ses '_task-' task '_' run '_8-tfr_desc-'];

% TFR file for region showing max power
p.tfr.tfr_grid_ts_file = @(sub, ses, run, task, region, cond) [p.directories.tfr_pre(sub, ses, run, task) 'TFR_maxpower' num2str(p.sourcemodel.grid.foi(1)) '_' num2str(p.sourcemodel.grid.foi(2)) 'Hz_' num2str(p.sourcemodel.grid.ts.foi(1)) '_' num2str(p.sourcemodel.grid.ts.foi(2)) 'Hz_' region '_' cond '.mat'];
p.tfr.tfr_grid_ts_stats_file = @(sub, ses, run, task, region, cond) [p.directories.tfr_pre(sub, ses, run, task) 'TFR_maxpower' num2str(p.sourcemodel.grid.foi(1)) '_' num2str(p.sourcemodel.grid.foi(2)) 'Hz_' num2str(p.sourcemodel.grid.ts.foi(1)) '_' num2str(p.sourcemodel.grid.ts.foi(2)) 'Hz_' region '_' cond '_stats.txt'];
p.tfr.tfr_grid_ts_fig = @(sub, ses, run, task, region, cond) [p.directories.tfr_pre(sub, ses, run, task) 'TFR_maxpower' num2str(p.sourcemodel.grid.foi(1)) '_' num2str(p.sourcemodel.grid.foi(2)) 'Hz_' num2str(p.sourcemodel.grid.ts.foi(1)) '_' num2str(p.sourcemodel.grid.ts.foi(2)) 'Hz_' region '_' cond '.png'];

% TFR file for region from a coordinate
p.tfr.tfr_coordinates_file = @(sub, ses, run, task, region, cond) [p.directories.tfr_pre(sub, ses, run, task) 'TFR_coordinates_' p.sourcemodel.coordinates.ve_atlas '_' num2str(p.sourcemodel.coordinates.foi(1)) '_' num2str(p.sourcemodel.coordinates.foi(2)) 'Hz_' region '_' cond '.mat'];
p.tfr.tfr_coordinates_stats_file = @(sub, ses, run, task, region, cond) [p.directories.tfr_pre(sub, ses, run, task) 'TFR_coordinates_' p.sourcemodel.coordinates.ve_atlas '_' num2str(p.sourcemodel.coordinates.foi(1)) '_' num2str(p.sourcemodel.coordinates.foi(2)) 'Hz_' region '_' cond '_stats.txt'];
p.tfr.tfr_coordinates_fig = @(sub, ses, run, task, region, cond) [p.directories.tfr_pre(sub, ses, run, task) 'TFR_coordinates_' p.sourcemodel.coordinates.ve_atlas '_' num2str(p.sourcemodel.coordinates.foi(1)) '_' num2str(p.sourcemodel.coordinates.foi(2)) 'Hz_' region '_' cond '.png'];

% Indices of top trials
p.tfr.tfr_grid_ts_ind_file = @(sub, ses, run, task, region, cond) [p.directories.tfr_pre(sub, ses, run, task) 'trialinds_maxpower' num2str(p.sourcemodel.grid.foi(1)) '_' num2str(p.sourcemodel.grid.foi(2)) 'Hz_' num2str(p.sourcemodel.grid.ts.foi(1)) '_' num2str(p.sourcemodel.grid.ts.foi(2)) 'Hz_' region '_' cond '.mat'];
p.tfr.tfr_coordinates_ind_file = @(sub, ses, run, task, region, cond) [p.directories.tfr_pre(sub, ses, run, task) 'trialinds_coordinates_' p.sourcemodel.coordinates.ve_atlas '_' num2str(p.sourcemodel.coordinates.foi(1)) '_' num2str(p.sourcemodel.coordinates.foi(2)) 'Hz_' region '_' cond '.mat'];

%% ENVELOPE 

p.envelope.ts_output_file = @(sub, ses, run, task, region, cond) [p.directories.tfr_pre(sub, ses, run, task) 'envelope_' region '_' cond '.mat'];

%% FOOOF

% FOOOF function
p.fooof.function = [p.pipeline_path '/functions/etc/custom/run_fooof.py'];

% Prefix
p.directories.fooof_pre = @(sub, ses, run, task) [p.directories.sub_dir(sub, ses) '/' sub '_' ses '_task-' task '_' run '_9-fooof_desc-'];

% Prefix, whole brain
p.directories.fooof_pre_whole_brain = @(sub, ses, run, task) [p.directories.sub_dir(sub, ses) '/' sub '_' ses '_task-' task '_' run '_9-fooof_desc-whole_brain'];

% FOOOF for the VE atlas
p.fooof.psd_mat_ve_coordinates = @(sub, ses, run, task) [p.directories.fooof_pre(sub, ses, run, task) 'VE_coordinates_' p.sourcemodel.coordinates.ve_atlas '_' num2str(p.sourcemodel.coordinates.foi(1)) '_' num2str(p.sourcemodel.coordinates.foi(2)) 'Hz_psd.mat'];

% FOOOF for the VE atlas: whole brain
p.fooof.psd_mat_ve_coordinates_whole_brain = @(sub, ses, run, task) [p.directories.fooof_pre_whole_brain(sub, ses, run, task) 'VE_coordinates_' p.sourcemodel.coordinates.ve_atlas '_' num2str(p.sourcemodel.coordinates.foi(1)) '_' num2str(p.sourcemodel.coordinates.foi(2)) 'Hz_psd.mat'];

