function opm_mri_mni_registrations(input_nii, output_dir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opm_registrations(p)
%
% Script to perform the co-registrations between MNI and subject MRI space
% required for the OPM pipeline. 
%
% Inputs:
%   input_nii       (string)    Input nifti (subject MRI)
%   template_nii    (string)    Template nifti (MNI)
%   output_dir      (string)    Output directory
%
% Created by: Marlee M Vandewouw
% Created on: February 7, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make registrations directory
output_dir_reg = [output_dir '/registrations'];
cmd = ['mkdir ' output_dir_reg];
system(cmd);

% First, perform brain extraction
cmd = ['antsBrainExtraction.sh -d 3 -a ' input_nii ' -e /usr/local/ANTs/install/data/MICCAI2012-Multi-Atlas-Challenge-Data/T_template0.nii.gz -m /usr/local/ANTs/install/data/MICCAI2012-Multi-Atlas-Challenge-Data/T_template0_BrainCerebellumProbabilityMask.nii.gz -f /usr/local/ANTs/install/data/MICCAI2012-Multi-Atlas-Challenge-Data/T_template0_BrainCerebellumRegistrationMask.nii.gz -o ' output_dir_reg '/T1_template_brain_'];
system(cmd);

% Now, perform the registration
cmd = ['antsRegistration --dimensionality 3 --float 0 --output [' output_dir_reg '/T1_template_registeredMNI152_,' output_dir_reg '/T1_template_registeredMNI152_Warped.nii.gz] --interpolation Linear --winsorize-image-intensities [0.005,0.995] --use-histogram-matching 0 --initial-moving-transform [/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz,' output_dir_reg '/T1_template_brain_BrainExtractionBrain.nii.gz,1] --transform Rigid[0.1] --metric MI[/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz,' output_dir_reg '/T1_template_brain_BrainExtractionBrain.nii.gz,1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox --transform Affine[0.1] --metric MI[/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz,' output_dir_reg '/T1_template_brain_BrainExtractionBrain.nii.gz,1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1  --smoothing-sigmas 3x2x1x0vox --transform SyN[0.1,3,0] --metric CC[/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz,' output_dir_reg '/T1_template_brain_BrainExtractionBrain.nii.gz,1,4] --convergence [100x70x50x20,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox'];
system(cmd);

% Convert to floats
cmd = ['mri_convert -odt float ' output_dir_reg '/T1_template_registeredMNI152_Warped.nii.gz ' output_dir_reg '/T1_template_registeredMNI152_Warped.nii.gz'];
system(cmd);
cmd = ['mri_convert -odt float ' output_dir_reg '/T1_template_registeredMNI152_1Warp.nii.gz ' output_dir_reg '/T1_template_registeredMNI152_1Warp.nii.gz'];
system(cmd);
cmd = ['mri_convert -odt float ' output_dir_reg '/T1_template_registeredMNI152_1InverseWarp.nii.gz ' output_dir_reg '/T1_template_registeredMNI152_1InverseWarp.nii.gz'];
system(cmd);

% Copy over the AAL atlas (in MNI space)
cmd = ['cp /d/mjt/9/users/marlee/atlases/aal/atlas/AAL_fix.nii ' output_dir '/'];
system(cmd);
cmd = ['gzip ' output_dir '/AAL_fix.nii'];
system(cmd);

% Convert to 1mm
cmd = ['flirt -in ' output_dir '/AAL_fix.nii.gz -ref /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz -applyisoxfm 1 -out ' output_dir '/AAL_fix.nii.gz -interp nearestneighbour'];
system(cmd);

% Apply the registrations
cmd = ['antsApplyTransforms -d 3 -i ' output_dir '/AAL_fix.nii.gz -r ' output_dir_reg '/T1_template_brain_BrainExtractionBrain.nii.gz -t [ ' output_dir_reg '/T1_template_registeredMNI152_0GenericAffine.mat, 1] -t ' output_dir_reg '/T1_template_registeredMNI152_1InverseWarp.nii.gz -n GenericLabel -o ' output_dir '/AAL90.nii.gz'];
system(cmd);

% Do the different resolutions
aal = ft_read_mri([output_dir '/AAL90.nii.gz']);
cfg_resample = [];
cfg_resample.method = 'nearest';
cfg_save = [];
cfg_save.parameter = 'anatomy';
cfg_save.filetype = 'nifti_gz';

% 1mm
cfg_resample.resolution = 1;
cfg_save.filename = [output_dir '/AAL90_1mm.nii.gz'];
aal_resample = ft_volumereslice(cfg_resample, aal);
ft_volumewrite(cfg_save, aal_resample);


% 2mm
cfg_resample.resolution = 2;
cfg_save.filename = [output_dir '/AAL90_2mm.nii.gz'];
aal_resample = ft_volumereslice(cfg_resample, aal);
ft_volumewrite(cfg_save, aal_resample);

% 5mm
cfg_resample.resolution = 5;
cfg_save.filename = [output_dir '/AAL90_5mm.nii.gz'];
aal_resample = ft_volumereslice(cfg_resample, aal);
ft_volumewrite(cfg_save, aal_resample);

% Remove the initial one
cmd = ['rm ' output_dir '/AAL_fix.nii.gz'];
system(cmd);

% Copy over the text file
cmd = ['cp /d/mjt/9/projects/OPM/opm_pipeline_templates/Adult/AAL90.txt ' output_dir '/'];
system(cmd);

% Split the parcels
cmd = ['mkdir ' output_dir '/split'];
system(cmd);
cmd = ['for i in {1..90}; do fslmaths ' output_dir '/AAL90.nii.gz -uthr ${i} -thr ${i} ' output_dir '/split/AAL90_${i}.nii.gz; done'];
system(cmd);
cmd = ['for i in {1..90}; do fslstats ' output_dir '/split/AAL90_${i}.nii.gz -c; done > ' output_dir '/AAL90_COG.txt'];
system(cmd);

% Remove the split directory
cmd = ['rm -r ' output_dir '/split'];
system(cmd);

% Copy over the FFG atlas (in MNI space)
cmd = ['cp /d/mjt/9/users/marlee/atlases/FFG/FFG_1mm.nii.gz ' output_dir '/FFG.nii.gz'];
system(cmd);

% Apply the registrations
cmd = ['antsApplyTransforms -d 3 -i ' output_dir '/FFG.nii.gz -r ' output_dir_reg '/T1_template_brain_BrainExtractionBrain.nii.gz -t [ ' output_dir_reg '/T1_template_registeredMNI152_0GenericAffine.mat, 1] -t ' output_dir_reg '/T1_template_registeredMNI152_1InverseWarp.nii.gz -n GenericLabel -o ' output_dir '/FFG.nii.gz'];
system(cmd);

% Do the different resolutions
aal = ft_read_mri([output_dir '/FFG.nii.gz']);
cfg_resample = [];
cfg_resample.method = 'nearest';
cfg_save = [];
cfg_save.parameter = 'anatomy';
cfg_save.filetype = 'nifti_gz';

% 1mm
cfg_resample.resolution = 1;
cfg_save.filename = [output_dir '/FFG_1mm.nii.gz'];
aal_resample = ft_volumereslice(cfg_resample, aal);
ft_volumewrite(cfg_save, aal_resample);

% 2mm
cfg_resample.resolution = 2;
cfg_save.filename = [output_dir '/FFG_2mm.nii.gz'];
aal_resample = ft_volumereslice(cfg_resample, aal);
ft_volumewrite(cfg_save, aal_resample);

% 5mm
cfg_resample.resolution = 5;
cfg_save.filename = [output_dir '/FFG_5mm.nii.gz'];
aal_resample = ft_volumereslice(cfg_resample, aal);
ft_volumewrite(cfg_save, aal_resample);

% Copy over the text file
cmd = ['cp /d/mjt/9/users/marlee/atlases/FFG/FFG.txt ' output_dir '/'];
system(cmd);

% Split the parcels
cmd = ['mkdir ' output_dir '/split'];
system(cmd);
cmd = ['for i in {1..2}; do fslmaths ' output_dir '/FFG.nii.gz -uthr ${i} -thr ${i} ' output_dir '/split/FFG_${i}.nii.gz; done'];
system(cmd);
cmd = ['for i in {1..2}; do fslstats ' output_dir '/split/FFG_${i}.nii.gz -c; done > ' output_dir '/FFG_COG.txt'];
system(cmd);

% Remove the split directory
cmd = ['rm -r ' output_dir '/split'];
system(cmd);



