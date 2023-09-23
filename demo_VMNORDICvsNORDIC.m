%% DEMO VM-NORDIC vs NORDIC, SIERO 2023
%  In the lines below, specify the file names for the magnitude, phase and 
%  mask nifti files included in the folder. 



ARG.DIRIN              = ['DEMO/files/'];
ARG.DIROUT             = ['DEMO/files/'];
fn_magn_in             = [ARG.DIRIN 'magn.nii.gz'];
fn_phase_in            = [ARG.DIRIN 'phase.nii.gz'];
fn_brainmask_in        = [ARG.DIRIN 'mask.nii.gz'];

ARG.temporal_phase     = 1;
ARG.phase_filter_width = 10;
ARG.mask               = 1; % = 1 the user provides a 3D brain mask. = 0 no mask is provided and VM_NORDIC uses a empty  mask
ARG.speed_factor       = 5; % integer between 1 and 10. Speeds up the VM NORDIC denoising
ARG.save_gfactor_map   = 1;
ARG.save_add_info      = 1;
ARG.gfactorcorr        = 1;
ARG.timepoints         = []; % timepoint selection to base similarity (voxel matching) on, if empty all timepoints are used
ARG.magnitude_only     = 0;
ARG.make_complex_nii   = 1 - ARG.magnitude_only ; % change when ARG.magnitude_only is ALTERED

% VM-NORDIC
fn_out                 = [ARG.DIROUT 'NORDIC_VM_'];
ARG.voxel_matching     = 1; % = 1 activate VM_NORDIC
[MAGN_VMND,PHASE_VMND,GFACTOR_VMND,ARG_VMND] = NIFTI_NORDIC_VM(fn_magn_in, fn_phase_in, fn_out, ARG, fn_brainmask_in);

% NORDIC
fn_out                 = [ARG.DIROUT 'NORDIC_'];
ARG.voxel_matching     = 0; % = 1 activate VM_NORDIC
[MAGN_ND,PHASE_ND,GFACTOR_ND,ARG_ND] = NIFTI_NORDIC_VM(fn_magn_in, fn_phase_in, fn_out, ARG, fn_brainmask_in);



