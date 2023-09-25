function [MAGN,PHASE,GFACTOR,ARG]=NIFTI_NORDIC_VM(fn_magn_in, fn_phase_in, fn_out, ARG, fn_brainmask_in)
% fMRI
% fn_magn_in = 'name.nii.gz';
% fn_phase_in = 'name2.nii.gz';
% fn_brainmask_in = 'name3.nii.gz';
% fn_out = ['NORDIC_' fn_magn_in(1:end-7)];
% ARG.voxel_matching = 1;
% ARG.mask = 1;
% ARG.speed_factor = 5;
% NIFTI_NORDIC_VM(fn_magn_in,fn_phase_in,fn_out,ARG,fn_brainmask_in)
%
%
% file_input assumes 4D data
%
%   OPTIONS
%   ARG.DIROUT    VAL=      string        Default is empty
%   ARG.noise_volume_last   VAL  = num  specifiec volume from the end of the series
%                                          0 default
%
%   ARG.factor_error        val  = num    >1 use higher noisefloor <1 use lower noisefloor, 1 default
%   ARG.temporal_phase      val = [1 2 3]  1 was default, 3 now in dMRI due tophase errors in some data
%   ARG.NORDIC              val = [0 1]    1 Default
%   ARG.MP                  val = [0 1 2]  1 NORDIC gfactor with MP estimation.
%                                          2 MP without gfactor correction
%                                          0 default
%   ARG.kernel_size_gfactor val = [val1 val2 val], defautl is [14 14 1]
%   ARG.kernel_size_PCA     val = [val1 val2 val], default is val1=val2=val3;
%                                                  ratio of 11:1 between spatial and temproal voxels
%   ARG.magnitude_only      val =[] or 1.  Using complex or magntiude only. Default is []
%                                          Function still needs two inputs but will ignore the second%
%   ARG.save_add_info       val =[0 1];  If it is 1, then an additonal matlab file is being saved with degress removed etc.
%                                         default is 0
%   ARG.make_complex_nii    if the field exist, then the phase is being saved in a similar format as the input phase
%
%   ARG.phase_filter_width  val = [1... 10]  Specifiec the width of the smoothing filter for the phase
%                                         default is now 3
%
%   ARG.save_gfactor_map   val = [1 2].  1, saves the RELATIVE gfactor, 2 saves the
%                                            gfactor and does not complete the NORDIC processing
%   ARG.voxel_matching     val = [0 1] 0, non-active, use default local patching 1, active
%
%   ARG.speed_factor       val = positive integer (suggested between 1 and 10). Skip val reference pixels to speed up denoising via pixel matching
%   Only works when ARG.voxel_matching = 1
%
%   ARG.mask            val = [0 1] 1, the user provides a 3D binary brain mask as last input  0, VMNORDIC uses an empty mask
%

% TODO
% Scaling relative to the width of the MP spectrum, if one wants to be conservative
%
% 4/15/21 swapped the uint16 and in16 for the phase
%
% VERSION 4/22/2021
% Copyright Board of Regents, University of Minnesota, 2022

if nargin == 4
    fn_brainmask_in = [];
end

if ~isfield(ARG,'noise_volume_last')
    ARG.noise_volume_last = 0; % 0 there is no noise volume {0 1 2 ...}, 1 if last volume is noise volume
end

if ~isfield(ARG,'factor_error')
    ARG.factor_error = 1.0; % error in gfactor estimatetion. >1 use higher noisefloor <1 use lower noisefloor
end

if ~isfield(ARG,'NORDIC') && ~isfield(ARG,'MP')
    ARG.NORDIC = 1; % threshold based on Noise
    ARG.MP = 0; % threshold based on Marchencko-Pastur
elseif ~isfield(ARG,'NORDIC') % MP selected
    if ARG.MP == 1
        ARG.NORDIC = 0;
    else
        ARG.NORDIC = 1;
    end
elseif ~isfield(ARG,'MP') % NORDIC selected
    if ARG.NORDIC == 1
        ARG.MP = 0;
    else
        ARG.MP = 1;
    end
end

if ~isfield(ARG,'NORDIC_patch_overlap')
    ARG.NORDIC_patch_overlap = 2; %
end

if ~isfield(ARG,'gfactor_patch_overlap')
    ARG.gfactor_patch_overlap = 2; %
end

if ~isfield(ARG,'kernel_size_gfactor')
    ARG.kernel_size_gfactor = []; % default is [14 14 1]
end

if ~isfield(ARG,'kernel_size_PCA')
    ARG.kernel_size_PCA = []; % default is 11:1 ratio
end

if ~isfield(ARG,'magnitude_only') % if legacy data
    ARG.magnitude_only = 0; %
end

if ~isfield(ARG,'save_gfactor_map') % save out a map of a relative gfactor
    ARG.save_gfactor_map = []; %
end

if ~isfield(ARG,'voxel_matching')
    ARG.voxel_matching = [];
end

if ~isfield(ARG,'mask')
    ARG.mask = 0; %
end

if ~isfield(ARG,'speed_factor')
    ARG.speed_factor = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ARG.magnitude_only == 0 % complex data input    
    try
        info_phase = niftiinfo(fn_phase_in);
        info_magn = niftiinfo(fn_magn_in);
    catch
        disp('The niftiinfo fails at reading the header');
    end
    
    I_M = abs(single(niftiread(fn_magn_in)));
    I_P = single(niftiread(fn_phase_in));
    info_phase.Datatype = class(I_P);
    info_magn.Datatype = class(I_M); 
    
    % Here, we combine magnitude and phase data into complex form
    fprintf('Phase should be -pi to pi...\n')
    
    % scale the phase
    phase_range = max(I_P(:));
    phase_range_min = min(I_P(:));       
    range_norm = phase_range - phase_range_min;
    range_center = (phase_range + phase_range_min)/range_norm * 1/2;
    I_P = (I_P./range_norm - range_center)*2*pi;
    II = I_M .* exp(1i*I_P); % construct compled data input II
    
    fprintf('Phase data range is %.2f to %.2f\n', min(I_P(:)), max(I_P(:)))
else % magnitude only input data
    try
        info_magn = niftiinfo(fn_magn_in);
    catch
        disp('The niftiinfo fails at reading the header');
    end
    
    I_M = abs(single(niftiread(fn_magn_in)));    
    info_magn.Datatype = class(I_M);

end

if ARG.magnitude_only == 1
    II = single(I_M);
end

TEMPVOL = abs(II(:,:,:,1));
ARG.ABSOLUTE_SCALE = min(TEMPVOL(TEMPVOL ~= 0));
II = II./ARG.ABSOLUTE_SCALE;

KSP2 = II;

if isempty(ARG.kernel_size_gfactor) || size(ARG.kernel_size_gfactor,2) < 4
    KSP2 = KSP2(:,:,:,1:min(90,end)); % takes first 90 volumes or less for gfactors estimation, but should be at least 30 volumes
else % fourth dimension in ARG.kernel_size_gfactor is the number of dynamics the gfactor estimation should consider
    KSP2 = KSP2(:,:,:,1:min(ARG.kernel_size_gfactor(4),end));
end

KSP2(isnan(KSP2)) = 0;
KSP2(isinf(KSP2)) = 0;
DIMS = size(KSP2);
% JCWS option to skip or perform gfactor estimation

if ARG.gfactorcorr == 1
    if isempty(ARG.kernel_size_gfactor)
        ARG.kernel_size_gfactor = [14 14 1];
    end
    QQ.KSP_processed = zeros(1,DIMS(1) - ARG.kernel_size_gfactor(1));
    ARG.patch_average = 0;
    ARG.patch_average_sub = ARG.gfactor_patch_overlap;
    
    ARG.LLR_scale = 0;
    ARG.NVR_threshold = 1;
    ARG.soft_thrs = 10; % MPPCa (When Noise varies), for gfactor estimation, %ARG.soft_thrs = []; % NORDIC (When noise is flat)
    
    KSP_recon = zeros(DIMS);
    KSP_weight =  zeros(DIMS(1:3));
    NOISE = KSP_weight; 
    Component_threshold = KSP_weight;
    energy_removed = KSP_weight;
    SNR_weight = KSP_weight;
    QQ.KSP_processed = zeros(1,DIMS(1) - ARG.kernel_size_gfactor(1));
    
    if ARG.patch_average == 0
        KSP_processed = zeros(size(QQ.KSP_processed));
        for nw1 = 2:max(1,floor(ARG.kernel_size_gfactor(1)/ARG.patch_average_sub))
            KSP_processed(1,nw1 : max(1,floor(ARG.kernel_size_gfactor(1)/ARG.patch_average_sub)):end) = 2;
        end
        KSP_processed(end) = 0; 
        QQ.KSP_processed = KSP_processed;
    end
  
    disp(QQ.KSP_processed)
    disp(num2str(size(QQ.KSP_processed,2)))    
    
    ARG.kernel_size = ARG.kernel_size_gfactor; % main function sub_LLR_Processing needs a ARG.kernel_size variable
    disp('estimating g-factor ...')    
    for n1 = 1:size(QQ.KSP_processed,2)% nl is the index that moves the patch around! for matrixsize in X = Nx  Nx - patchsize steps
        [KSP_recon,~,KSP_weight,NOISE,Component_threshold,energy_removed,SNR_weight] = sub_LLR_Processing(KSP_recon,KSP2,ARG,n1,QQ,KSP_weight,NOISE,Component_threshold,energy_removed,SNR_weight) ; % save all files
    end
    
    KSP_recon = KSP_recon./repmat((KSP_weight),[1 1 1 DIMS(4)]); % normalize by KSP_weight, which is duplicated N timepoints of reduced KSP input
    
    ARG.NOISE = sqrt(NOISE./KSP_weight);
    ARG.Component_threshold = Component_threshold./KSP_weight;
    ARG.energy_removed = energy_removed./KSP_weight;
    ARG.SNR_weight = SNR_weight./KSP_weight;
    
    disp('completed estimating g-factor')    
    ARG.gfactor = ARG.NOISE;     
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (ARG.save_gfactor_map == 1) || (ARG.save_gfactor_map == 2)
        
        gfactor_IMG = abs( ARG.gfactor); 
        gfactor_IMG(isnan(gfactor_IMG)) = 0;
        
        if strcmp(info_magn.Datatype,'uint16')
            gfactor_IMG = uint16(gfactor_IMG);
        elseif strcmp(info_magn.Datatype,'int16')
            gfactor_IMG = int16(gfactor_IMG);
        else
            gfactor_IMG = single(gfactor_IMG);
        end
        
        GFACTOR = gfactor_IMG; % for function output parse
        niftiwrite(gfactor_IMG,[fn_out 'gfactor.nii'],'Compressed',true)
        
        if ARG.save_gfactor_map == 2 % only estimates gfactor and return omitting NORDIC
            return
        end
    end
elseif ARG.gfactorcorr == 0
    disp('Not performing gfactor estimation and normalization, set it to 1s')
    ARG.gfactor = ones(DIMS(1:3));
    GFACTOR = ARG.gfactor;   
end

KSP2 = II;% after gfactor estimation replace KSP2 with original input again

% gfactor normalization, this also makes the noise have standard deviation of 1 ! therefor eigenvalue threshold estimation from random matrices with zero mean, and 1 std normal distributed noise canbe used!
KSP2 = bsxfun(@rdivide, KSP2, ARG.gfactor);

if ARG.noise_volume_last > 0
    KSP2_NOISE = KSP2(:,:,:,end + 1 - ARG.noise_volume_last); % extract last volume as NOISE volume (ARG.noise_volume_last = 1 then)
end

KSP2(isnan(KSP2)) = 0;
KSP2(isinf(KSP2)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ARG.noise_volume_last > 0 % extract noise value from noise scan
    tmp_noise = KSP2_NOISE;
    tmp_noise(isnan(tmp_noise)) = 0; % cleanup
    tmp_noise(isinf(tmp_noise)) = 0; % cleanup
    ARG.measured_noise = std(tmp_noise(tmp_noise ~= 0 )); % sqrt(2) for real and complex JWCS: WHAT TO DO WITH 'std' from complex noise or magnitude niuse????
else
    ARG.measured_noise = 1; % IF COMPLEX DATA   , CHECK!!!!
end

if ~isfield(ARG,'use_magn_for_gfactor') && (isempty(ARG.magnitude_only) || ARG.magnitude_only == 0) %% WOULD THIS BE THE ISSUE & replaced by |
    ARG.measured_noise = ARG.measured_noise/sqrt(2); % scale noise with 1/sqrt(2), only for complex data
    disp(num2str(ARG.measured_noise))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine Patch size
if isempty(ARG.kernel_size_PCA)
    ARG.kernel_size = repmat(round((DIMS(4)*11)^(1/3)),1,3);  % use spatial temporal dimension size ratio of 11:1: MxMxM=11*Ntimepoints = M = cuberoot (11xNtimepoints)
else
    ARG.kernel_size = ARG.kernel_size_PCA ;
end

if DIMS(3) <= ARG.kernel_size(3) % Number of slices is less than cubic kernel
    ARG.kernel_size = repmat(round((DIMS(4)*11/DIMS(3) )^(1/2)),1,2);
    ARG.kernel_size(3) = DIMS(3);
end
ARG.kernel_size_total = ARG.kernel_size(2)*ARG.kernel_size(1)*ARG.kernel_size(3);
disp(['Patch size for standard NORDIC = ' num2str(ARG.kernel_size) ' = ' num2str(ARG.kernel_size_total)]);

QQ.KSP_processed = zeros(1,DIMS(1)-ARG.kernel_size(1));
ARG.patch_average = 0;
ARG.patch_average_sub = ARG.NORDIC_patch_overlap ;

if ARG.MP > 0 % for MPPCA
    ARG.kernel_size = [7 7 7];
    ARG.patch_average_sub = 7;
    ARG.soft_thrs = 10; % MPPCa (When Noise varies)
end
ARG.soft_thrs = []; % NORDIC (When noise is flat)
ARG.LLR_scale = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eigenvalue estimation of 10 random matrices, zero mean normal distributes

ARG.NVR_threshold = 0;
for ntmp = 1:10
    [~,S,~] = svd(randn(prod(ARG.kernel_size),DIMS(4)));
    ARG.NVR_threshold = ARG.NVR_threshold+S(1,1); % store (sum) 10 highest eigenvalues, later divide by to to obtain mean
end

% divide by 10, and some correction factor due to bias in g-factor estimation
if ARG.magnitude_only ~= 1 % 4/29/2021
    ARG.NVR_threshold = ARG.NVR_threshold/10*sqrt(2)* ARG.measured_noise*ARG.factor_error; % sqrt(2) due to complex, ARG.factor_error due to understimate of g-factor
else
    ARG.NVR_threshold = ARG.NVR_threshold/10*ARG.measured_noise*ARG.factor_error; % no sqrt(2) due to magnitude, ARG.factor_error due to understimate of g-factor
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ARG.voxel_matching ~= 1 % standard NORDIC
    KSP_recon = zeros(DIMS);  
    KSP_weight =  zeros(DIMS(1:3));
    NOISE = KSP_weight;
    Component_threshold = KSP_weight;
    energy_removed = KSP_weight;
    SNR_weight = KSP_weight;
    
    QQ.KSP_processed = zeros(1,DIMS(1) - ARG.kernel_size(1));
    
    if ARG.patch_average == 0
        KSP_processed = zeros(size(QQ.KSP_processed));
        for nw1 = 2:max(1,floor(ARG.kernel_size(1)/ARG.patch_average_sub))
            KSP_processed(1,nw1 : max(1,floor(ARG.kernel_size(1)/ARG.patch_average_sub)):end) = 2;
        end
        KSP_processed(end) = 0; % disp
        QQ.KSP_processed = KSP_processed;
    end
    disp('starting NORDIC ...')
    
    size(QQ.KSP_processed,2)
    for n1 = 1:size(QQ.KSP_processed,2)
        [KSP_recon,~,KSP_weight,NOISE,Component_threshold,energy_removed,SNR_weight] = sub_LLR_Processing(KSP_recon,KSP2,ARG,n1,QQ,KSP_weight,NOISE,Component_threshold,energy_removed,SNR_weight) ; % save all files
    end
    KSP_recon = KSP_recon./repmat((KSP_weight),[1 1 1 DIMS(4)]); % Assumes that the combination is with N instead of sqrt(N). Works for NVR not MPPCA
    ARG.NOISE = sqrt(NOISE./KSP_weight);
    ARG.Component_threshold = Component_threshold./KSP_weight;
    ARG.energy_removed = energy_removed./KSP_weight;
    ARG.SNR_weight = SNR_weight./KSP_weight;
    IMG_denoised = KSP_recon;
    disp('completing NORDIC ...')
    
else % VM-NORDIC
    if ARG.mask == 1
        brain_mask = single(niftiread(fn_brainmask_in));
    else
        brain_mask = [];
    end
    
    disp('starting non-local voxel matching VM-NORDIC...')    
    [IMG_denoised,ARG] = voxel_matchNLLR(KSP2,ARG,brain_mask);    
    disp('completing NORDIC with non-local voxel matching (VM_NORDIC)...')
    
    % write visiting map
    niftiwrite(single(ARG.weight_image),[fn_out 'visitingmap.nii'],'Compressed',true)
end

if isfield(ARG,'save_residual_matlab')
    if ARG.save_residual_matlab == 1
        Residual = KSP2-KSP_recon;
        save([ fn_out 'RESIDUAL.mat' ],'Residual','-v7.3')
    end
end

% gfactor renormalization
IMG_denoised = bsxfun(@times, IMG_denoised, ARG.gfactor);

IMG_denoised = IMG_denoised.*ARG.ABSOLUTE_SCALE;
IMG_denoised(isnan(IMG_denoised)) = 0;

if ARG.make_complex_nii == 1
    IMG_denoised_M = abs(IMG_denoised);
    IMG_denoised_M(isnan(IMG_denoised_M)) = 0;
   
    if strcmp(info_magn.Datatype,'uint16')
        IMG_denoised_M = uint16(IMG_denoised_M);
    elseif strcmp(info_magn.Datatype,'int16')
        IMG_denoised_M = int16(IMG_denoised_M);
    else
        IMG_denoised_M = single(IMG_denoised_M);
    end
    
    MAGN = IMG_denoised_M;
    niftiwrite(MAGN,[fn_out 'magn.nii'],info_magn,'Compressed',true)
    
    IMG_denoised_P = angle(IMG_denoised);
    IMG_denoised_P = (IMG_denoised_P/(2*pi) + range_center)*range_norm;
    
    if strcmp(info_phase.Datatype,'uint16')
        IMG_denoised_P = uint16(IMG_denoised_P);
    elseif strcmp(info_phase.Datatype,'int16')
        IMG_denoised_P = int16(IMG_denoised_P);
    else
        IMG_denoised_P = single((IMG_denoised_P));
    end
    PHASE = IMG_denoised_P;
    niftiwrite(PHASE,[fn_out 'phase.nii'],info_phase,'Compressed',true)
    
elseif ARG.make_complex_nii == 0
    
    IMG_denoised_M = abs(IMG_denoised);
    IMG_denoised_M(isnan(IMG_denoised_M)) = 0;

    if strcmp(info_magn.Datatype,'uint16')
        IMG_denoised_M = uint16(IMG_denoised_M);
    elseif strcmp(info_magn.Datatype,'int16')
        IMG_denoised_M = int16(IMG_denoised_M);
    else
        IMG_denoised_M = single(IMG_denoised_M);
    end
    
    MAGN = IMG_denoised_M;
    PHASE = [];
    niftiwrite(MAGN,[fn_out 'magn.nii'],info_magn,'Compressed',true)
end

if isfield(ARG,'save_add_info')
    if ARG.save_add_info == 1
        disp('saving additional info')
        save([fn_out 'ARG.mat' ],'ARG','-v7.3')
    end
end
return
end

function [KSP_recon,KSP2,KSP2_weight,NOISE, Component_threshold,energy_removed,SNR_weight] = sub_LLR_Processing(KSP_recon,KSP2,ARG,n1,QQ,KSP2_weight,NOISE,Component_threshold,energy_removed,SNR_weight)

Iwindowx = 1:ARG.kernel_size(1)+(n1-1); %patch selection window indices range in X direction

if QQ.KSP_processed(1,n1) ~= 1 && QQ.KSP_processed(1,n1) ~= 3 % not being processed also not completed yet
    if QQ.KSP_processed(1,n1) == 2 % processed but not added.
        try  % try to load otherwise go to next slice
         load([ARG.filename  'slice' num2str(n1)  '.mat'],'DATA_full2')
        catch
          QQ.KSP_processed(1,n1)=0;  % identified as bad file and being identified for reprocessing
          return
        end
        QQ.KSP_processed(1,n1) = 0; %         
    end    
    if QQ.KSP_processed(1,n1) ~= 2        
        QQ.KSP_processed(1,n1) = 1 ; % STARTING
        KSP2a = KSP2(Iwindowx,:,:,:); % select kernelsize(1), KSP2a is a patchblock confined to kernelsize X and crossing all Y and Z, timepoints range
        lambda = ARG.LLR_scale*ARG.NVR_threshold;
        disp(['lambla threshold = ' num2str(lambda)])        
        if ARG.patch_average == 0
            KSP2_weight_patchtmp = KSP2_weight(Iwindowx,:,:,:);
            NOISE_patchtmp = NOISE(Iwindowx,:,:,:);
            Component_threshold_patchtmp = Component_threshold(Iwindowx,:,:,:);
            energy_removed_patchtmp = energy_removed(Iwindowx,:,:,:);
            SNR_weight_tmp = SNR_weight(Iwindowx,:,:,:);
            
            patch_avg = 1;
            [KSP2_patchtmp_update,KSP2_weight_patchtmp, NOISE_patchtmp, Component_threshold_patchtmp, energy_removed_tmp, SNR_weight_tmp] = ...
                subfunction_loop_for_NVR_avg_update(KSP2a,ARG.kernel_size(2),ARG.kernel_size(3),lambda,patch_avg,ARG.soft_thrs,KSP2_weight_patchtmp,ARG,NOISE_patchtmp,Component_threshold_patchtmp,energy_removed_patchtmp,SNR_weight_tmp);
            
            KSP2_weight(Iwindowx,:,:,:) = KSP2_weight_patchtmp;
            NOISE(Iwindowx,:,:,:) = NOISE_patchtmp;
            Component_threshold(Iwindowx,:,:,:) = Component_threshold_patchtmp;
            energy_removed(Iwindowx,:,:,:) = energy_removed_tmp;
            SNR_weight(Iwindowx,:,:,:) = SNR_weight_tmp;
        end        
    end    
    KSP_recon(Iwindowx,1:size(KSP2_patchtmp_update,2),:,:) = KSP_recon(Iwindowx,1:size(KSP2_patchtmp_update,2),:,:) + KSP2_patchtmp_update;
    QQ.KSP_processed(1,n1) = 3 ;
end
return
end

function [KSP2_patchtmp_update, KSP2_weight,NOISE,KSP2_tmp_update_threshold,energy_removed,SNR_weight] = subfunction_loop_for_NVR_avg_update(KSP2a,w2,w3,lambda,patch_avg, soft_thrs,KSP2_weight,ARG,NOISE,KSP2_tmp_update_threshold,energy_removed,SNR_weight)
DIMS_block = size(KSP2a); % KSP2a is a patchblock confined to kernelsize X and crossing all Y and Z, timepoints range
KSP2_patchtmp_update = zeros(DIMS_block);
patch_scale = 1;

Ipatchstepy = max(1,floor(w2/ARG.patch_average_sub));
Ipatchmaxy = DIMS_block(2) - w2 + 1;

Ipatchstepz = max(1,floor(w3/ARG.patch_average_sub));
Ipatchmaxz = DIMS_block(3) - w3+1;

for n2 = [1 : Ipatchstepy : Ipatchmaxy, Ipatchmaxy] % move patch around in Y direction in steps of kernelsize y (w2)
    for n3 = [1 : Ipatchstepz : Ipatchmaxz, Ipatchmaxz] % move patch around in Z direction in steps of kernelsize z (w3)
        Iwindowy = (1:w2)+(n2-1); %  moving kernel range and index in y
        Iwindowz = (1:w3)+(n3-1); %  moving kernel range and index  in z
       
        % this is the patch: dimensions kernelsize(1) x kernelsize(2) x kernelsize(3) x timepoints
        KSP2_patchtmp = KSP2a(:, Iwindowy, Iwindowz,:);
        tmp1 = reshape(KSP2_patchtmp,[],DIMS_block(4)); % this is the Casorati matrix of the patch
        
        [U,S,V] = svd(tmp1,'econ');
        S = diag(S);
        
        [idx] = sum(S < lambda);
        if isempty(soft_thrs)% NORDIC denoising, thresholding eigenvalues            
            energy_scrub = sqrt(sum(S)).\sqrt(sum(S(S < lambda)));
            S(S < lambda) = 0;
            t = idx;
        elseif soft_thrs ~= 10
            S = S - lambda*soft_thrs;
            S(S < 0) = 0;
            energy_scrub = 0;
            t = 1;
        elseif soft_thrs == 10 % USING MPPCA
            centering = 0;
            MM = size(tmp1,1); % Correction for some zero entries
            if MM > 0
                NNN = size(tmp1,2);
                R = min(MM, NNN);
                scaling = (max(MM, NNN) - (0:R-centering-1)) / NNN;
                scaling = scaling(:);
                vals = S;
                vals = (vals).^2 / NNN;
                % First estimation of Sigma^2; Eq 1 from ISMRM presentation
                csum = cumsum(vals(R-centering:-1:1));
                cmean = csum(R-centering:-1:1)./(R-centering:-1:1)';
                sigmasq_1 = cmean./scaling;
                % Second estimation of Sigma^2; Eq 2 from ISMRM presentation
                gamma = (MM - (0:R-centering-1)) / NNN;
                rangeMP = 4*sqrt(gamma(:));
                rangeData = vals(1:R-centering) - vals(R-centering);
                sigmasq_2 = rangeData./rangeMP;
                t = find(sigmasq_2 < sigmasq_1, 1);
                idx = size(S(t:end),1) ;
                energy_scrub = sqrt(sum(S)).\sqrt(sum(S(t:end)));
                S(t:end) = 0;
            else % all zero entries
                t = 1;
                energy_scrub = 0;
                sigmasq_2 = 0;
            end
        else
            S(max(1,end - floor(idx*soft_thrs)):end) = 0;
        end
        
        tmp1 = U*diag(S)*V';
        tmp1 = reshape(tmp1,size(KSP2_patchtmp)); % reshape denoised by SVT tmp1 back to a 4D patch
        
        if patch_scale == 1
        else
            patch_scale = size(S,1)-idx;
        end
        
        if isempty(t)
            t = 1;
        end % threshold removed all.
        
        if patch_avg == 1
            KSP2_patchtmp_update(:,Iwindowy, Iwindowz,:) = KSP2_patchtmp_update(:,Iwindowy, Iwindowz,:) + patch_scale*tmp1;
            KSP2_weight(:,Iwindowy, Iwindowz,:) = KSP2_weight(:,Iwindowy, Iwindowz,:) + patch_scale;
            KSP2_tmp_update_threshold(:,Iwindowy, Iwindowz,:) = KSP2_tmp_update_threshold(:,Iwindowy, Iwindowz,:) + idx;
            energy_removed(:,Iwindowy, Iwindowz,:) = energy_removed(:,Iwindowy, Iwindowz,:) + energy_scrub;
            SNR_weight(:,Iwindowy, Iwindowz,1) = SNR_weight(:,Iwindowy, Iwindowz,1) + S(1)./S(max(1,t-1));
            if soft_thrs == 10
                NOISE(1:DIMS_block(1),Iwindowy, Iwindowz,1) = NOISE(1:DIMS_block(1),Iwindowy, Iwindowz,1) + sigmasq_2(t);
            end
        end
    end
end
return
end

