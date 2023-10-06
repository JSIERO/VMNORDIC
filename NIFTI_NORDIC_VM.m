function [MAGN,PHASE,GFACTOR,ARG]=NIFTI_NORDIC_VM(ARG, fn_magn_in, fn_phase_in, fn_out, fn_brainmask_in)
%% Voxel Matching NORDIC: VM-NORDIC denoising for fMRI - Siero, Nigi 2023 
%% based on Nigi, Siero 'Improved NORDIC denoising for submillimetre BOLD fMRI using patch formation via non-local pixels similarity - pixel-matching (PM) NORDIC' ISMRM Toronto programnr 480
% also includes a lean version of standard NORDIC and MPPCALL: (based on Moeller NI2021, Vizioli NatComm2021,https://github.com/SteenMoeller/NORDIC_Raw 
% 
% example settings commonly used for VMNORDIOC
% ARG.fn_magn_in            = 'magn.nii.gz';
% ARG.fn_phase_in           = 'phase.nii.gz';
% ARG.fn_brainmask_in       = 'mask.nii.gz';
% ARG.fn_out                = [ARG.DIROUT 'NORDIC_VM_Monly_'];
% ARG.save_gfactor_map      = 1; %saves gfactor to nifti
% ARG.save_add_info         = 1; % save matlab ARG struct to ARG.mat containing all setting
% ARG.gfactorcorr           = 1; % perform gfactor estimation
% ARG.mask                  = 1; % = 1 the user provides a 3D brain mask. = 0 no mask is provided and VM_NORDIC uses a empty mask
% ARG.speed_factor          = 5; % integer between 1 and 10. Speeds up the denoising
% ARG.VMNORDIC              = 1; % = 1 activate VM_NORDIC, is 0 then use standard NORDIC
% ARG
% ARG.magnitude_only        = 1; % use only magnitude data
% ARG.timepoints            = []; % array containing timepoint indices on which the VMNORDIC similarity metric should be computed on (ie using only the fMRI baseline period as to avoid periods with a strong stimulus)
% [MAGN,PHASE,GFACTOR,ARG]	= NIFTI_NORDIC_VM(fn_magn_in,fn_phase_in,fn_out,ARG,fn_brainmask_in);
%
%  file_input assumes 4D NIFTIdata
%
%  OPTIONS
%  ARG.DIROUT               val =   string    Default is empty: output folder
%  ARG.noise_volume_last    val = 0 1  noise volume assumed at the end of the timeseries series: 0 = default (ie no noise scan), 1 = noise can is last volume%
%  ARG.factor_error         val = num  >1 use higher noisefloor <1 use lower noisefloor, 1 default
%  ARG.VMNORDIC             val = [0 1] 1, use VM NORDIC, 0 use standard NORDIC ( SVT threshold determined by MP estimation, yielding the gfactor, or use sigma from noise scan (ARG.noise_volume_last=1)
%  ARG.MP                   val = [0 1] 1, MPPCA denoising, using Marchencko-Pastur SVT threshold and noise estimation per patch, in a 7 x 7 x 7 kernel with patch_overlap=1
%  ARG.kernel_size_gfactor  val = [val1 val2 val], defautl is [14 14 1]
%  ARG.kernel_size_PCA      val = [val1 val2 val], default is val1=val2=val3; ratio of 11:1 between spatial and temproal voxels
%  ARG.magnitude_only       val = [0 1] using complex or magnitude only. Default is 0 using complex data
%  ARG.save_add_info        val = [0 1] 1, an additonal matlab file is being saved with degrees removed etc: default is 0
%  ARG.save_gfactor_map     val = [1 2] 1, saves the RELATIVE gfactor, 2 saves the gfactor and does not complete the NORDIC processing%
%  FOR VM NORDIC:           only used works when ARG.VMNORDIC = 1
%    ARG.speed_factor       val = positive integer (suggested between 1 and 10). Skip val reference pixels to speed up denoising via pixel matching, default is 1, 5 is a good number
%    ARG.mask               val = [0 1] 1, the user provides a 3D binary brain mask (NIFTI as last input; 0, VMNORDIC uses an mask with 1s of entire FOV
%    ARG.timepoints         val = [1:val1, or val1:val2], array containing timepoint indices on which the VMNORDIC similarity metric should be computed on (ie using only the fMRI baseline period as to avoid periods with a strong stimulus)
if nargin == 4
  fn_brainmask_in = []; % will then use mask with 1s for entire FOV in VMNORDIC
end

if ~isfield(ARG,'noise_volume_last')
  ARG.noise_volume_last = 0; % 0 there is no noise volume, 1 if last volume is noise volume
end

if ~isfield(ARG,'factor_error')
  ARG.factor_error = 1.0; % error in gfactor estimatetion. >1 use higher noisefloor <1 use lower noisefloor
end

if ~isfield(ARG,'MP')
  ARG.MP = 0; % threshold based on Marchencko-Pastur noise estimation
end

if ~isfield(ARG,'VMNORDIC')
  ARG.VMNORDIC = 0;
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

if ~isfield(ARG,'magnitude_only') % if complex data
  ARG.magnitude_only = 0; %
end

if ~isfield(ARG,'save_gfactor_map') % save out a map of a relative gfactor
  ARG.save_gfactor_map = []; %
end

if ~isfield(ARG,'mask')
  ARG.mask = 0; %
end

if ~isfield(ARG,'speed_factor')
  ARG.speed_factor = 1;
end

if ~isfield(ARG,'timepoints')
ARG.timepoints     = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ARG.magnitude_only == 0 % complex data input  
  try
    info_phase = niftiinfo(fn_phase_in);
    info_magn = niftiinfo(fn_magn_in);
  catch
    warning('The niftiinfo fails at reading the header');
    return
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
  fprintf('Phase data range is %.2f to %.2f\n', min(I_P(:)), max(I_P(:)))
  I_P = (I_P./range_norm - range_center)*2*pi;  
  
  IMG = I_M .* exp(1i*I_P); % construct complex data input IMG  
  
else % magnitude only input data
  try
    info_magn = niftiinfo(fn_magn_in);
  catch
    warning('The niftiinfo fails at reading the header');
    return
  end  
  I_M = abs(single(niftiread(fn_magn_in)));  
  info_magn.Datatype = class(I_M);
  
  IMG = single(I_M);
end

DIMS_IMG = size(IMG);
ARG.DIMS_IMG = DIMS_IMG; 
% JCWS option to skip or perform gfactor estimation
if ARG.gfactorcorr == 1
  
  if isempty(ARG.kernel_size_gfactor) || size(ARG.kernel_size_gfactor,2) < 4
    IMG_for_gfactor = IMG(:,:,:,1:min(90,end)); % takes first 90 volumes or less for gfactors estimation, but should be at least 30 volumes
  else % fourth dimension in ARG.kernel_size_gfactor is the number of dynamics the gfactor estimation should consider
    IMG_for_gfactor = IMG(:,:,:,1:min(ARG.kernel_size_gfactor(4),end));
  end
  DIMS_IMG_for_gfactor = size(IMG_for_gfactor);
  
  IMG_for_gfactor(isnan(IMG_for_gfactor)) = 0;
  IMG_for_gfactor(isinf(IMG_for_gfactor)) = 0;
  
  if isempty(ARG.kernel_size_gfactor)
    ARG.kernel_size_gfactor = [14 14 1];
  end
  QQ.IMG_processed = zeros(1,DIMS_IMG_for_gfactor(1) - ARG.kernel_size_gfactor(1));
  ARG.patch_overlap = ARG.gfactor_patch_overlap;
  
  ARG.LLR_scale = 0;
  ARG.SVT_lambda = 1;
  ARG.soft_thrs = 10; % MPPCa (when Noise varies), for noise sigma and thus gfactor estimation,

  QQ.IMG_processed = zeros(1,DIMS_IMG_for_gfactor(1) - ARG.kernel_size_gfactor(1));
  
  IMG_processed = zeros(size(QQ.IMG_processed));
  for nw1 = 2:max(1,floor(ARG.kernel_size_gfactor(1)/ARG.patch_overlap))
      IMG_processed(1,nw1 : max(1,floor(ARG.kernel_size_gfactor(1)/ARG.patch_overlap)):end) = 2;
  end
  IMG_processed(end) = 0;
  QQ.IMG_processed = IMG_processed;
    
  ARG.kernel_size_gfactor_total = prod(ARG.kernel_size_gfactor);
  disp(['Patch size for gfactor estimation = ' num2str(ARG.kernel_size_gfactor) ' = ' num2str(ARG.kernel_size_gfactor_total)]);
  
  ARG.kernel_size = ARG.kernel_size_gfactor; % main function sub_LLR_Processing needs a ARG.kernel_size variable
  disp('estimating g-factor ...')  
  for n1 = 1:size(QQ.IMG_processed,2)% nl is the index that moves the patch around! for matrixsize in X = Nx Nx - patchsize steps
    [~,~,IMG_weight,NOISE,~,~,~] = sub_LLR_Processing(IMG_for_gfactor, ARG, n1, QQ);
  end
  disp('completed estimating g-factor')  
  
  ARG.NOISE = sqrt(NOISE./IMG_weight);  
  ARG.gfactor = ARG.NOISE;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (ARG.save_gfactor_map == 1) || (ARG.save_gfactor_map == 2)    
    gfactor_IMG = abs(ARG.gfactor); 
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
  ARG.gfactor = ones(DIMS_IMG(1:3));
  GFACTOR = ARG.gfactor;  
end

% gfactor normalization, this also makes the noise have standard deviation of 1 ! therefor eigenvalue threshold estimation from random matrices with zero mean, and 1 std normal distributed noise canbe used!
 %IMG = bsxfun(@rdivide, IMG, ARG.gfactor);

IMG(isnan(IMG)) = 0;% cleanup
IMG(isinf(IMG)) = 0;% cleanup

if ARG.noise_volume_last == 1
  IMG_NOISE = IMG(:,:,:,end); % extract last volume as NOISE volume (ARG.noise_volume_last = 1 then)
  ARG.measured_noise = std(IMG_NOISE(IMG_NOISE ~= 0 )); % sqrt(2) for real and complex JWCS: WHAT TO DO WITH 'std' from complex noise or magnitude niuse????
else
  ARG.measured_noise = 1; % IF COMPLEX DATA  , CHECK!!!!
end

% %%% CHECK THIS JCWS %% %
if ARG.magnitude_only == 0 % for complex data
  ARG.measured_noise = ARG.measured_noise/sqrt(2); % scale noise with 1/sqrt(2), only for complex data
  disp(num2str(ARG.measured_noise))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine Patch size
if isempty(ARG.kernel_size_PCA)
  ARG.kernel_size = repmat(round((DIMS_IMG(4)*11)^(1/3)),1,3); % use spatial temporal dimension size ratio of 11:1: MxMxM=11*Ntimepoints = M = cuberoot (11xNtimepoints)
else
  ARG.kernel_size = ARG.kernel_size_PCA ;
end

if DIMS_IMG(3) <= ARG.kernel_size(3) % Number of slices is less than cubic kernel
  ARG.kernel_size = repmat(round((DIMS_IMG(4)*11/DIMS_IMG(3) )^(1/2)),1,2);
  ARG.kernel_size(3) = DIMS_IMG(3);
end
ARG.kernel_size_total = prod(ARG.kernel_size);
disp(['Patch size for standard NORDIC = ' num2str(ARG.kernel_size) ' = ' num2str(ARG.kernel_size_total)]);

QQ.IMG_processed = zeros(1,DIMS_IMG(1)-ARG.kernel_size(1));

ARG.patch_overlap = ARG.NORDIC_patch_overlap ;

if ARG.MP == 1 % using standard MPPCA with kernel 7x7x7 kernel and 7 patch overlap
  ARG.kernel_size = [7 7 7];
  ARG.patch_overlap = 7;
  ARG.soft_thrs = 10; % MPPCa (When Noise varies) SVT threshold and sigma estimation per patch
elseif ARG.MP == 0
   ARG.soft_thrs = []; % NORDIC (When noise is flat), do SVT using predetermined lambda
end
ARG.LLR_scale = 1;
   ARG.soft_thrs = 10; % NORDIC (When noise is flat), do SVT using predetermined lambda

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eigenvalue threshold SVT_lambda estimation of 10 random matrices, zero mean normal distributes
gamma=312/3375
lambda_min=sqrt((1+sqrt(gamma))^2*3375)
lambda_plus=sqrt((1-sqrt(gamma))^2*3375)
ARG.SVT_lambda = 0;
Niter=10;
for ntmp = 1:Niter
  [~,S,~] = svd(randn(prod(ARG.kernel_size),DIMS_IMG(4)));
  ARG.SVT_lambda = ARG.SVT_lambda + S(1,1)/Niter; % store Niter highest eigenvalues, and divide by Niter to obtain mean lambda
end

% %% JCWS CHECK THIS
% correction factor due complex data and any bias in g-factor estimation
if ARG.magnitude_only == 0 % complex data
  ARG.SVT_lambda = ARG.SVT_lambda*sqrt(2)*ARG.measured_noise*ARG.factor_error; % sqrt(2) due to complex data, ARG.factor_error correction for any bias in g-factor
else
  ARG.SVT_lambda = ARG.SVT_lambda*ARG.measured_noise*ARG.factor_error; % no sqrt(2) due to magnitude, ARG.factor_error correction for any bias in g-factor
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ARG.VMNORDIC == 0 % standard NORDIC
  
  QQ.IMG_processed = zeros(1,DIMS_IMG(1) - ARG.kernel_size(1));
  
  IMG_processed = zeros(size(QQ.IMG_processed));
  for nw1 = 2:max(1,floor(ARG.kernel_size(1)/ARG.patch_overlap))
    IMG_processed(1,nw1 : max(1,floor(ARG.kernel_size(1)/ARG.patch_overlap)):end) = 2;
  end
  IMG_processed(end) = 0; 
  QQ.IMG_processed = IMG_processed;
  
  disp('starting NORDIC ...')  
  for n1 = 1:size(QQ.IMG_processed,2)
    [IMG_denoised,~,IMG_weight,NOISE,Component_threshold,energy_removed,SNR_weight] = sub_LLR_Processing(IMG, ARG, n1, QQ) ;
  end
  disp('completing NORDIC ...')
  
  IMG_denoised = IMG_denoised./repmat((IMG_weight),[1 1 1 DIMS_IMG(4)]); % Assumes that the combination is with N instead of sqrt(N). Works for NVR not MPPCA
  
  ARG.NOISE = sqrt(NOISE./IMG_weight);
  ARG.Component_threshold = Component_threshold./IMG_weight;
  ARG.energy_removed = energy_removed./IMG_weight;
  ARG.SNR_weight = SNR_weight./IMG_weight;
  
else % VM-NORDIC  
  disp('starting non-local voxel matching VM-NORDIC...')  
  [IMG_denoised,ARG] = voxel_matchNLLR(IMG, ARG);  
  disp('completing NORDIC with non-local voxel matching (VM_NORDIC)...')
  
  % write visiting map
  niftiwrite(single(ARG.visitingmap),[fn_out 'visitingmap.nii'],'Compressed',true)
end

% gfactor renormalization
%IMG_denoised = bsxfun(@times, IMG_denoised, ARG.gfactor);

IMG_denoised(isnan(IMG_denoised)) = 0;

if ARG.magnitude_only == 0 % complex data
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
  
elseif ARG.magnitude_only == 1
  
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
end

function [IMG_denoised,IMG,IMG_weight,NOISE, Component_threshold,energy_removed,SNR_weight] = sub_LLR_Processing(IMG, ARG, n1, QQ)
DIMS = size(IMG);
IMG_denoised = zeros(DIMS); 
IMG_weight = zeros(DIMS(1:3));
NOISE = zeros(DIMS(1:3));
Component_threshold = zeros(DIMS(1:3));
energy_removed = zeros(DIMS(1:3));
SNR_weight = zeros(DIMS(1:3));  

Iwindowx = 1:ARG.kernel_size(1)+(n1-1); %patch selection window indices range in X direction

if QQ.IMG_processed(1,n1) ~= 1 && QQ.IMG_processed(1,n1) ~= 3 % not being processed also not completed yet
  if QQ.IMG_processed(1,n1) == 2 % processed but not added.
     % %% CHECK ALL THIS JCWS !!
    try % try to load otherwise go to next slice
     load([ARG.filename 'slice' num2str(n1) '.mat'],'DATA_full2')
    catch
     QQ.IMG_processed(1,n1)=0; % identified as bad file and being identified for reprocessing
     return
    end
    QQ.IMG_processed(1,n1) = 0; %     
  end  
  if QQ.IMG_processed(1,n1) ~= 2    
    QQ.IMG_processed(1,n1) = 1 ; % STARTING
    IMG_block = IMG(Iwindowx,:,:,:); % select kernelsize(1), IMG_block is a patchblock confined to kernelsize X and crossing all Y and Z, timepoints range
    lambda = ARG.LLR_scale*ARG.SVT_lambda;
    disp(['lambla threshold = ' num2str(lambda)])    
    
    IMG_weight_patchtmp = IMG_weight(Iwindowx,:,:,:);
    NOISE_patchtmp = NOISE(Iwindowx,:,:,:);
    Component_threshold_patchtmp = Component_threshold(Iwindowx,:,:,:);
    energy_removed_patchtmp = energy_removed(Iwindowx,:,:,:);
    SNR_weight_tmp = SNR_weight(Iwindowx,:,:,:);
    
    patch_avg = 1;
    [IMG_block_denoised,IMG_weight_patchtmp, NOISE_patchtmp, Component_threshold_patchtmp, energy_removed_tmp, SNR_weight_tmp] = ...
      subfunction_loop_for_NVR_avg_update(IMG_block,ARG.kernel_size(2),ARG.kernel_size(3),lambda,patch_avg,ARG.soft_thrs,IMG_weight_patchtmp,ARG,NOISE_patchtmp,Component_threshold_patchtmp,energy_removed_patchtmp,SNR_weight_tmp);
    
    IMG_weight(Iwindowx,:,:,:) = IMG_weight_patchtmp;
    NOISE(Iwindowx,:,:,:) = NOISE_patchtmp;
    Component_threshold(Iwindowx,:,:,:) = Component_threshold_patchtmp;
    energy_removed(Iwindowx,:,:,:) = energy_removed_tmp;
    SNR_weight(Iwindowx,:,:,:) = SNR_weight_tmp;
        
  end  
  IMG_denoised(Iwindowx,1:size(IMG_block_denoised,2),:,:) = IMG_denoised(Iwindowx,1:size(IMG_block_denoised,2),:,:) + IMG_block_denoised;
  QQ.IMG_processed(1,n1) = 3 ;
end
end

function [IMG_block_denoised, IMG_weight,NOISE,IMG_tmp_update_threshold,energy_removed,SNR_weight] = subfunction_loop_for_NVR_avg_update(IMG_block,kernel_sizeY,kernel_sizeZ,lambda,patch_avg, soft_thrs,IMG_weight,ARG,NOISE,IMG_tmp_update_threshold,energy_removed,SNR_weight)
DIMS_block = size(IMG_block); % IMG_block is a patchblock confined to kernelsize X and crossing all Y and Z, timepoints range
IMG_block_denoised = zeros(DIMS_block);
patch_scale = 1;

Ipatchstepy = max(1,floor(kernel_sizeY/ARG.patch_overlap));
Ipatchmaxy = DIMS_block(2) - kernel_sizeY + 1;

Ipatchstepz = max(1,floor(kernel_sizeZ/ARG.patch_overlap));
Ipatchmaxz = DIMS_block(3) - kernel_sizeZ + 1;

for n2 = [1 : Ipatchstepy : Ipatchmaxy, Ipatchmaxy] % move patch around in Y direction in steps of kernelsize y (kernel_sizeY)
  for n3 = [1 : Ipatchstepz : Ipatchmaxz, Ipatchmaxz] % move patch around in Z direction in steps of kernelsize z (w3)    
    Iwindowy = (1:kernel_sizeY)+(n2-1); % moving kernel range and index in y
    Iwindowz = (1:kernel_sizeZ)+(n3-1); % moving kernel range and index in z
    
    % this is the patch: dimensions kernelsize(1) x kernelsize(2) x kernelsize(3) x timepoints
    IMG_patch = IMG_block(:, Iwindowy, Iwindowz,:);
    IMG_patch_casorati = reshape(IMG_patch,[],DIMS_block(4)); % this is the Casorati matrix of the patch
    
    [U,S,V] = svd(IMG_patch_casorati,'econ'); % this can go faster, see voxel_matchNLLR.m code using 'eig.m' function
    S = diag(S);    
    [idx] = sum(S < lambda);
    
    if isempty(soft_thrs)% NORDIC denoising,  SVT thresholding eigenvalues      
      energy_scrub = sqrt(sum(S)).\sqrt(sum(S(S < lambda)));
      S(S < lambda) = 0; %SVT
      t = idx;
    elseif soft_thrs ~= 10 % use soft thresholding on the eigenvalues using soft_thrs*labmda 
      S(S < lambda*soft_thrs) = 0; %SVT
      energy_scrub = 0;
      t = 1;
    elseif soft_thrs == 10 % using Marchencko-Pastur sigma estimation 
      centering = 0;
      MM = size(IMG_patch_casorati,1); % Correction for some zero entries
      if MM > 0
        NNN = size(IMG_patch_casorati,2);
        R = min(MM, NNN);
        scaling = (max(MM, NNN) - (0:R-centering-1)) / NNN;
        scaling = scaling(:);
        eigenvals = S;
        eigenvals = (eigenvals).^2 / NNN;
        % First estimation of Sigma^2; Eq 1 from ISMRM presentation
        csum = cumsum(eigenvals(R-centering:-1:1));
        cmean = csum(R-centering:-1:1)./(R-centering:-1:1)';
        sigmasq_1 = cmean./scaling;
        % Second estimation of Sigma^2; Eq 2 from ISMRM presentation
        gamma = (MM - (0:R-centering-1)) / NNN;
        rangeMP = 4*sqrt(gamma(:));
        rangeData = eigenvals(1:R-centering) - eigenvals(R-centering);
        sigmasq_2 = rangeData./rangeMP;
        t = find(sigmasq_2 < sigmasq_1, 1);
        idx = size(S(t:end),1) ;
        energy_scrub = sqrt(sum(S)).\sqrt(sum(S(t:end)));
        S(t:end) = 0; % SVT
      else % all zero entries
        t = 1;
        energy_scrub = 0;
        sigmasq_2 = 0;
      end    
    end
    
    IMG_patch_denoised = U*diag(S)*V';
    IMG_patch_denoised = reshape(IMG_patch_denoised,size(IMG_patch)); % reshape denoised by SVT IMG_patch_denoised back to a 4D patch
    
    if patch_scale == 1
    else
      patch_scale = size(S,1) - idx;
    end
    
    if isempty(t)
      t = 1;
    end % threshold removed all.
    
    if patch_avg == 1
      IMG_block_denoised(:,Iwindowy, Iwindowz,:) = IMG_block_denoised(:,Iwindowy, Iwindowz,:) + patch_scale*IMG_patch_denoised;
      
      IMG_weight(:,Iwindowy, Iwindowz,:) = IMG_weight(:,Iwindowy, Iwindowz,:) + patch_scale;
      IMG_tmp_update_threshold(:,Iwindowy, Iwindowz,:) = IMG_tmp_update_threshold(:,Iwindowy, Iwindowz,:) + idx;
      energy_removed(:,Iwindowy, Iwindowz,:) = energy_removed(:,Iwindowy, Iwindowz,:) + energy_scrub;
      SNR_weight(:,Iwindowy, Iwindowz,1) = SNR_weight(:,Iwindowy, Iwindowz,1) + S(1)./S(max(1,t-1));
      if soft_thrs == 10
        NOISE(1:DIMS_block(1),Iwindowy, Iwindowz,1) = NOISE(1:DIMS_block(1),Iwindowy, Iwindowz,1) + sigmasq_2(t);
      end
    end
  end
end
end

