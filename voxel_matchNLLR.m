function [KSP2a_RECON,ARG] = voxel_matchNLLR(KSP2a, ARG, mask) %make mask = 1s out of the function if not provided

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MASK

noisy_masked  = zeros(size(KSP2a));
imdim         = size(KSP2a);

if ARG.mask == 1
    mask(mask==0) = NaN;    
    noisy_masked = KSP2a.*mask;
elseif ARG.mask == 0 % no user supplied mask, greate empty mask of size data  
    mask          = ones(imdim(1:3));        
    noisy_masked = KSP2a.*mask;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE CHUNKS

imdim         = size(noisy_masked);
Nvoxels_brain = length(find(~isnan(noisy_masked(:))));
max_voxels    = 5e6;
step          = round(Nvoxels_brain/max_voxels);  % step between slices per chunk
chunk_size    = floor(size(noisy_masked,3)/step); % amount of slices per chunk

% ensure chunks are of similar sizes
if chunk_size > round(imdim(3)*0.6) && chunk_size < imdim(3)
    while chunk_size > round(imdim(3)*0.6)
        max_voxels = max_voxels - 5e5;
        step       = round(Nvoxels_brain/max_voxels);
        chunk_size = floor(size(noisy_masked,3)/step);
    end
end

len    = zeros(step,1);
slices = zeros(step,chunk_size);
images = cell(1,step);

% make chunks
for s = 1:step
    img         = noisy_masked(:,:,s:step:end-(step-s),:);
    slices(s,:) = s:step:size(noisy_masked,3)-(step-s);
    Nvoxels_img = reshape(img(:,:,:,1),prod(imdim(1:2))*chunk_size,1);
    len(s,:)    = length(Nvoxels_img(~isnan(Nvoxels_img(:))));
    images{1,s} = img;
end

% if some slices are left out, reinclude them in smaller chunks
if imdim(3)/step ~= round(imdim(3)/step)
    left          = abs(chunk_size*step - imdim(3));
    [~,short_idx] = mink(len,left);
    missing       = find(~ismember(1:imdim(3),slices(:)));
    for i = 1:left
        s_idx              = short_idx(i);
        short              = images{s_idx};
        short(:,:,end+1,:) = noisy_masked(:,:,missing(end-(i-1)),:);
        images{1,s_idx}    = short;
        slices(short_idx(i),chunk_size+1) = missing(end-(i-1));
    end
end

% release some memory
clear noisy_masked KSP2a mask short img

KSP2a_RECON = zeros(imdim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS INDIVIDUAL CHUNKS

for s = 1:step % denoise one chunk per time
    
    chunk  = images{s};
    matdim = size(chunk);
    
    timeseries = reshape(chunk,prod(matdim(1:3)),matdim(4));
    brain_vox  = timeseries(~isnan(timeseries(:,1)),:);
    brain_idx  = find(~isnan(timeseries(:,1)));
    
    % timepoints to use for similarity estimation
    if isempty(ARG.timepoints)
        tp = logical(ones(1,matdim(4)));

    else
        tp = ismember(1:matdim(4),ARG.timepoints);     
    end
    
    if ARG.magnitude_only == 0
        tp    = repmat(tp,1,2);
        brain = cat(2,real(brain_vox),imag(brain_vox));
    else
        brain = brain_vox;
    end
    
    references = brain(1:ARG.speed_factor:end,:);
    
    nTimepoints_forsimilarity = length(tp);
    if s == 1
        if nTimepoints_forsimilarity > 40
            patch_size1 = nTimepoints_forsimilarity*3; % initial (overestimated) guess for patch size
        else                         % correction for datasets with short timeseries
            patch_size1 = ceil(nTimepoints_forsimilarity*(((length(brain)/(round(length(brain),1,'significant')*100))*nTimepoints_forsimilarity)^(-2.7)));
        end
        
        disp(['Estimating similarities chunk ' num2str(s) ' of ' num2str(step)]);
        distances         = knnsearch(abs(brain(:,tp)),abs(references(:,tp)),'K',patch_size1,'Distance','euclidean');
        similars_idx_full = brain_idx(distances);
        patch_size2       = tune_patchsize(patch_size1,timeseries,similars_idx_full,matdim,ARG); % finetune patch size
        similars_idx_cut  = similars_idx_full(:,1:patch_size2);
    else
        disp(['Estimating similarities chunk ' num2str(s) ' of ' num2str(step)]);
        distances         = knnsearch(abs(brain(:,tp)),abs(references(:,tp)),'K',patch_size2,'Distance','euclidean'); %Indexes in VEC
        similars_idx_cut  = brain_idx(distances); %Indexes in complete image\
    end
    disp(['Tuned patchsize = ' num2str(patch_size2)])
    if ARG.speed_factor > 1 % reinclude omitted voxels due to speeding up
        miss = ismember(brain_idx,similars_idx_cut);
        zrs  = find(miss==0);
        disp([num2str((length(zrs)/size(brain,1))*100) '% of voxels to re-include due to acceleration'])
        
        brain_iz_vals = timeseries(brain_idx(zrs),:);
        
        if ARG.magnitude_only == 0
            brain_iz = cat(2,real(brain_iz_vals),imag(brain_iz_vals));
        else
            brain_iz = brain_iz_vals;
        end
        
        distances_iz    = knnsearch(abs(brain(:,tp)),abs(brain_iz(:,tp)),'K',patch_size2,'Distance','euclidean'); %Indexes in VEC
        similars_iz     = brain_idx(distances_iz);
        similars_all    = cat(1,similars_idx_cut,similars_iz);
        
        if isfield(ARG,'save_add_info')
            re_included                    = timeseries(:,1);
            re_included(brain_idx(zrs(:))) = mean(abs(re_included(brain_idx(zrs(:)))),'omitnan')*5;
            ARG.re_included                = reshape(re_included,size(chunk,1,2,3));
        end
    else
        similars_all = similars_idx_cut;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % THRESHOLD
    
    thr = 0;
    for i = 1:10
        if patch_size2 < matdim(4)
            A = randn(patch_size2,matdim(4));
        else
            A = randn(matdim(4),patch_size2);
        end
        E   = A*A';
        St  = eig(E);
        thr = thr + St(end,end);
    end
    thr = thr/10;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SVT
    
    patches = reshape(timeseries(similars_all,:),[],patch_size2,matdim(4)); % all patches
    img_vec = zeros(size(timeseries));
    wei_vec = zeros(size(timeseries,1),1);
    
    if patch_size2 < matdim(4) % transpose patch if necessary, just for code speed
        for i = 1:size(patches,1)
            
            Y     = squeeze(patches(i,:,:));
            M     = size(Y,1);
            E     = Y*Y';             % covariance matrix
            [U,S] = eig(E);
            N     = length(S(S>thr)); % SVT
            UL    = U(:,M-N+1:end);
            YL    = UL*(UL'*Y);       % low-rank representation of Y
            
            img_vec(similars_all(i,1:patch_size2),:) = img_vec(similars_all(i,1:patch_size2),:)+YL;
            wei_vec(similars_all(i,1:patch_size2),:) = wei_vec(similars_all(i,1:patch_size2),:)+ones(patch_size2, 1);
        end
    else
        for i = 1:size(patches,1)
            
            Y     = (squeeze(patches(i,:,:)))';
            M     = size(Y,1);
            E     = Y*Y';             % covariance matrix
            [U,S] = eig(E);
            N     = length(S(S>thr)); % SVT
            UL    = U(:,M-N+1:end);
            YL    = UL*(UL'*Y);       % low-rank representation of Y
            
            img_vec(similars_all(i,1:patch_size2),:) = img_vec(similars_all(i,1:patch_size2),:)+YL';
            wei_vec(similars_all(i,1:patch_size2),:) = wei_vec(similars_all(i,1:patch_size2),:)+ones(patch_size2, 1);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RECONSTRUCT DENOISED IMAGE
    
    img_out   = reshape(img_vec,size(chunk));
    wei_out   = repmat(reshape(wei_vec,matdim(1:3)),[1 1 1 matdim(4)]);
    img_final = img_out./wei_out;
    
    % release some more memory
    clear img_out wei_out distances similars_all similars_idx_full ...
        similars_idx_cut similars_iz img_vec wei_vec patches;
    
    position = nonzeros(slices(s,:));
    KSP2a_RECON(:,:,position,:) = img_final;
    clear img_final;
end


