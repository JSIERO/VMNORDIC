function [IMGDATA_denoised, ARG] = voxel_matchNLLR(IMGDATA, ARG) %make mask = 1s out of the function if not provided

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MASK

IMGDATA_masked  = zeros(size(IMGDATA));
DIMS_IMG         = size(IMGDATA);
ARG.visitingmap = zeros(DIMS_IMG(1:3)); % store the visiting map

if ARG.mask == 1
    mask = single(niftiread(ARG.fn_brainmask_in)); % load mask nifti
    mask(mask==0) = NaN;
    IMGDATA_masked = IMGDATA.*mask;
elseif ARG.mask == 0 % no user supplied mask, greate empty mask of size data
    mask          = ones(DIMS_IMG(1:3));
    IMGDATA_masked = IMGDATA.*mask;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE CHUNKS

Nvoxels_brain = length(find(~isnan(IMGDATA_masked(:))));
max_voxels    = 5e6;
step          = round(Nvoxels_brain/max_voxels);  % step between slices per chunk
chunk_size    = floor(size(IMGDATA_masked,3)/step); % amount of slices per chunk

% ensure chunks are of similar sizes
if chunk_size > round(DIMS_IMG(3)*0.6) && chunk_size < DIMS_IMG(3)
    while chunk_size > round(DIMS_IMG(3)*0.6)
        max_voxels = max_voxels - 5e5;
        step       = round(Nvoxels_brain/max_voxels);
        chunk_size = floor(size(IMGDATA_masked,3)/step);
    end
end

len    = zeros(step,1);
slices = zeros(step,chunk_size);
images = cell(1,step);

% make chunks
for s = 1:step
    img         = IMGDATA_masked(:,:,s:step:end-(step-s),:);
    slices(s,:) = s:step:size(IMGDATA_masked,3)-(step-s);
    Nvoxels_img = reshape(img(:,:,:,1),prod(DIMS_IMG(1:2))*chunk_size,1);
    len(s,:)    = length(Nvoxels_img(~isnan(Nvoxels_img(:))));
    images{1,s} = img;
end

% if some slices are left out, reinclude them in smaller chunks
if DIMS_IMG(3)/step ~= round(DIMS_IMG(3)/step)
    left          = abs(chunk_size*step - DIMS_IMG(3));
    [~,short_idx] = mink(len,left);
    missing       = find(~ismember(1:DIMS_IMG(3),slices(:)));
    for i = 1:left
        s_idx              = short_idx(i);
        short              = images{s_idx};
        short(:,:,end+1,:) = IMGDATA_masked(:,:,missing(end-(i-1)),:);
        images{1,s_idx}    = short;
        slices(short_idx(i),chunk_size+1) = missing(end-(i-1));
    end
end

% release some memory
clear IMGDATA_masked IMGDATA mask short img

IMGDATA_denoised = zeros(DIMS_IMG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS INDIVIDUAL CHUNKS

for s = 1:step % denoise one chunk per time
    
    chunk  = images{s};
    matdim = size(chunk);
    disp(['size chunk(' num2str(s) ')  = ' num2str(matdim)]);
    chunk_all_timeseries = reshape(chunk,prod(matdim(1:3)),matdim(4));
    chunk_brain_timeseries  = chunk_all_timeseries(~isnan(chunk_all_timeseries(:,1)),:); % select brain (non NaN) signal chunk_all_timeseries from chunk
    Idx_chunk_brain_timeseries  = find(~isnan(chunk_all_timeseries(:,1))); % store linear indices of the brain (non NaN) signal chunk_all_timeseries from chunk
    
    % timepoints to use for similarity estimation
    if isempty(ARG.timepoints)
        tp = true(1,matdim(4));
        
    else
        tp = ismember(1:matdim(4),ARG.timepoints);
    end
    
    if ARG.magnitude_only == 0 % concatenate real and imaginary data in 1 vector, i.e. double the timepoints)
        tp    = repmat(tp,1,2);
        chunk_brain_timeseries = cat(2,real(chunk_brain_timeseries),imag(chunk_brain_timeseries));
    end
    
    % select reference voxels in chunk, jumping every ARG.speed_factor for speed
    chunk_brain_timeseries_refs = chunk_brain_timeseries(1:ARG.speed_factor:end,:);
    
    nTimepoints_forsimilarity = length(tp);
    
    if s == 1
        if nTimepoints_forsimilarity >= 40 || nTimepoints_forsimilarity <= 500
            patchsize_initial = nTimepoints_forsimilarity*11; % initial large guess for patch size, using the NORDIC M(spatial) x Q(time) = 11 ratio
        elseif nTimepoints_forsimilarity > 500
            patchsize_initial = min(nTimepoints_forsimilarity*11, 11^3); % initial guess for patch size for very long timeseries (max a 11^3 spatial patchsize is used as in NORDIC)
        elseif nTimepoints_forsimilarity < 40 % patch_size_initial  correction for datasets with small amount of timepoints < 40
            patchsize_initial = ceil(nTimepoints_forsimilarity*(((size(chunk_brain_timeseries,1))/(round(size(chunk_brain_timeseries,1),1,'significant')*100))*nTimepoints_forsimilarity)^(-2.7));
        end
        
        ARG.patchsize_initial_VMNORDIC = patchsize_initial;
        disp(['Initial patchsize = ' num2str(patchsize_initial) ' = ' num2str(round(patchsize_initial^(1/3),1)) ' ^(1/3)'])
        
        disp(['Estimating similarities chunk ' num2str(s) ' of ' num2str(step)]); % make sure we work on magnitude data to bade the similarity one: taking abs()
        
        % matrix of indices of neighbours (amount is the patch_size_initial) for each reference point (length array references)
        Idx_similar_neighbours = knnsearch(abs(chunk_brain_timeseries(:,tp)), abs(chunk_brain_timeseries_refs(:,tp)), 'K', patchsize_initial, 'Distance', 'euclidean');
        Idx_similar_full = Idx_chunk_brain_timeseries(Idx_similar_neighbours); % find the corresponding indices of the chunk dataset
        
        % finetune patch size
        patchsize_tuned  = tune_patchsize(patchsize_initial, chunk_all_timeseries, Idx_similar_full, matdim, ARG);
        Idx_similar_patchsize_tuned  = Idx_similar_full(:,1:patchsize_tuned);  % select only the first patchsize tuned amount of similar voxels
        
        ARG.patch_size_tuned_VMNORDIC = patchsize_tuned;
        disp(['Tuned patchsize = ' num2str(patchsize_tuned) ' = ' num2str(round(patchsize_tuned^(1/3),1)) ' ^(1/3)'])
        
    elseif s > 1
        disp(['Estimating similarities chunk ' num2str(s) ' of ' num2str(step)]);
        Idx_similar_neighbours       = knnsearch(abs(chunk_brain_timeseries(:,tp)),abs(chunk_brain_timeseries_refs(:,tp)),'K',patchsize_tuned,'Distance','euclidean'); %Indexes in VEC
        Idx_similar_patchsize_tuned  = Idx_chunk_brain_timeseries(Idx_similar_neighbours); %Indexes in complete image\
    end
    
    
    if ARG.speed_factor > 1 % reinclude omitted voxels due to speeding up
        missed = ismember(Idx_chunk_brain_timeseries,Idx_similar_patchsize_tuned);
        Idx_missed  = find(missed==0);
        disp([num2str((length(Idx_missed)/size(chunk_brain_timeseries,1))*100) '% of voxels to re-include due to acceleration'])
        
        chunk_missed_timeseries = chunk_all_timeseries(Idx_chunk_brain_timeseries(Idx_missed),:);
        
        if ARG.magnitude_only == 0  % concatenate real and imaginary data in 1 vector, i.e. double the timepoints)
            chunk_missed_timeseries = cat(2,real(chunk_missed_timeseries),imag(chunk_missed_timeseries));
        end
        
        % matrix of indices of neighbours (amount is the patch_size_initial) for each reference point (length array references)
        Idx_similar_neighbours_iz      = knnsearch(abs(chunk_brain_timeseries(:,tp)),abs(chunk_missed_timeseries(:,tp)),'K',patchsize_tuned,'Distance','euclidean'); %Indexes in VEC
        Idx_similar_neighbours_full_iz = Idx_chunk_brain_timeseries(Idx_similar_neighbours_iz);  % find the corresponding indices of the chunk dataset
        if size(Idx_similar_neighbours_full_iz,2) == 1 % transpose indices array, as when only 1 voxel is found the above line flips the output
            Idx_similar_neighbours_full_iz = Idx_similar_neighbours_full_iz';
        end
        Idx_similar_neighbours_all     = cat(1,Idx_similar_patchsize_tuned,Idx_similar_neighbours_full_iz);
        
        if isfield(ARG,'save_add_info')
            re_included                    = chunk_all_timeseries(:,1);
            re_included(Idx_chunk_brain_timeseries(Idx_missed(:))) = mean(abs(re_included(Idx_chunk_brain_timeseries(Idx_missed(:)))),'omitnan')*5;
            ARG.re_included                = reshape(re_included,size(chunk,1,2,3));
        end
    else
        Idx_similar_neighbours_all = Idx_similar_patchsize_tuned;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % THRESHOLD
    
    thr = 0;
    for i = 1:10
        if patchsize_tuned < matdim(4)
            A = randn(patchsize_tuned,matdim(4));  % generate 10-times random matrices with zero-mean normal distributed noise
        else
            A = randn(matdim(4),patchsize_tuned);
        end
        E   = A*A';
        St  = eig(E);  % perform  eigenvalue decomposition
        thr = thr + St(end,end);
    end
    thr = thr/10;    % take the average of thresholds
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SVT
    
    patches = reshape(chunk_all_timeseries(Idx_similar_neighbours_all,:),[],patchsize_tuned,matdim(4)); % all patches
    img_vec = zeros(size(chunk_all_timeseries));
    wei_vec = zeros(size(chunk_all_timeseries,1),1);
    
    if patchsize_tuned < matdim(4) % transpose patch if necessary, just for code speed
        for i = 1:size(patches,1)
            
            Y     = squeeze(patches(i,:,:));
            M     = size(Y,1);
            E     = Y*Y';             % covariance matrix
            [U,S] = eig(E);           % perform  eigenvalue decomposition
            N     = length(S(S>thr)); % SVT
            UL    = U(:,M-N+1:end);
            YL    = UL*(UL'*Y);       % low-rank representation of Y
            
            img_vec(Idx_similar_neighbours_all(i,1:patchsize_tuned),:) = img_vec(Idx_similar_neighbours_all(i,1:patchsize_tuned),:)+YL;
            wei_vec(Idx_similar_neighbours_all(i,1:patchsize_tuned),:) = wei_vec(Idx_similar_neighbours_all(i,1:patchsize_tuned),:)+ones(patchsize_tuned, 1);
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
            
            img_vec(Idx_similar_neighbours_all(i,1:patchsize_tuned),:) = img_vec(Idx_similar_neighbours_all(i,1:patchsize_tuned),:)+YL';
            wei_vec(Idx_similar_neighbours_all(i,1:patchsize_tuned),:) = wei_vec(Idx_similar_neighbours_all(i,1:patchsize_tuned),:)+ones(patchsize_tuned, 1);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RECONSTRUCT DENOISED IMAGE
    
    img_out   = reshape(img_vec,size(chunk));
    wei_out   = repmat(reshape(wei_vec,matdim(1:3)),[1 1 1 matdim(4)]);
    img_final = img_out./wei_out;
    
    position = nonzeros(slices(s,:));
    IMGDATA_denoised(:,:,position,:) = img_final;
    ARG.visitingmap(:,:,position) = squeeze(wei_out(:,:,:,1));
    % release some more memory
    clear img_out wei_out distances Idx_similar_neighbours_all similars_idx_full ...
        similars_idx_cut similars_iz img_vec wei_vec patches img_final
    
end


