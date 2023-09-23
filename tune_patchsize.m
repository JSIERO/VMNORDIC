function [voxel_number] = tune_patchsize(patchsize_initial, timeseries, Idx_similar_full, matdim,ARG)
disp('Finetuning patch size...');
%finding optimal patch size based maximum tSNR and patchsize_initial
range  = (0.05:0.05:1);         % fractions of the patchsize_initial to try out
fract  = 1/50*ARG.speed_factor; % fraction of data to use for tests 
iter   = 5;
PNsT   = zeros(1,1);
count0 = 1;

for f = 1:iter
    disp(['test ' num2str(count0) ' of ' num2str(iter)]);
    
    [TSNRs,PNs]    = finetune(patchsize_initial,timeseries,Idx_similar_full,range,fract,matdim);
    
    [~,T]          = max(TSNRs);    
    PNsT(count0,:) = PNs(T);
    count0         = count0+1;
end

voxel_number = max(PNsT); % take highest patch size as a small sample leads to underestimation of optimal patch size
return

function [TSNRs,PNs] = finetune(patchsize_initial,timeseries,Idx_similar_full,range,fract,matdim)

TSNRs = zeros(1,1);
PNs   = TSNRs;
count = 1;
for tt = 1:length(range)    
    
    PN = round(patchsize_initial*range(tt)); % test patch number   
    Idx_sample = datasample(Idx_similar_full,round(size(Idx_similar_full,1)*fract),1); % take a random sample of the data --> change this!
    Idx_sample = Idx_sample(:,1:PN);
    
    thr   = 0;
    for i = 1:10
        if PN < matdim(4)
            A   = randn(PN,matdim(4)); % generate 10-times random matrices with zero-mean normal distributed noise
        else
            A   = randn(matdim(4),PN); 
        end
        E   = A*A';
        St  = eig(E); % perform  eigenvalue decomposition
        thr = thr+St(end,end); % largest eigenvalue is the SVT threshold
    end
    thr     = thr/10; % take the average of thresholds
    
    Y_hat     = zeros(size(timeseries));
    W_hat     = zeros(size(timeseries,1),1);
    patches   = reshape(timeseries(Idx_sample,:),[], PN, matdim(4));
    
    if PN < matdim(4)
        for i = 1:size(patches,1)
            
            Y       = squeeze(patches(i,:,:));
            M       = size(Y,1);
            E       = Y*Y';
            [U,S]   = eig(E);
            N       = length(S(S>thr)); % number of components to keep ,SVT 
            UL      = U(:,M-N+1:end);
            YL      = UL*(UL'*Y); % reconstruct after SVT matrix Y
            
            Y_hat(Idx_sample(i,:),:) = Y_hat(Idx_sample(i,:),:)+YL;
            W_hat(Idx_sample(i,:),:) = W_hat(Idx_sample(i,:),:)+ones(PN, 1); % weight (visiting) vector, ie how many times a voxel has been SVTed
        end
    else % flip matrix Y when number of timepoints too large
        for i = 1:size(patches,1)    
            
            Y       = (squeeze(patches(i,:,:)))';
            M       = size(Y,1);
            E       = Y*Y';
            [U,S]   = eig(E);
            N       = length(S(S>thr));
            UL      = U(:,M-N+1:end);
            YL      = UL*(UL'*Y);
            
            Y_hat(Idx_sample(i,:),:) = Y_hat(Idx_sample(i,:),:)+YL';
            W_hat(Idx_sample(i,:),:) = W_hat(Idx_sample(i,:),:)+ones(PN, 1);
        end
    end
    
    final_samp     = abs(Y_hat./repmat(W_hat,[1 matdim(4)])); % average denoised Y hat by the weight (visiting) vector
    fin            = final_samp(~isnan(final_samp(:,1)),:);
    TSNRs(count,:) = mean(mean(fin,2)./std(fin,0,2),'omitnan'); % compute tSNR
    PNs(count,:)   = PN;
    
    if tt ~= 1 && TSNRs(end,:) < TSNRs(end-1,:)        
        return  % stop loop if maximum is found
    end   
    
    count = count+1;
    
end
return


