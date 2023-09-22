function [voxel_number] = tune_patchsize(pixnum, timeseries, similars, matdim,ARG)
disp('Finetuning patch size...');

range  = (0.05:0.05:1);         % factors to multiply to initial guess
fract  = 1/50*ARG.speed_factor; % fraction of data to use for tests 
iter   = 5;
PNsT   = zeros(1,1);
count0 = 1;

for f = 1:iter
    disp(['test ' num2str(count0) ' of ' num2str(iter)]);
    
    [TSNRs,PNs]    = finetune(pixnum,timeseries,similars,range,fract,matdim);
    
    [~,T]          = max(TSNRs);    
    PNsT(count0,:) = PNs(T);
    count0         = count0+1;
end

voxel_number = max(PNsT); % take highest patch size as a small sample leads to underestimation of optimal patch size
return

function [TSNRs,PNs] = finetune(pixnum,all_voxels,similars,range,fract,matdim)

TSNRs = zeros(1,1);
PNs   = TSNRs;
count = 1;
for tt = 1:length(range)    
    
    PN          = round(pixnum*range(tt)); % test patch number   
    near_sample = datasample(similars,round(size(similars,1)*fract),1); % take a random sample of the data --> change this!
    near_sample = near_sample(:,1:PN);
    
    thr   = 0;
    for i = 1:10
        if PN < matdim(4)
            A   = randn(PN,matdim(4));
        else
            A   = randn(matdim(4),PN);
        end
        E   = A*A';
        St  = eig(E);
        thr = thr+St(end,end);
    end
    thr     = thr/10;
    
    Y_hat     = zeros(size(all_voxels));
    W_hat     = zeros(size(all_voxels,1),1);
    patches   = reshape(all_voxels(near_sample,:),[],PN,matdim(4));
    
    if PN < matdim(4)
        for i = 1:size(patches,1)
            
            Y       = squeeze(patches(i,:,:));
            M       = size(Y,1);
            E       = Y*Y';
            [U,S]   = eig(E);
            N       = length(S(S>thr));
            UL      = U(:,M-N+1:end);
            YL      = UL*(UL'*Y);
            
            Y_hat(near_sample(i,1:PN),:) = Y_hat(near_sample(i,1:PN),:)+YL;
            W_hat(near_sample(i,1:PN),:) = W_hat(near_sample(i,1:PN),:)+ones(PN, 1);
        end
    else
        for i = 1:size(patches,1)    
            
            Y       = (squeeze(patches(i,:,:)))';
            M       = size(Y,1);
            E       = Y*Y';
            [U,S]   = eig(E);
            N       = length(S(S>thr));
            UL      = U(:,M-N+1:end);
            YL      = UL*(UL'*Y);
            
            Y_hat(near_sample(i,1:PN),:) = Y_hat(near_sample(i,1:PN),:)+YL';
            W_hat(near_sample(i,1:PN),:) = W_hat(near_sample(i,1:PN),:)+ones(PN, 1);
        end
    end
    
    final_samp     = abs(Y_hat./repmat(W_hat,[1 matdim(4)]));
    fin            = final_samp(~isnan(final_samp(:,1)),:);
    TSNRs(count,:) = mean(mean(fin,2)./std(fin,0,2),'omitnan');
    PNs(count,:)   = PN;
    
    if tt ~= 1 && TSNRs(end,:) < TSNRs(end-1,:)
        return  % stop loop if maximum is found
    end   
    
    count = count+1;
    
end
return


