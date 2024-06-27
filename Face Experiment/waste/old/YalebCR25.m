close all
clear 
addpath(genpath('C:\Users\Pastoral\Desktop\Copy_SSC_ADMM_v1.1'));
addpath(genpath('C:\Users\Pastoral\Desktop\Unlabeled_PCA_Code'));
addpath(genpath('C:\Users\Pastoral\Desktop\Pursuit\Functions'));

load('YaleBCrop025');
num_groups = 3;

h = 192; w = 168; 
hh = 48; ww = 42;
sp = 1;

N = 64*num_groups;

outliers_num = 20; % number of outliers out of 64 points
rrank = 20;
v_patches = 48;
h_patches = 42;
% perm_pattern = 0 for fully shuffling; 1 for partially shuffling.
perm_flag = 1;
detect_flag = 2;

X_raw = zeros(hh*ww,N);
X_tilde = zeros(hh*ww,N);
outliers_ID = [];


for i = 1:num_groups
    X_raw(:,64*(i-1)+1:64*i) = Y(:,:,sp-1+i);
    
    shuffled_ids = randperm(64);
    outliers_id = shuffled_ids(1:outliers_num);
    if length(setdiff(outliers_id, 46))== outliers_num
        outliers_id(1) = 46;
    end
    if length(setdiff(outliers_id, 23))== outliers_num
        outliers_id(2) = 23;
    end
    if length(setdiff(outliers_id, 54))== outliers_num
        outliers_id(3) = 54;
    end
    if length(setdiff(outliers_id, 12))== outliers_num
        outliers_id(4) = 12;
    end
    [X] = permute_face(X_raw(:,(i-1)*64+1:i*64), outliers_id, hh, ww, v_patches, h_patches, perm_flag);
    X_tilde(:,(i-1)*64+1:i*64) = X;
    outliers_ID = [outliers_ID, sort(outliers_id)+64*(i-1)];

end

inliers_ID = setdiff(1:N,outliers_ID);
X_std = X_tilde./vecnorm(X_tilde,2);


%SSC_outliers_detection
affine = false ;
[Inliers_id, Outliers_id] = Outlier_Detect(X_tilde, affine, 2, 1e10);

err_in = length(setdiff(Inliers_id,inliers_ID))/length(Inliers_id);
err_out = length(setdiff(Outliers_id,outliers_ID))/length(Outliers_id);
str = ['err_in = ' num2str(err_in) ', err_out = ' num2str(err_out)];
disp(str);

%Parameters for SSC
Label_gt = zeros(N,1);
for i = 1:num_groups
    Label_gt((i-1)*64+1:i*64) = i;
end
t = Label_gt(Inliers_id);
Y = X_tilde(:,Inliers_id);

r = 0; affine = false; outlier = true; rho = 1;alpha = 20;
[missrate_in,~,~,ggrps] = SSC(Y,r,affine,alpha,outlier,rho,t);
str = ['missrate_in = ' num2str(missrate_in)];
disp(str);

label_hat = [];
for k = 1:num_groups
    idx = find(ggrps == k);
    idxx{k} = Inliers_id(idx)';
    X_sel{1,k} = X_tilde(:,idxx{k}); 
    X_sel_std{1,k} = X_std(:,idxx{k});
    % The following procedure is designed for the final assessment.
    label_hat = [label_hat; idxx{k}, k*ones(length(idxx{k}),1)];
end
Q = sortrows(label_hat);
Label_sel = Q(:,2);

Basis_estimated = zeros(hh*ww, rrank, num_groups);
Basis_estimated_std = zeros(hh*ww, rrank, num_groups);
for k = 1:num_groups
        [U,~,~] = svd(X_sel{1,k});
        [UU,~,~] = svd(X_sel_std{1,k});
       Basis_estimated(:,:,k) = U(:,1:rrank) ;
       Basis_estimated_std(:,:,k) = UU(:,1:rrank) ;
end

if (perm_flag == 0)
    [X_solved,~, Label_unsel] = Solve_all(X_tilde, Outliers_id, Basis_estimated);
else
    [X_solved,~, Label_unsel] = Solve_partial(X_tilde, Outliers_id, Basis_estimated);
end

IDX = zeros(64*num_groups,1);
for i = 1:size(X_tilde,2)
    if (ismember(i,Inliers_id) == 1)
        IDX(i,1) = Label_sel(Inliers_id == i);
    elseif (ismember(i,Outliers_id) ==1)
        IDX(i,1) = Label_unsel(Outliers_id == i);
    end
end

MISSRATE = sum(IDX(:) ~= Label_gt(:))/N;
missrate_out = sum(IDX(Outliers_id) ~= Label_gt(Outliers_id))/length(Label_gt(Outliers_id));
[err_ratio] = Evaluate(X_solved, X_tilde);
[err_ratio_in] = Evaluate(X_solved(:,inliers_ID), X_tilde(:,inliers_ID));
[err_ratio_out] = Evaluate(X_solved(:,outliers_ID),X_tilde(:,outliers_ID));
str = ['missrate_out = ' num2str(missrate_out)];
disp(str);


