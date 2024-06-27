close all
clear 

addpath(genpath('..\Copy_SSC_ADMM_v1.1'));
addpath(genpath('..\Functions'));

% load('1R2RC_truth.mat');
load('1R2RC_g13_truth.mat');
% load('1RT2RCR_g12_truth.mat');
% load('articulated_g12_truth.mat');
% load('dancing_truth.mat')
% load('kanatani3_truth.mat')
n = max(s); % s来源于load来的数据，为每个点所属的cluster类别(allpoints_class_gt)，则n为num_groups
N = size(x,2);
F = size(x,3); % x来源于load来的数据，为3*num_points*num_frames
D = 2*F;
% 下述代码reshape复合permute，只取每个三维点的前两个维度，并将每个点的所有时刻（frame）的这两个维度拼接成大向量，
% 即维度为上述D = 2*F
X_gt = reshape(permute(x(1:2,:,:),[1 3 2]),D,N); 

rng(42);
%generate the corrupted data
shuffled_ratio = 0.4 ;outlier_ratio = 0.4;
num_shuffled = fix(N*outlier_ratio);
outliers_ID = randperm(N,num_shuffled);
inliers_ID = setdiff(1:N,outliers_ID);

X_tilde = X_gt;
for i = 1:num_shuffled
    j = outliers_ID(i);
    X_tilde(:,j) = shuffle(X_gt(:,j),shuffled_ratio);
end

% detect_flag = 1 ;affine2 = true; alpha2 = 1e10;
% [Inliers_id, Outliers_id] = Outlier_Detect_v2(cnormalize(X_tilde), detect_flag, alpha2, -1, affine2);
alpha = 5;
lambda = 0.95;
thres = 1e-3;
[Inliers_id, Outliers_id] = Outlier_Detect_v3(X_tilde, alpha, lambda, thres);
remark(Inliers_id, Outliers_id, inliers_ID, outliers_ID);

%Parameters for SSC
r = 0; affine = true; outlier = true; rho = 0.7; alpha = 20;
[missrate_in,~,~,ggrps] = SSC(X_tilde(:,Inliers_id),r,affine,alpha,outlier,rho,s(Inliers_id));
fprintf('\tmissrate_in = %.4f\n\n', missrate_in);

label_hat = zeros(length(Inliers_id), 2);
X_selected_inlier_groups = cell(1, n);
index = 1;
for k = 1:n
    kth_group_ids = Inliers_id(ggrps == k);
    X_selected_inlier_groups{k} = X_tilde(:,kth_group_ids); 
    % The following procedure is designed for the final assessment.
    num_kth_group_ids = length(kth_group_ids);
    label_hat(index:index+num_kth_group_ids-1, :) = [reshape(kth_group_ids, num_kth_group_ids, 1), k*ones(num_kth_group_ids, 1)];
    index = index + num_kth_group_ids;
end
Q = sortrows(label_hat);
inlier_class_solved = Q(:,2);

U_esti = zeros(D, 4, n);
for k = 1:n
   [U,~,~] = svd(X_selected_inlier_groups{1,k});
   U_esti(:,:,k) = U(:,1:4);
end





















