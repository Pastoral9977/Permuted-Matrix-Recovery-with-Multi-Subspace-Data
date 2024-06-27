close all
clear 
addpath(genpath('..\Copy_SSC_ADMM_v1.1'));
addpath(genpath('..\Functions'));


% load('1R2RC_truth.mat');
load('1R2RC_g12_truth.mat');
% load('1RT2RCR_g12_truth.mat');
% load('articulated_g12_truth.mat');
% load('dancing_truth.mat')
%  load('two_cranes_g23_truth.mat')
% load('two_cranes_g13_truth.mat')
% load('1R2RC_g13_truth.mat');
n = max(s);
N = size(x,2);
F = size(x,3);
D = 2*F;
X_raw = reshape(permute(x(1:2,:,:),[1 3 2]),D,N);

%generate the corrupted data
shuffled_ratio = 0.4 ;outlier_ratio = 0.5;
num_shuffled = fix(N*outlier_ratio);
outliers_ID = randperm(N,num_shuffled);
inliers_ID = setdiff(1:N,outliers_ID);

X_tilde = X_raw;
for i = 1:num_shuffled
    j = outliers_ID(i);
    X_tilde(:,j) = shuffle(X_raw(:,j),shuffled_ratio);
end


X_std = X_tilde./vecnorm(X_tilde,2);

affine = true;alpha = 10^9;

detect_flag = 2;
[Inliers_id, Outliers_id, Z] = Outlier_Detect_v2(X_std, affine, detect_flag, alpha);
value = Z;
T = kmeans(value,2);
I = find(T == 1); J = find(T == 2);
I_mean = mean(value(I));J_mean = mean(value(J));

err_in = length(setdiff(Inliers_id,inliers_ID))/length(Inliers_id);
err_out = length(setdiff(Outliers_id,outliers_ID))/length(Outliers_id);
remark(Inliers_id, Outliers_id, inliers_ID, outliers_ID);

Z_s = Z([inliers_ID,outliers_ID]);
value = vecnorm(Z_s,1);
figure
plot(value);grid on; hold on;
line([length(inliers_ID) length(inliers_ID)],[0,max(I_mean,J_mean)*1.3],'linestyle','--', 'Color','r', 'LineWidth', 1.5);
title('1-norm of the consequential coefficients for clustering')
text(0.4*length(inliers_ID),max(I_mean,J_mean)*0.5,'inliers','color','red','FontSize',16);
text(N-0.6*length(outliers_ID),max(I_mean,J_mean)*0.5,'outliers','color','red','FontSize',16);






