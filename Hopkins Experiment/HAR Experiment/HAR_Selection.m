close all
clear 
addpath(genpath('C:\Users\Pastoral\Desktop\Copy_SSC_ADMM_v1.1'));
addpath(genpath('C:\Users\Pastoral\Desktop\Unlabeled_PCA_Code'));
addpath(genpath('C:\Users\Pastoral\Desktop\Pursuit\Functions'));

load('HARtrain.mat');

[n,N] = size(Xtrain);
%generate the corrupted data
shuffled_ratio = 0.45 ;outlier_ratio = 0.5;
num_shuffled = fix(N*outlier_ratio);
outliers_ID = randperm(N,num_shuffled);
inliers_ID = setdiff(1:N,outliers_ID);
X_tilde = Xtrain ;
detect_flag = 0;

for i = 1:num_shuffled
    j = outliers_ID(i);
    X_tilde(:,j) = shuffle(Xtrain(:,j),shuffled_ratio);
end

X_std = X_tilde;

%Parameters for SSC
r = 0; affine = false; rho = 1;alpha = 10^9;
Z = admmLasso_mat_func(X_std, affine, alpha); 


value2 = vecnorm(Z,1);
T = kmeans(value2',2);
I = find(T == 1); J = find(T == 2);
I_mean = mean(value2(:,I));
J_mean = mean(value2(:,J));
if (I_mean > J_mean)
    Inliers_id = J;Outliers_id = I;
else
    Inliers_id = I;Outliers_id = J;
end

err_in = length(setdiff(Inliers_id,inliers_ID))/length(Inliers_id);
err_out = length(setdiff(Outliers_id,outliers_ID))/length(Outliers_id);

Z_s = Z(:,[inliers_ID,outliers_ID]);
value = vecnorm(Z_s,1);
figure
plot(value);grid on; hold on;
line([length(inliers_ID) length(inliers_ID)],[0,max(I_mean,J_mean)+3],'linestyle','--', 'Color','r', 'LineWidth', 1.5);
title('1-norm of the consequential coefficients for clustering')
text(0.4*length(inliers_ID),1.5,'inliers','color','red','FontSize',16);
text(N-0.6*length(outliers_ID),1.5,'outliers','color','red','FontSize',16);




