close all
clear 
addpath(genpath('C:\Users\Pastoral\Desktop\Copy_SSC_ADMM_v1.1'));
addpath(genpath('C:\Users\Pastoral\Desktop\Unlabeled_PCA_Code'));
addpath(genpath('C:\Users\Pastoral\Desktop\Pursuit\Functions'));

load('HARtrain.mat');

[n,N] = size(Xtrain);
rank = 12;num_groups = max(ytrain);
%generate the corrupted data
shuffled_ratio = 0.4;outlier_ratio = 0.5;
num_shuffled = fix(N*outlier_ratio);
outliers_ID = sort(randperm(N,num_shuffled));
inliers_ID = setdiff(1:N,outliers_ID);
X_tilde = Xtrain;

for i = 1:num_shuffled
    j = outliers_ID(i);
    X_tilde(:,j) = shuffle(Xtrain(:,j),shuffled_ratio);
end

X_std = X_tilde;

%Parameters for SSC
affine = false; rho = 1;alpha = 10^9;
%Parameters for DPCP
c = n-rank;budget = 10000;epsilon_J = 1e-9;maxIter = 1000;delta = 1e-9;

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

%SSC
r = 0; alpha2 = 2000; affine2 = true; rho2 = 0.7;
outlier = false;
Y = X_tilde(:,Inliers_id);
t = ytrain(Inliers_id)';
[missrate_in,C,grps,ggrps] = SSC(Y,r,affine2,alpha2,outlier,rho2,t);


% Pre-setting
X_sel = cell(1,num_groups);idxx = cell(1,num_groups); 
Basis_estimated = zeros(n, rank, num_groups);label_hat = [];
% record the estimated labels of the selected data and
% group the selected data.
for k = 1:num_groups
    idx = find(ggrps == k);
    idxx{k} = Inliers_id(idx);
    X_sel{1,k} = X_tilde(:,idxx{k});
    label_hat = [label_hat; idxx{k}, k*ones(length(idxx{k}),1)];
end
Q = sortrows(label_hat);
Label_sel = Q(:,2);

for k = 1:num_groups
    [U, ~, ~] = svd(X_sel{1,k});
    Basis_estimated(:,:,k) = U(:,1:rank);
end

if (shuffled_ratio >= 0.5)
        [X_solved,~, Label_unsel] = Solve_all(X_tilde, Outliers_id, Basis_estimated);
else
        [X_solved,~, Label_unsel] = Solve_partial(X_tilde, Outliers_id, Basis_estimated);
end

IDX_tol = zeros(size(X_tilde,2),1);

for i = 1:size(X_tilde,2)
    if (ismember(i,Inliers_id) == 1)
        IDX_tol(i,1) = Label_sel(Inliers_id == i);
    elseif (ismember(i,Outliers_id) ==1)
        IDX_tol(i,1) = Label_unsel(Outliers_id == i);
    end
end



MISSRATE = sum(IDX_tol(:) ~= ytrain(:))/length(ytrain);
err_ratio_out= Evaluate(X_solved(:,outliers_ID), Xtrain(:,outliers_ID));








