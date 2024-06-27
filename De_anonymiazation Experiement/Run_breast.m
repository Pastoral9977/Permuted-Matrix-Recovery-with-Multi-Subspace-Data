clear;format('compact');
close all
addpath(genpath('..\Copy_SSC_ADMM_v1.1'));
addpath(genpath('..\Functions'));
% addpath(genpath('C:\Users\Pastoral\Desktop\Copy_SSC_ADMM_v1.1'));
% addpath(genpath('C:\Users\Pastoral\Desktop\Unlabeled_PCA_Code'));
% addpath(genpath('C:\Users\Pastoral\Desktop\Pursuit\Functions'));
%% Data Loading
load('breast_benign.mat');
data_gt = benign;
perm_flag = 0;

%% Data Inspection
Vmean = mean(data_gt);  % Matlab中mean函数对矩阵作用是对每列求均值，即(m,n)->(1,n);
data_central = data_gt - Vmean;
Vnorm = vecnorm(data_central); % Matlab中vecnorm函数对矩阵作用同样是作用于列，即(m,n)->(1,n);
X = data_central./Vnorm;

m = size(data_gt,1);
n = size(data_gt,2);
id_out = 2:2:n;
id_in = setdiff(1:n, id_out);

r = 4;
all_ids = 1:n;
outlier_ratio = 0.5;
shuffled_ratio_list = 0.1:0.1:1;
trials = 2;
LRWC_name_list = {'LSR','MUL'};

ll = length(LRWC_name_list);
num_groups = 3;

Err3 = ones(length(shuffled_ratio_list), trials, ll);
Err_tilde = zeros(length(shuffled_ratio_list), trials);
for t_id = 1:trials
    %% data corruption
    X_tilde = zeros(m,n,length(shuffled_ratio_list));
    for shuffled_ratio = shuffled_ratio_list
        s_id = fix(shuffled_ratio*10);
        [X_tilde(:,:,s_id)] = permute_corruption(X, id_out, shuffled_ratio);
        [Err_tilde(s_id,t_id)] = ComputeErr(X, X_tilde(:,:,s_id), Vnorm, Vmean);
    end
    %% UPCA
    for l_id = 1:ll
        LRWC_name = LRWC_name_list{l_id};
        for shuffled_ratio = shuffled_ratio_list
            s_id = fix(shuffled_ratio*10);
            X_tilde_mat = X_tilde(:,:,s_id);

            if (l_id == 1)
                c = m-r;
                budget = 1000;
                epsilon_J = 1e-9;
                maxIter = 10000;
                delta = 1e-9;
                [B_IRLS] = DPCP_IRLS(X_tilde_mat, c, delta, maxIter,epsilon_J,budget);
                X_hat = US_mat_v1(X_tilde_mat, B_IRLS, all_ids, perm_flag);  
                [Err3(s_id,t_id,l_id)] = ComputeErr(X,  X_hat, Vnorm, Vmean);
            end

            if(l_id == 2)
                [Inliers_id, Outliers_id] = Outlier_Detect_v3(X_tilde_mat);
                remark(Inliers_id, Outliers_id, id_in, id_out);
                [U_esti_1] = sub_esti(X_tilde_mat(:,Inliers_id),1,r);
                [U_esti_2] = sub_esti(X_tilde_mat(:,Inliers_id),2,r);
                [U_esti_3] = sub_esti(X_tilde_mat(:,Inliers_id),3,r);
                X_hat_1 = Solve_partial_LSR(X_tilde_mat,all_ids, U_esti_1);
                X_hat_2 = Solve_partial_LSR(X_tilde_mat,all_ids, U_esti_2);
                X_hat_3 = Solve_partial_LSR(X_tilde_mat,all_ids, U_esti_3);
                [Err3(s_id,t_id,l_id)] = ComputeErr(X,  X_hat_1, Vnorm, Vmean);
                [Err3(s_id,t_id,3)] = ComputeErr(X,  X_hat_2, Vnorm, Vmean);
                [Err3(s_id,t_id,4)] = ComputeErr(X,  X_hat_3, Vnorm, Vmean);
            end
        end
    end
end

Err_LSR     = Err3(:,:,1)';
mean_LSR    = mean(Err_LSR);
std_LSR     = std(Err_LSR);
Err_MUL_1   = Err3(:,:, 2)';
mean_MUL_1  = mean(Err_MUL_1);
std_MUL_1   = std(Err_MUL_1);
Err_MUL_2   = Err3(:,:, 3)';
mean_MUL_2  = mean(Err_MUL_2);
std_MUL_2   = std(Err_MUL_2);
Err_MUL_3   = Err3(:,:,4)';
mean_MUL_3  = mean(Err_MUL_3);
std_MUL_3   = std(Err_MUL_3);
mean_tilde  = mean(Err_tilde, 2);
std_tilde   = std(Err_tilde, 0, 2);

run plotErr.m
