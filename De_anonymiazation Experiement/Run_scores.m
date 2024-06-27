clear;format('compact');
close all
addpath(genpath('..\Copy_SSC_ADMM_v1.1'));
addpath(genpath('..\Functions'));
%% Data Loading
load('scores.mat');
data_raw = scoreM;
perm_flag = 0;

%% Data Inspection
Vmean = mean(data_raw);
data_central = data_raw - Vmean;
Vnorm = vecnorm(data_central);
X = data_central./Vnorm;

m = size(data_raw,1);
n = size(data_raw,2);
id_out = [8:14];        
id_in = setdiff(1:n, id_out);

r = 3;
O_id = 1:n;
shuffled_ratio_list = 0.1:0.1:1;
trials = 3;
LRWC_name_list = {'LSR','MUL'};

num_groups = 3;
for t_id = 1:trials
    %% data corruption
    X_tilde = zeros(m,n,length(shuffled_ratio_list));
    for shuffled_ratio = shuffled_ratio_list
        s_id = fix(shuffled_ratio*10);
        [X_tilde(:,:,s_id)] = permute_corruption(X, id_out, shuffled_ratio);
        [Err_tilde(s_id,t_id)] = ComputeErr(X, X_tilde(:,:,s_id), Vnorm, Vmean);
    end
    %% UPCA
    for l_id = 1:length(LRWC_name_list)
        LRWC_name = LRWC_name_list{l_id};
            for shuffled_ratio = shuffled_ratio_list
                s_id = fix(shuffled_ratio*10);
                X_tilde_mat = X_tilde(:,:,s_id);
                
                if (l_id == 1)
                    c = m-r;
                    budget = 1000;
                    epsilon_J = 1e-9;
                    maxIter = 3000;
                    delta = 1e-9;
                    [B_IRLS] = DPCP_IRLS(X_tilde_mat, c, delta, maxIter,epsilon_J,budget);
                    X_hat = US_mat_v1(X_tilde_mat, B_IRLS, O_id, perm_flag);  
                    [Err3(s_id,t_id,l_id)] = ComputeErr(X,  X_hat, Vnorm, Vmean);
                end
                
                if(l_id == 2)
                    affine2 = false ;detect_flag = 0;validation = -1;
                    alpha2 = 500;
%                     alpha2 = 1e9;
                    [Inliers_id, Outliers_id] = Outlier_Detect_v2(X_tilde_mat,detect_flag,alpha2,validation,affine2);
                    remark(Inliers_id, Outliers_id, id_in, id_out);
                    [Basis_estimated_1] = sub_esti(X_tilde_mat(:,Inliers_id),1,r);
                    [Basis_estimated_2] = sub_esti(X_tilde_mat(:,Inliers_id),2,r);
                    [Basis_estimated_3] = sub_esti(X_tilde_mat(:,Inliers_id),3,r);
                     X_hat_1 = Solve_partial_LSR(X_tilde_mat,Outliers_id, Basis_estimated_1);
                     X_hat_2 = Solve_partial_LSR(X_tilde_mat,Outliers_id, Basis_estimated_2);
                     X_hat_3 = Solve_partial_LSR(X_tilde_mat,Outliers_id, Basis_estimated_3);
                     [Err3(s_id,t_id,l_id)] = ComputeErr(X,  X_hat_1, Vnorm, Vmean);
                     [Err3(s_id,t_id,3)] = ComputeErr(X,  X_hat_2, Vnorm, Vmean);
                     [Err3(s_id,t_id,4)] = ComputeErr(X,  X_hat_3, Vnorm, Vmean);
                end
                
            end
    end
end

Err_LSR = Err3(:,:,1)';
mean_LSR = mean(Err_LSR);
std_LSR = std(Err_LSR);
Err_MUL_1 = Err3(:,:,2)';
mean_MUL_1 = mean(Err_MUL_1);
std_MUL_1 = std(Err_MUL_1);
Err_MUL_2 = Err3(:,:,3)';
mean_MUL_2 = mean(Err_MUL_2);
std_MUL_2 = std(Err_MUL_2);
Err_MUL_3 = Err3(:,:,4)';
mean_MUL_3 = mean(Err_MUL_3);
std_MUL_3 = std(Err_MUL_3);
mean_tilde = mean(Err_tilde');
std_tilde = std(Err_tilde');

run plotErr.m
