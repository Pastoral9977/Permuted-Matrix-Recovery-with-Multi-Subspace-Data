clear;
clc;
format('compact');
close all
addpath(genpath('..\Functions'));
%% Data Loading
load('scores.mat');
data_gt = scoreM;
perm_flag = 0;

%% Data Inspection
Vmean = mean(data_gt);
data_central = data_gt - Vmean;
Vnorm = vecnorm(data_central);
X = data_central./Vnorm;

n = size(data_gt,1);
m = size(data_gt,2);
id_out = [2:2:m];        
id_in = setdiff(1:m, id_out);

r = 4;
all_id = 1:m;
shuffled_ratio_list = 0.1:0.1:1;
trials = 5;

num_groups = 3;
seed = randi(1000000);
rng(seed);
Err_tilde = ones(length(shuffled_ratio_list),trials);
Err3 = ones(length(shuffled_ratio_list), trials, 7);
for t_id = 1:trials
    %% data corruption
    X_tilde = zeros(n,m,length(shuffled_ratio_list));
    for s_id = 1:length(shuffled_ratio_list)
        shuffled_ratio = shuffled_ratio_list(s_id);
        [X_tilde(:,:,s_id)] = permute_corruption(X,id_out,shuffled_ratio,seed*t_id*s_id);
        [Err_tilde(s_id,t_id)] = ComputeErr(X, X_tilde(:,:,s_id), Vnorm, Vmean);
    end

    for s_id = 1:length(shuffled_ratio_list)
        shuffled_ratio = shuffled_ratio_list(s_id);
        disp(repmat('==', 1, 26))
        fprintf('Trial = %d, Shuffled Ratio: %.2f\n', t_id, shuffled_ratio);
        disp(repmat('==', 1, 26))
        X_tilde_mat = X_tilde(:,:,s_id);

      %% UPCA
        fprintf('---UPCA---\n')
        c = n-r;
        budget = 1000;
        epsilon_J = 1e-9;
        maxIter = 10000;
        delta = 1e-9;
        [B_IRLS] = DPCP_IRLS(X_tilde_mat,c,delta,maxIter,epsilon_J,budget);
        X_hat = US_mat_v1(X_tilde_mat, B_IRLS, all_id, perm_flag);
        err_upca = ComputeErr(X,  X_hat, Vnorm, Vmean);
        fprintf('err_upca = %.4f\n',err_upca)
        [Err3(s_id,t_id,1)] = err_upca;

       %% PMSDR 
        [Inliers_id, Outliers_id] = outlier_detection_with_outlier_num_known(X_tilde_mat, 7);
        [U_esti_1] = sub_esti(X_tilde_mat(:,Inliers_id),1,r);
        [U_esti_2] = sub_esti(X_tilde_mat(:,Inliers_id),2,r);
        [U_esti_3] = sub_esti(X_tilde_mat(:,Inliers_id),3,r);
         X_hat_1 = Solve_partial_LSR(X_tilde_mat,Outliers_id, U_esti_1, false);
         X_hat_2 = Solve_partial_LSR(X_tilde_mat,Outliers_id, U_esti_2, false);
         X_hat_3 = Solve_partial_LSR(X_tilde_mat,Outliers_id, U_esti_3, false);
         err_pmsdr_1 = ComputeErr(X,  X_hat_1, Vnorm, Vmean);
         err_pmsdr_2 = ComputeErr(X,  X_hat_2, Vnorm, Vmean);
         err_pmsdr_3 = ComputeErr(X,  X_hat_3, Vnorm, Vmean);
         [Err3(s_id,t_id,2)] = err_pmsdr_1;
         [Err3(s_id,t_id,3)] = err_pmsdr_2;
         [Err3(s_id,t_id,4)] = err_pmsdr_3;
         fprintf('err_pmsdr_1 = %.4f\n',err_pmsdr_1)
         fprintf('err_pmsdr_2 = %.4f\n',err_pmsdr_2)
         fprintf('err_pmsdr_3 = %.4f\n',err_pmsdr_3)

       %% RPCA
         fprintf('---RPCA---')
         lambdas_rpca = [0.2/sqrt(size(X_tilde_mat,1)), 0.2/sqrt(size(X_tilde_mat,2)),  ...
                        1/sqrt(size(X_tilde_mat,1)),   1/sqrt(size(X_tilde_mat,2)),    ...
                        5/sqrt(size(X_tilde_mat,1)),   5/sqrt(size(X_tilde_mat,2)),   ];
         err_rpca = inf; a = nan;
         for kk = 1:length(lambdas_rpca)
             lambda = lambdas_rpca(kk);
             mu = 10*lambda;
             X_rpca = RobustPCA(X_tilde_mat, lambda, mu, 1e-6);
             err = ComputeErr(X,  X_rpca, Vnorm, Vmean);
             if err < err_rpca
                 err_rpca = err;
                 a = lambda;
             end
         end
         fprintf('lambda_rpca = %.4f\n', a)
         fprintf('err_rpca = %.4f\n',err_rpca)
         [Err3(s_id,t_id,5)] = err_rpca;

       %% RKPCA
        fprintf('---RKPCA---')
        ker.type='rbf';
        ker.par=0;
        options.p=1;
        err_rkpca = inf; a = nan;
        if shuffled_ratio < 0.5
            lambdas_rkpca = 0.5:0.05:0.95;
        else
            lambdas_rkpca = 0.1:0.05:0.5;
        end
        for kk = 1:length(lambdas_rkpca)
            lambda = lambdas_rkpca(kk);
            X_rkpca = RKPCA_PLMAdSS(X_tilde_mat,ker,lambda,options);
            err = ComputeErr(X, X_rkpca, Vnorm, Vmean);
            if err < err_rkpca
                err_rkpca = err;
                a = lambda;
            end
        end
        fprintf('lambda_rkpca = %.2f\n', a)
        fprintf('err_rkpca = %.4f\n',err_rkpca)
        [Err3(s_id,t_id,6)] = err_rkpca;

       %% SSC
        fprintf('---SSC---')
        affines_ssc = [false, true];
        err_ssc = inf; a = nan; b = nan;
        if shuffled_ratio < 0.5
            alphas_ssc = [100,200,500];
        else
            alphas_ssc = [5,20,100];
        end
        for kk = 1:length(alphas_ssc)
            alpha = alphas_ssc(kk);
            for jj = 1:length(affines_ssc)
                affine = affines_ssc(jj);
                CMat = admmOutlier_mat_func(X_tilde_mat,affine,alpha);
                N = size(X_tilde_mat,2);
                C = CMat(1:N,:);
                X_ssc = X_tilde_mat*C;
                err = ComputeErr(X, X_ssc, Vnorm, Vmean);
                if err < err_ssc
                    err_ssc = err;
                    a = alpha;
                    b = affine;
                end
            end
        end
        fprintf('alpha = %d, affine = %d\n', a, b)
        fprintf('err_ssc = %.4f\n',err_ssc)
        [Err3(s_id,t_id,7)] = err_ssc;
                 
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

Err_rpca = Err3(:,:,5)';
mean_rpca = mean(Err_rpca);
std_rpca = std(Err_rpca);

Err_rkpca = Err3(:,:,6)';
mean_rkpca = mean(Err_rkpca);
std_rkpca = std(Err_rkpca);

Err_ssc = Err3(:,:,7)';
mean_ssc = mean(Err_ssc);
std_ssc = std(Err_ssc);

mean_tilde = mean(Err_tilde');
std_tilde = std(Err_tilde');

run plotErr.m
