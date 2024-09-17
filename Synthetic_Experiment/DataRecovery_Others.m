close all
clear 
addpath(genpath('..\Functions'));
addpath(genpath(pwd));


store = true;
D = 50;
V = 120;
num_groups_set = [2,3,5,8,10];
outlier_ratio = 0.4;

rn_ratio_set = 0.02:0.02:0.5;
shuffled_ratio_set = 0.1:0.1:0.6;
c = 0.1;
verbose = false;
tol = 1e-5;
trials_for_sub = 3;

num_sr = length(shuffled_ratio_set);
num_rank = length(rn_ratio_set);

for g_idx = 1: length(num_groups_set)
    num_groups = num_groups_set(g_idx);
    recovery_error_RKPCA_mat = ones(num_rank, num_sr);
    recovery_error_PMSDR_RKPCA_mat = ones(num_rank, num_sr);
    recovery_error_RPCA_mat = ones(num_rank, num_sr);
    recovery_error_PMSDR_RPCA_mat = ones(num_rank, num_sr);
    recovery_error_SSC_mat = ones(num_rank, num_sr);
    recovery_error_PMSDR_SSC_mat = ones(num_rank, num_sr);
    
    for r_idx = 1:num_rank
        seed = randi(100000);
        rng(seed)
        rn_ratio = rn_ratio_set(r_idx);
        rrank = min(int64(D * rn_ratio), fix(V*(1-outlier_ratio)*0.9));
        [X_gt, U_gt, allpoints_class_gt] = Generate_data(D, V, num_groups, rrank);

        for sr_idx = 1:num_sr
            rng(seed)
            shuffled_ratio = shuffled_ratio_set(sr_idx);
            disp(repmat('==', 1, 55))
            fprintf('Trial: %d/%d, D: %d, V: %d, Groups: %d, Rank: %d, Shuffled Ratio: %.2f, Outlier Ratio: %.2f, seed: %d\n', ...
                            (g_idx-1)*num_rank*num_sr+(r_idx-1)*num_sr+sr_idx, length(num_groups_set)*num_rank*num_sr, D, V, num_groups, rrank, shuffled_ratio, outlier_ratio, seed);
            disp(repmat('==', 1, 55))
            X_tilde = permute_data(X_gt, shuffled_ratio, outlier_ratio, allpoints_class_gt);
            
           %% 
            fprintf('Part1: Complete Recovery and Assessment\n');
            
           %% rkpca
            rng(seed)
            lambdas_rkpca = [0.05:0.05:0.5];
            fprintf('\tStep1: checking best parameters for whole rkpca\n');
            ker.type='rbf';
            ker.par=0;
            options.p=1;
            options.c=c;
            options.verbose=verbose;
            err_rkpca = inf;
            display_progress()
            for ii = 1:length(lambdas_rkpca)
                par = lambdas_rkpca(ii);
                X_RKPCA = RKPCA_PLMAdSS(X_tilde,ker,par,options);
                Err = Evaluate(X_RKPCA, X_gt);
                if Err<err_rkpca
                    err_rkpca = Err;
                    recovery_error_RKPCA_mat(r_idx, sr_idx) = Err;
                    lambda_rkpca = par;
                end
                update_progress(ii, length(lambdas_rkpca))
            end
            fprintf('\t\trecovery_err_rkpca = %.4f\n', err_rkpca)
            fprintf('\t\tlambda_rkpca = %.2f\n', lambda_rkpca)
            
           %% rpca
            rng(seed)
            err_rpca = inf;
            lambdas_rpca = [0.2/sqrt(D), 0.2/sqrt(V), 0.2/sqrt(num_groups*V), ...
                            1/sqrt(D),   1/sqrt(V),   1/sqrt(num_groups*V),...
                            5/sqrt(D),   5/sqrt(V),   5/sqrt(num_groups*V)];
            fprintf('\tStep2: checking best parameters for whole rpca\n');
            display_progress()
            
%             for pp = 1:
            for ii = 1:length(lambdas_rpca)
                par = lambdas_rpca(ii);
                mu = 10*par;
                [X_RPCA, S] = RobustPCA(X_tilde, par, mu, tol);
                Err = Evaluate(X_RPCA, X_gt);
                if Err<err_rpca
                    err_rpca = Err;
                    recovery_error_RPCA_mat(r_idx, sr_idx) = Err;
                    lambda_rpca = par;
                end
                update_progress(ii, length(lambdas_rpca))
            end
            fprintf('\t\tlambda_rpca = %.2f\n', lambda_rpca)
            fprintf('\t\trecovery_err_rpca = %.4f\n', err_rpca)
            
           %% ssc
            rng(seed)
            alphas = [5,20,100,500];
            affines = [false, true];
            outliers = [false, true];
            err_ssc = inf;
            fprintf('\tStep3: checking best parameters for whole ssc\n');
            display_progress()
            for ii = 1:length(alphas)
                par1 = alphas(ii);
                for jj = 1:length(affines)
                    par2 = affines(jj);
                    for kk = 1:length(outliers)
                        par3 = outliers(kk);
                        if par3
                            CMat = admmOutlier_mat_func(X_tilde,par2,par1);
                        else
                            CMat = admmLasso_mat_func(X_tilde,par2,par1);
                        end
                        N = size(X_tilde,2);
                        C = CMat(1:N,:);
                        X_SSC = X_tilde*C;
                        Err = Evaluate(X_SSC, X_gt);
                        if Err<err_ssc
                            err_ssc = Err;
                            recovery_error_SSC_mat(r_idx, sr_idx) = Err;
                            alpha = par1;
                            affine = par2;
                            outlier = par3;
                        end
                        update_progress(kk+(jj-1)*length(affines)+(ii-1)*length(affines)*length(outliers), length(affines)*length(alphas)*length(outliers))
                    end
                end
            end
            fprintf('\t\talpha_ssc = %.2f, affine_ssc = %d, outlier = %d\n', alpha, affine, outlier)
            fprintf('\t\trecovery_err_ssc = %.4f\n', err_ssc)     
            
           %% 
            fprintf('Part2: Sub Recovery and Assessment (with gt subspace cluster)\n');
            
           %% rkpca PMSDR
            ker.type='rbf';
            ker.par=0;
            options.p=1;
            options.c=c;
            options.verbose=verbose;
            X_PMSDR_RKPCA = zeros(size(X_tilde));
            errs_for_sub = ones(1, trials_for_sub);
            for pp = 1:trials_for_sub
                rng(seed*pp)
                X_tilde = permute_data(X_gt, shuffled_ratio, outlier_ratio, allpoints_class_gt);
                for ii = 1:num_groups
                    X_rkpca = RKPCA_PLMAdSS(X_tilde(:, (ii-1)*V+1:ii*V),ker,lambda_rkpca,options); 
                    X_PMSDR_RKPCA(:, allpoints_class_gt==ii) = X_rkpca;
                end
                Err = Evaluate(X_PMSDR_RKPCA, X_gt);
                fprintf('\t\t\trecovery_err_rkpca_PMSDR_trial = %.4f\n', Err)
                errs_for_sub(pp) = Err;
            end
            recovery_error_PMSDR_RKPCA_mat(r_idx, sr_idx) = median(errs_for_sub);
            fprintf('\t\trecovery_err_rkpca_PMSDR = %.4f\n', median(errs_for_sub))
            
           %% rpca PMSDR
            X_PMSDR_RPCA = zeros(size(X_tilde));
            errs_for_sub = ones(1, trials_for_sub);
            for pp = 1:trials_for_sub
                rng(seed*pp)
                X_tilde = permute_data(X_gt, shuffled_ratio, outlier_ratio, allpoints_class_gt);
                for ii = 1:num_groups
                    mu = 10*lambda_rpca;
                    X_rpca = RobustPCA(X_tilde(:, (ii-1)*V+1:ii*V), lambda_rpca, mu, tol);
                    X_PMSDR_RPCA(:, allpoints_class_gt==ii) = X_rpca;
                end
                Err = Evaluate(X_PMSDR_RPCA, X_gt);
                fprintf('\t\t\trecovery_err_rpca_PMSDR_trial = %.4f\n', Err)
                errs_for_sub(pp) = Err;
            end
            recovery_error_PMSDR_RPCA_mat(r_idx, sr_idx) = median(errs_for_sub);
            fprintf('\t\trecovery_err_rpca_PMSDR = %.4f\n', median(errs_for_sub))
            
           %% ssc PMSDR
            X_PMSDR_SSC = zeros(size(X_tilde));
            errs_for_sub = ones(1, trials_for_sub);
            for pp = 1:trials_for_sub
                rng(seed*pp)
                X_tilde = permute_data(X_gt, shuffled_ratio, outlier_ratio, allpoints_class_gt);
                for ii = 1:num_groups
                    XX = X_tilde(:, allpoints_class_gt==ii);
                    if outlier
                        CMat = admmOutlier_mat_func(XX,affine,alpha);
                    else
                        CMat = admmLasso_mat_func(XX,affine,alpha);
                    end
                    N = size(XX,2);
                    C = CMat(1:N,:);
                    X_PMSDR_SSC(:, allpoints_class_gt==ii) = XX*C;
                end
                Err = Evaluate(X_PMSDR_SSC, X_gt);
                fprintf('\t\t\trecovery_err_rpca_PMSDR_trial = %.4f\n', Err)
                errs_for_sub(pp) = Err;
            end
            recovery_error_PMSDR_SSC_mat(r_idx, sr_idx) = median(errs_for_sub);
            fprintf('\t\trecovery_err_ssc_PMSDR = %.4f\n', median(errs_for_sub))
        end
    end
    
end
save('RESULT_Others.mat');
run DataRecovery_Others_plot.m