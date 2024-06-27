close all
clear 
addpath(genpath('..\Copy_SSC_ADMM_v1.1'));
addpath(genpath('..\Functions'));

rng(21)
store = true;
D = 50;
V = 120;
num_groups_set = 1; % 不同的num_groups取值
outlier_ratio = 0.6;

rn_ratio_set = 0.02:0.02:0.5;
shuffled_ratio_set = 0.1:0.1:0.6;

num_sr = length(shuffled_ratio_set);
num_rank = length(rn_ratio_set);
SNR = 40;

RESULT = struct();

RESULT.recover_error_singlegroup = cell(1, length(num_groups_set));


diary('DataRecovery_single_subspace.txt');
for g_idx = 1: length(num_groups_set)
    num_groups = num_groups_set(g_idx);
    recover_error_singlesubspace_mat = ones(num_rank, num_sr);
    
    for r_idx = 1:num_rank
        rn_ratio = rn_ratio_set(r_idx);
        rrank = min(int64(D * rn_ratio), fix(V*(1-outlier_ratio)*0.9));
        [X_gt, U_gt, allpoints_class_gt] = Generate_data(D, V, num_groups, rrank);

        for sr_idx = 1:num_sr
            shuffled_ratio = shuffled_ratio_set(sr_idx);
            
            disp(repmat('==', 1, 46))
            fprintf('Trial: %d/%d, D: %d, V: %d, Groups: %d, Rank: %d, Shuffled Ratio: %.2f, Outlier Ratio: %.2f\n', ...
                            (g_idx-1)*num_rank*num_sr+(r_idx-1)*num_sr+sr_idx, length(num_groups_set)*num_rank*num_sr, D, V, num_groups, rrank, shuffled_ratio, outlier_ratio);
            disp(repmat('==', 1, 46))
            
            M_gt = reshape(X_gt, [D, V, num_groups]);
            [X_tilde, outliers_ID, inliers_ID] = Generate_observed_data_noisy(M_gt, outlier_ratio, shuffled_ratio, SNR);
            
            %% Recovery and Assessment
            fprintf('\tPart1: Recovery and Assessment\n');
            X_solved = Solve_single_subspace(X_tilde, outliers_ID, U_gt);
            [recover_err_gtbases_singlesubspace, ~] = EvaluateRefined_v3(X_solved, X_gt, outliers_ID, outliers_ID, num_groups, rrank); 
            recover_error_singlesubspace_mat(r_idx, sr_idx) = recover_err_gtbases_singlesubspace;
            fprintf('\t\tSingle subspace recover error with ground truth bases = %0.4f\n', recover_err_gtbases_singlesubspace);
            
        end
    end
    
    %% Store the RESULT
    if store
        valid_shuffled_ratio_set = arrayfun(@(x) matlab.lang.makeValidName(sprintf('SRatio_%.1f', x)), shuffled_ratio_set, 'UniformOutput', false);
        valid_rn_ratio_set = arrayfun(@(x) matlab.lang.makeValidName(sprintf('RRatio_%.2f', x)), rn_ratio_set, 'UniformOutput', false);

        RESULT.recover_error_singlegroup{g_idx} = array2table(recover_error_singlesubspace_mat,...
            'VariableNames', valid_shuffled_ratio_set, ...
            'RowNames', valid_rn_ratio_set);
    end
end
diary off;
save('RESULT_single_subspace.mat', 'RESULT');

function X_solved = Solve_single_subspace(X_tilde, outliers_ID, U_gt)
    n = length(outliers_ID);
    X_solved = X_tilde;
    for i = 1:n
        id = outliers_ID(i);
        y = X_tilde(:, id);
        c_hat = LSR(U_gt, y);
        X_solved(:,id) = U_gt * c_hat;
    end
end
