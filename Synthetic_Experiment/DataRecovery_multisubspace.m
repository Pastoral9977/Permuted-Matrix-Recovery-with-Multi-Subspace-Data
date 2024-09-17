close all
clear 
addpath(genpath(pwd));
addpath(genpath('..\Functions'));

use_diary = true;
broadcast = ~use_diary;
store = true;
D = 50;
V = 120;
num_groups_set = [2,3,5,8,10]; 

outlier_ratio = 0.6;

rn_ratio_set = 0.02:0.02:0.5;
shuffled_ratio_set = 0.1:0.1:0.6;

num_sr = length(shuffled_ratio_set);
num_rank = length(rn_ratio_set);

RESULT = struct();
RESULT.missrate_out_allgroup = cell(1, length(num_groups_set));
RESULT.missrate_out_gtbases_allgroup = cell(1, length(num_groups_set));
RESULT.recover_error_allgroup = cell(1, length(num_groups_set));
RESULT.recover_error_gtbases_allgroup = cell(1, length(num_groups_set));
RESULT.missrate_in_allgroup = cell(1, length(num_groups_set));
RESULT.undetected_outlier_ratio_allgroup = cell(1, length(num_groups_set));

if exist('DataRecovery_Group.txt', 'file') == 2
    delete('DataRecovery_Group.txt');
end
if use_diary
    diary('DataRecovery_Group.txt');
end

for g_idx = 1: length(num_groups_set)
    num_groups = num_groups_set(g_idx);
    
    missrate_out_gtbases_mat = ones(num_rank, num_sr);
    missrate_out_mat = ones(num_rank, num_sr);
    recover_error_gtbases_mat = ones(num_rank, num_sr);
    recover_error_mat = ones(num_rank, num_sr);
    recover_error_detectedpart_mat = ones(num_rank, num_sr);
    missrate_in_mat = ones(num_rank, num_sr);
    undetected_outlier_ratio_mat = ones(num_rank, num_sr);
    
    
    for r_idx = 1:num_rank
        rn_ratio = rn_ratio_set(r_idx);
        rrank = min(int64(D * rn_ratio), fix(V*(1-outlier_ratio)*0.9));
        [X_gt, U_gt, allpoints_class_gt] = Generate_data(D, V, num_groups, rrank);

        for sr_idx = 1:num_sr
            seed = randi(100000);
            rng(seed)
            shuffled_ratio = shuffled_ratio_set(sr_idx);
            
            disp(repmat('==', 1, 52))
            fprintf('Trial: %d/%d, D: %d, V: %d, Groups: %d, Rank: %d, Shuffled Ratio: %.2f, Outlier Ratio: %.2f, seed: %d\n', ...
                            (g_idx-1)*num_rank*num_sr+(r_idx-1)*num_sr+sr_idx, length(num_groups_set)*num_rank*num_sr, D, V, num_groups, rrank, shuffled_ratio, outlier_ratio, seed);
            disp(repmat('==', 1, 52))
            
            M_gt = reshape(X_gt, [D, V, num_groups]);
            [X_tilde, outliers_ID, inliers_ID] = Generate_observed_noisy_data(M_gt, outlier_ratio, shuffled_ratio);
            
            %% Ground Truth Recovery and Assessment
            fprintf('\tPart1: Ground Truth Recovery and Assessment\n');
            [X_solved_gtbases, outlier_class_solved_gtbases] = Solve_partial_LSR_(X_tilde, outliers_ID, U_gt, broadcast);
            allpoints_class_solved_gtbases = zeros(V * num_groups, 1);
            for ii = 1:V * num_groups
                if (ismember(ii, outliers_ID) == 1)
                    allpoints_class_solved_gtbases(ii) = outlier_class_solved_gtbases(outliers_ID == ii);
                end
            end
            missrate_out_gtbases = sum(allpoints_class_solved_gtbases(outliers_ID) ~= allpoints_class_gt(outliers_ID)) / length(outliers_ID);
            missrate_out_gtbases_mat(r_idx, sr_idx) = missrate_out_gtbases;
            fprintf('\t\tmissrate_out with ground truth bases = %0.4f\n', missrate_out_gtbases);
            [recover_err_gtbases, recover_err_gtbases_detectedpart] = EvaluateRefined(X_solved_gtbases, X_gt, outliers_ID, outliers_ID, num_groups, rrank);  % Noted that a same number of outlier in each group is needed
            recover_error_gtbases_mat(r_idx, sr_idx) = recover_err_gtbases;
            fprintf('\t\trecover error with ground truth bases = %0.4f\n', recover_err_gtbases);
            
            
           %% Outlier Detection
            fprintf('\tPart2: Reconstructed Outlier Detection\n');
            alpha = 5;
            lambda = 0.95; 
            thres = 5e-4;
            [Inliers_id, Outliers_id] = outlier_detection_with_outlier_num_known(X_tilde, length(outliers_ID), alpha, lambda);
            remark(Inliers_id, Outliers_id, inliers_ID, outliers_ID);       
            
            
           %% SSC Subspace Clustering and Inliers Classification
            fprintf('\tPart3: Reconstructed SSC Subspace Clustering and Inliers Classification\n');
            Y = X_tilde(:, Inliers_id);
            t = allpoints_class_gt(Inliers_id);
            r = 0; affine = false; rho = 1; alpha = 20; outlier = true;
            [missrate_in, ggrps] = SSC(Y, r, affine, alpha, outlier, rho, t, false);
            missrate_in_mat(r_idx, sr_idx) = missrate_in;
            fprintf(['\t\tmissrate_in = ' num2str(missrate_in) '\n']);
            
           %% Reconstructed Recovery and Assessment
            fprintf('\tPart3: Reconstructed Recovery and Assessment\n');
            label_hat = zeros(length(Inliers_id), 2);
            X_selected_inlier_groups = cell(1, num_groups);
            index = 1;
            for k = 1:num_groups
                kth_group_ids = Inliers_id(ggrps == k);
                X_selected_inlier_groups{k} = X_tilde(:, kth_group_ids); 
                num_kth_group_ids = length(kth_group_ids);
                label_hat(index:index + num_kth_group_ids - 1, :) = [kth_group_ids', k * ones(num_kth_group_ids, 1)];
                index = index + num_kth_group_ids;
            end
            label_hat = sortrows(label_hat);
            inlier_class_solved = label_hat(:, 2);
            
            U_esti = zeros(D, rrank, num_groups);
            for k = 1:num_groups
                [U, ~, ~] = svd(X_selected_inlier_groups{k});
                U_esti(:, :, k) = U(:, 1:rrank);
            end
            
            % Estimated Basis Recover
            [X_solved_esti, outlier_class_solved] = Solve_partial_LSR_(X_tilde, Outliers_id, U_esti, broadcast);
            allpoints_class_solved = zeros(size(X_tilde, 2), 1); 
            for ii = 1:size(X_tilde, 2)
                if (ismember(ii, Inliers_id) == 1)
                    allpoints_class_solved(ii, 1) = inlier_class_solved(Inliers_id == ii);
                elseif (ismember(ii, Outliers_id) == 1)
                    allpoints_class_solved(ii, 1) = outlier_class_solved(Outliers_id == ii);
                end
            end
            undetected_outlier_num = length(setdiff(outliers_ID, Outliers_id));
            detected_true_Outliers_id = setdiff(Outliers_id, inliers_ID);
            missclassified_outlier_num = sum(allpoints_class_solved(detected_true_Outliers_id) ~= allpoints_class_gt(detected_true_Outliers_id));
            missrate_out_reconstructed = (missclassified_outlier_num) / (length(outliers_ID) - undetected_outlier_num);
            undetected_outlier_num_ratio = undetected_outlier_num / length(outliers_ID);
            missrate_out_mat(r_idx, sr_idx) = missrate_out_reconstructed;
            undetected_outlier_ratio_mat(r_idx, sr_idx) = undetected_outlier_num_ratio;
            fprintf('\t\tmissrate_out with reconstructed bases = %0.4f, \n\t\t\twith undetected_outlier_num_ratio = %.4f\n', missrate_out_reconstructed, undetected_outlier_num_ratio);
            [recover_err, recover_err_detectedpart] = EvaluateRefined(X_solved_esti, X_gt, outliers_ID, Outliers_id, num_groups, rrank); % Noted that a same number of outlier in each group is needed
            recover_error_mat(r_idx, sr_idx) = recover_err;
            recover_error_detectedpart_mat(r_idx, sr_idx) = recover_err_detectedpart;
            fprintf('\t\trecover error with reconstructed bases = %0.4f, \n\t\t\twhere successfully detected part recover error = %.4f\n\n', recover_err, recover_err_detectedpart);
        end
    end
    
    %% the RESULT
    if store
        valid_shuffled_ratio_set = arrayfun(@(x) matlab.lang.makeValidName(sprintf('SRatio_%.1f', x)), shuffled_ratio_set, 'UniformOutput', false);
        valid_rn_ratio_set = arrayfun(@(x) matlab.lang.makeValidName(sprintf('RRatio_%.2f', x)), rn_ratio_set, 'UniformOutput', false);

        RESULT.missrate_out_allgroup{g_idx} = array2table(missrate_out_mat,...
            'VariableNames', valid_shuffled_ratio_set, ...
            'RowNames', valid_rn_ratio_set);

        RESULT.missrate_out_gtbases_allgroup{g_idx} = array2table(missrate_out_gtbases_mat, ...
            'VariableNames', valid_shuffled_ratio_set, ...
            'RowNames', valid_rn_ratio_set);

        RESULT.recover_error_allgroup{g_idx} = array2table(recover_error_mat, ...
            'VariableNames', valid_shuffled_ratio_set, ...
            'RowNames', valid_rn_ratio_set);

        RESULT.recover_error_gtbases_allgroup{g_idx} = array2table(recover_error_gtbases_mat,...
            'VariableNames', valid_shuffled_ratio_set, ...
            'RowNames', valid_rn_ratio_set);

        RESULT.recover_error_detectedpart_allgroup{g_idx} = array2table(recover_error_detectedpart_mat,...
            'VariableNames', valid_shuffled_ratio_set, ...
            'RowNames', valid_rn_ratio_set);

        RESULT.missrate_in_allgroup{g_idx} = array2table(missrate_in_mat,...
            'VariableNames', valid_shuffled_ratio_set, ...
            'RowNames', valid_rn_ratio_set);

        RESULT.undetected_outlier_ratio_allgroup{g_idx} = array2table(undetected_outlier_ratio_mat,...
            'VariableNames', valid_shuffled_ratio_set, ...
            'RowNames', valid_rn_ratio_set);
    end
end
if use_diary
    diary off;
end
save('RESULT_Group.mat', 'RESULT');
run DataRecovery_multisubspace_plot.m