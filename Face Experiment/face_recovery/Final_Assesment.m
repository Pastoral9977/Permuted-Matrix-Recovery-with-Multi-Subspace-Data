%% Final Assesment for outlier belongings
allpoints_class_solved = zeros(sum(n), 1);
for i = 1:size(X_tilde, 2)
    if (ismember(i, Inliers_id) == 1)
        allpoints_class_solved(i) = inlier_class_solved(Inliers_id == i);
    elseif (ismember(i, Outliers_id) == 1)
        allpoints_class_solved(i) = outlier_class_solved(Outliers_id == i);
    end
end

undetected_outlier_num = length(setdiff(outliers_ID, Outliers_id));
detected_true_Outliers_id = setdiff(Outliers_id, inliers_ID);
missclassified_outlier_num = sum(allpoints_class_solved(detected_true_Outliers_id) ~= allpoints_class_gt(detected_true_Outliers_id));
missrate_out = (missclassified_outlier_num+undetected_outlier_num) / length(outliers_ID);
undetected_outlier_num_ratio = undetected_outlier_num/length(outliers_ID);
fprintf('\tmissrate_out with reconstructed bases = %0.4f, where undetected_outlier_num_ratio = %.4f\n\n', missrate_out, undetected_outlier_num_ratio);

MISSRATE = sum(allpoints_class_solved(:) ~= allpoints_class_gt(:)) / N;
[err_ratio_out] = Evaluate(X_solved(:, outliers_ID), X_tilde(:, outliers_ID));