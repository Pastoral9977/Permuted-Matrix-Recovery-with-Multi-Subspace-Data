if (shuffled_ratio >= 0.6)
    [X_solved,~, outlier_class_solved] = Solve_all(X_tilde, Outliers_id, U_esti);
else
    [X_solved,~, outlier_class_solved] = Solve_partial_LSR(X_tilde, Outliers_id, U_esti);
end

allpoints_class_solved = zeros(N,1);
for p = 1:size(X_tilde,2)
    if (ismember(p,Inliers_id) == 1)
        allpoints_class_solved(p,1) = inlier_class_solved(Inliers_id == p);
    elseif (ismember(p,Outliers_id) == 1)
        allpoints_class_solved(p,1) = outlier_class_solved(Outliers_id == p);
    end
end  

% Calculate missrate for reconstructed bases
undetected_outlier_num = length(setdiff(outliers_ID, Outliers_id));
detected_true_Outliers_id = setdiff(Outliers_id, inliers_ID);
missclassified_outlier_num = sum(allpoints_class_solved(detected_true_Outliers_id) ~= s(detected_true_Outliers_id));
missrate_out = (missclassified_outlier_num+undetected_outlier_num) / length(outliers_ID);
undetected_outlier_num_ratio = undetected_outlier_num / length(outliers_ID);
[err_ratio_out] = Evaluate(X_solved(:,outliers_ID),X_gt(:,outliers_ID));

fprintf('\tmissrate_out  = %0.4f, where undetected_outlier_num_ratio = %.4f\n', missrate_out, undetected_outlier_num_ratio);
fprintf('\terr_ratio_out = %.4f\n\n', err_ratio_out);
