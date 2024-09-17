function [err_ratio, err_ratio_detectedout] = EvaluateRefined(X_solved, X_gt,  outliers_ID, Outliers_id, num_groups, rrank)

[D,N] = size(X_solved);
V = N/num_groups;

num_outlier_each = length(outliers_ID)/num_groups;
detected_ids = intersect(outliers_ID, Outliers_id);
M_gt = reshape(X_solved, [D, V, num_groups]);
M_solved_allout = reshape(X_solved(:, outliers_ID), [D, num_outlier_each, num_groups]);

for i = 1:num_groups
    [U,~,~] = svd(M_gt(:,:,i),'econ');
    B = U(:,1:rrank);
    P = B*B';
    A = P * M_solved_allout(:,:,i);
    M_solved_allout(:,:,i) = A ./ vecnorm(A);
end
X_solved_refine = reshape(M_solved_allout, [D, length(outliers_ID)]);
X_solved_refine_detectedout = X_solved_refine(:, ismember(outliers_ID, detected_ids));
X_gt_refine = X_gt(:, outliers_ID);
X_gt_refine_detectedout = X_gt(:, detected_ids);

err_refine = norm(X_solved_refine - X_gt_refine, 'fro');
err_ratio = err_refine / norm(X_gt_refine, 'fro');
err_refine_detectedout = norm(X_solved_refine_detectedout - X_gt_refine_detectedout, 'fro');
err_ratio_detectedout = err_refine_detectedout / norm(X_gt_refine_detectedout, 'fro');

end
