%% Inliers Detection
global DATA LABEL

fprintf('Inliers detection...\n')

% SSC_outliers_detection
% alpha = 1e10;
% detect_flag = 2; %---%
% [Inliers_id, Outliers_id] = Outlier_Detect_v2(DATA.X_tilde, detect_flag, alpha);

tic;
[LABEL.Inliers_id, LABEL.Outliers_id] = outlier_detection(DATA.X_tilde);
remark(LABEL.Inliers_id, LABEL.Outliers_id, LABEL.inliers_ID, LABEL.outliers_ID);
fprintf('\t...Over, took %.3fs.\n\n', toc)