%% Basis reconstrution
global DATA LABEL BASES INIT PIC
tic;
% Parameters for SSC
fprintf('Basis reconstrution...\n')
fprintf('\trank of each estimated space = %d\n', INIT.rrank)
fprintf('\t-->Now run SSC algorithm<--\n')
% fprintf('\t\talpha = %d\n', alpha)
r = 0; affine = false; outlier = true; rho = 1; 
if DATA.num_groups<8
    alpha = 25;
else
    alpha = 60;
end
Y = DATA.X_tilde(:, LABEL.Inliers_id);
t = LABEL.allpoints_class_gt(LABEL.Inliers_id);
[missrate_in, ggrps] = SSC(Y, r, affine, alpha, outlier, rho, t);
fprintf(['\t-->missrate_in = ' num2str(missrate_in) '<--\n']);

% 设定实验算法的基本参数
num_inliers = length(LABEL.Inliers_id);
label_hat = zeros(num_inliers, 2);
X_selected_inlier_groups = cell(1, DATA.num_groups); 
%X_selected_inlier_groups的顺序决定Basis的顺序，进而决定每个solved outlier class
%的取值方式，而solved outlier class需要与gt outlier class做比较，而gt outlier class
%取值来源于load data时候每个subspace的顺序，故而若要最后的评估有意义，
%需要X_selected_inlier_groups的顺序依次沿着load data时候每个subspace的顺序来select，
%这便需要下述循环中ggrps == k所代表的确实是第k个subspace的数据。

index = 1;
for k = 1:DATA.num_groups
    kth_group_ids = LABEL.Inliers_id(ggrps == k);
    X_selected_inlier_groups{k} = DATA.X_tilde(:, kth_group_ids); 
    
    % The following procedure is designed for the final assessment.
    num_kth_group_ids = length(kth_group_ids);
    label_hat(index:index+num_kth_group_ids-1, :) = [kth_group_ids', k*ones(num_kth_group_ids, 1)];
    
    index = index + num_kth_group_ids;
end
DATA.X_selected_inlier_groups = X_selected_inlier_groups;
Q = sortrows(label_hat);
% LABEL.inlier_class_solved = Q(:, 2);
LABEL.inlier_class_solved = Q;

U_esti = zeros(PIC.height*PIC.width, INIT.rrank, DATA.num_groups);
for k = 1:DATA.num_groups
    [U, ~, ~] = svd(DATA.X_selected_inlier_groups{k});
    U_esti(:, :, k) = U(:, 1:INIT.rrank);
end
BASES.U_esti = U_esti;
fprintf('\tOver, costing %0.2fs\n\n', toc)



