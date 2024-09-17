%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi Subspace Face Recovery
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clear
close all

use_SSC = true;
addpath(genpath('..\Functions'));
addpath(genpath(pwd));

%% load images
use_downsampled_version = true;
if ~use_downsampled_version
    suffix = '';
    n = 192;
    m = 168;
    patch = 6;
else
    suffix = '_light';
    n = 48;
    m = 42;
    np = 8;
    mp = 7;
    patch = 1;
    missing_ratio = 0.0;
    shuffled_ratio = 0.4; 
    seed = 6218348;
end
a = 8;
b = true;
normalization = false;

load('YaleB_cell.mat')
images = YaleB_cell{1};
images = DownsamplePicture(images, 192, 168, n, m);

%% originate images example
fprintf('images\n')
fprintf("\trank = %d, size = (%d %d)\n", rank(images), size(images))
X_gt = images;
id = 23;

%% permute images example
outliers_num = 16;
specific_ids = [46, 23, 54, 12];
outliers_id = randperm(size(images,2), outliers_num);
for j = 1:length(specific_ids)
    if length(setdiff(outliers_id, specific_ids(j))) == outliers_num
        outliers_id(j) = specific_ids(j);
    end
end
perm_images = permute_face_randomly(images, outliers_id, n, m, np, mp, missing_ratio, shuffled_ratio, seed);
fprintf('perm_images\n')
fprintf("\trank = %d, size = (%d %d)\n", rank(perm_images), size(perm_images))
X_tilde = perm_images;

%% Robust PCA Recovery
ids = 2:10; 
num_groups = length(ids)+1;
for ii = 1:length(ids)
    id = ids(ii);
    load('YaleB_cell.mat')
    images = YaleB_cell{id};
    images = DownsamplePicture(images, 192, 168, n, m);
    outliers_id = randperm(size(images,2), outliers_num);
    specific_ids = [46, 23, 54, 12];
    for j = 1:length(specific_ids)
        if length(setdiff(outliers_id, specific_ids(j))) == outliers_num
            outliers_id(j) = specific_ids(j);
        end
    end
    perm_images = permute_face_randomly(images, outliers_id, n, m, np, mp, missing_ratio, shuffled_ratio, seed*(ii+1));
    perm_img = reshape(perm_images(:, id), n, m);
    X_gt = [X_gt, images];
    X_tilde = [X_tilde, perm_images];
end
fprintf('X_gt\n')
fprintf("\trank = %d, size = (%d %d)\n", rank(X_gt), size(X_gt))
fprintf('X_tilde\n')
fprintf("\trank = %d, size = (%d %d)\n", rank(X_tilde), size(X_tilde))

%% Single SSC Recovery
affine = b; 
alpha = a;
if ~normalization
    CMat = admmOutlier_mat_func(X_tilde,affine,alpha);
    N = size(X_tilde,2);
    C = CMat(1:N,:);
    X_SSC = X_tilde*C;
else
    Vnorms = vecnorm(X_tilde);
    CMat = admmOutlier_mat_func(X_tilde./Vnorms, affine, alpha);
    N = size(X_tilde,2);
    C = CMat(1:N,:);
    X_SSC = (X_tilde*C).*Vnorms;
end

reco_img = reshape(X_SSC(:,23),n,m);
fprintf('PMSDR recovered version: X_SSC\n')
fprintf("\trank = %d\n", rank(X_SSC))

%% SSC + PMSDR 4-stage Pipeline
if use_SSC
    addpath(genpath('..\Functions'));
    addpath(genpath(pwd));

    [Inliers_id, Outliers_id] = outlier_detection(X_tilde);

    r = 0; affine = false; outlier = true; rho = 1; alpha = 50;
    allpoints_class_gt = [];
    for ii = 1:num_groups
        allpoints_class_gt = [allpoints_class_gt; ones(64,1)*ii];
    end
    [missrate_in, ggrps] = SSC(X_tilde(:, Inliers_id), r, affine, alpha, outlier, rho, allpoints_class_gt(Inliers_id));

    label_hat = zeros(length(Inliers_id), 2);
    X_selected_inlier_groups = cell(1, num_groups);
    index = 1;
    for k = 1:num_groups
        kth_group_ids = Inliers_id(ggrps == k);
        X_selected_inlier_groups{k} = X_tilde(:, kth_group_ids); 
        num_kth_group_ids = length(kth_group_ids);
        label_hat(index:index+num_kth_group_ids-1, :) = [kth_group_ids', k*ones(num_kth_group_ids, 1)];
        index = index + num_kth_group_ids;
    end
    label_hat = sortrows(label_hat);
    inlier_class_solved = label_hat(:, 2);
    
    rrank = 8;
    U_esti = zeros(n*m, rrank, num_groups);
    for k = 1:num_groups
        [U, ~, ~] = svd(X_selected_inlier_groups{k});
        U_esti(:, :, k) = U(:, 1:rrank);
    end

    [X_PMSDR, outlier_class_solved] = Solve_partial_LSR(X_tilde, Outliers_id, U_esti);
    fprintf('PMSDR recovered version: X_PMSDR\n')
    fprintf("\trank = %d\n", rank(X_PMSDR))
    N = size(X_tilde,2);
    allpoints_class_solved = zeros(N, 1);
    for i = 1:N
        if (ismember(i, Inliers_id) == 1)
            allpoints_class_solved(i) = inlier_class_solved(Inliers_id == i); 
        elseif (ismember(i, Outliers_id) == 1)
            allpoints_class_solved(i) = outlier_class_solved(Outliers_id == i);
        end
    end

    X_PMSDR_SSC = zeros(size(X_tilde));
    for id = 1:num_groups
        XX = X_tilde(:, allpoints_class_solved==id);
        affine = b; 
        alpha = a;
        if ~normalization
            CMat = admmOutlier_mat_func(XX,affine,alpha);
            N = size(XX,2);
            C = CMat(1:N,:);
            X_PMSDR_SSC(:, allpoints_class_solved==id) = XX*C;
        else
            vnorms = vecnorm(XX);
            CMat = admmOutlier_mat_func(XX./vnorms,affine,alpha);
            N = size(XX,2);
            C = CMat(1:N,:);
            X_PMSDR_SSC(:, allpoints_class_solved==id) = (XX*C).*vnorms;
        end
    end
    fprintf('PMSDR_SSC recovered version: X_PMSDR_SSC\n')
    fprintf("\trank = %d\n", rank(X_PMSDR_SSC))
end

run figure_plot_ssc












