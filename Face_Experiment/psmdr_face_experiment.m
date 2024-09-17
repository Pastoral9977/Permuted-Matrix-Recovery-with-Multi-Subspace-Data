close all
clear
%% load data
addpath(genpath('..\Functions'));
addpath(genpath(pwd));
fprintf('Data Loading...')
tic
load('YaleB_cell');
h = 192; 
w = 168; 
fprintf('Over, costing %0.2fs\n\n', toc)
missing_ratio  = 0.0;
shuffled_ratio = 0.4;
% seed = randi(10000000)
% seed=21;
% seed = 7925820;6138427;
seed = 6218348;
rng(seed)
%% initialization
num_groups = 10; 
start_person_id = 1; 
n = [];
for ii = 1:num_groups
    [~, n(ii)] = size(YaleB_cell{start_person_id-1+ii});
end
N = sum(n); 
outliers_num = 16; 
id = 23;

hh = 48; % height
ww = 42; % width
v_patches = 8; % v_patches: number of vertical blocks in one 'big column', #rows
h_patches = 7; % h_patches: number of horizontal blocks in one 'big row', #cols
perm_flag = 1; % 0 for fully shuffling; 1 for partially shuffling.
X_gt = zeros(hh*ww, N); 
X_tilde = zeros(hh*ww, N); 

%% corruption
tic;
fprintf('Data Construction(downsampling & corrpution)...')
nn = [0, n];
outliers_ID = [];
for ii = 1:num_groups
    X_gt(:, sum(nn(1:ii))+1:sum(nn(1:ii+1))) = DownsamplePicture(YaleB_cell{start_person_id-1+ii}, h, w, hh, ww);
%     outliers_id = randperm(n(ii), outliers_num);
%     specific_ids = [46, 23, 54, 12];
%     for j = 1:length(specific_ids)
%         if length(setdiff(outliers_id, specific_ids(j))) == outliers_num
%             outliers_id(j) = specific_ids(j);
%         end
%     end
    outliers_id = [1:4:60, id];

    [X] = permute_face_randomly(X_gt(:, sum(nn(1:ii))+1:sum(nn(1:ii+1))), outliers_id, hh, ww, v_patches, h_patches, missing_ratio, shuffled_ratio, seed*ii);
    X_tilde(:, sum(nn(1:ii))+1:sum(nn(1:ii+1))) = X;
    outliers_ID = [outliers_ID, sort(outliers_id)+sum(nn(1:ii))];
end
inliers_ID = setdiff(1:N, outliers_ID);
fprintf('Over, costing %0.2fs\n', toc)
fprintf('\t-->Data Information:<--\n')
    fprintf(...
    ['\tnum_groups = ' num2str(num_groups) '\n\tnumber of total samples = ' num2str(size(X_tilde,2)) ...
    '\n\tdimension of each sample = ' num2str(size(X_tilde,1)) '\n\tnumber of outliers = ' num2str(outliers_num*num_groups) ...
    '\n\toutlier_mod = ' num2str(perm_flag) '\n\n']...
    )

%% Basis reconstrution
fprintf('Basis reconstrution...')
rrank = 8; 
fprintf('\n\trank of each estimated space = %d\n', rrank)
[Inliers_id, Outliers_id] = outlier_detection(X_tilde);
remark(Inliers_id, Outliers_id, inliers_ID, outliers_ID);

% Parameters for SSC
allpoints_class_gt = zeros(N, 1);
for ii = 1:num_groups
    allpoints_class_gt(sum(nn(1:ii))+1:sum(nn(1:ii+1))) = ii;
end

r = 0; affine = false; outlier = true; rho = 1; alpha = 50;
fprintf('\t-->Now run SSC algorithm<--\n')
[missrate_in, ggrps] = SSC(X_tilde(:, Inliers_id), r, affine, alpha, outlier, rho, allpoints_class_gt(Inliers_id));
fprintf(['\t-->missrate_in = ' num2str(missrate_in) '<--\n']);

X_gt = X_gt./max(X_gt);
X_tilde = X_tilde./max(X_tilde);
label_hat = zeros(length(Inliers_id), 2);
X_selected_inlier_groups = cell(1, num_groups);
index = 1;
for k = 1:num_groups
    kth_group_ids = Inliers_id(ggrps == k);
    X_selected_inlier_groups{k} = X_tilde(:, kth_group_ids); 
    
    % The following procedure is designed for the final assessment
    num_kth_group_ids = length(kth_group_ids);
    label_hat(index:index+num_kth_group_ids-1, :) = [kth_group_ids', k*ones(num_kth_group_ids, 1)];
    
    index = index + num_kth_group_ids;
end
label_hat = sortrows(label_hat);
inlier_class_solved = label_hat(:, 2);

U_esti = zeros(hh*ww, rrank, num_groups);
for k = 1:num_groups
    [U, ~, ~] = svd(X_selected_inlier_groups{k});
    U_esti(:, :, k) = U(:, 1:rrank);
end
fprintf('\n')

%% Matrix reconstrution
tic
fprintf('Matrix reconstrution...\n')
[X_solved, outlier_class_solved] = Solve_partial_LSR(X_tilde, Outliers_id, U_esti);
fprintf('\tOver, costing %0.2fs\n', toc)

%% Final Assesment for outlier belongings
allpoints_class_solved = zeros(sum(n), 1);
for ii = 1:size(X_tilde, 2)
    % Inliers_id与inlier_class_solved长度相同，前者是id，后者是对应的class，
    % Inliers_id == i找出i所在Inliers_id的位置，再带入inlier_class_solved中
    % 找出i所对应的class。下述Outliers_id与outlier_class_solved同理。前者
    % 由SSC聚类得到class，后者由自研法方法解决。
    if (ismember(ii, Inliers_id) == 1)
        allpoints_class_solved(ii) = inlier_class_solved(Inliers_id == ii); 
    elseif (ismember(ii, Outliers_id) == 1)
        allpoints_class_solved(ii) = outlier_class_solved(Outliers_id == ii);
    end
end

MISSRATE = sum(allpoints_class_solved(:) ~= allpoints_class_gt(:)) / N;

undetected_outlier_num = length(setdiff(outliers_ID, Outliers_id));
detected_true_Outliers_id = setdiff(Outliers_id, inliers_ID);
missclassified_outlier_num = sum(allpoints_class_solved(detected_true_Outliers_id) ~= allpoints_class_gt(detected_true_Outliers_id));
missrate_out_reconstructed = (missclassified_outlier_num) / (length(outliers_ID)-undetected_outlier_num);
undetected_outlier_num_ratio = undetected_outlier_num/length(outliers_ID);
fprintf('\tmissrate_out with reconstructed bases = %0.4f, where undetected_outlier_num_ratio = %.4f\n\n', missrate_out_reconstructed, undetected_outlier_num_ratio);

[err_ratio_out] = Evaluate(X_solved(:, outliers_ID), X_tilde(:, outliers_ID));

%% figures
close all;
g = figure;
g.Units = 'centimeters';
shrink = 0.9;
g.Position = [4 3.1 3.5*num_groups*shrink 13*shrink];
rows = 3;
for ii  = 1:num_groups
    id = 23+(ii-1)*64;
    img = reshape(X_gt(:, id),hh, ww);
    perm_img = reshape(X_tilde(:, id),hh, ww);
    reco_img = reshape(X_solved(:,id),hh, ww);
    img = img/max(img, [], 'all');
    perm_img = perm_img/max(perm_img, [], 'all');
    reco_img = reco_img/max(reco_img, [], 'all');

    subplot(rows,num_groups,ii)
    imshow(img)
    title('Original')

    subplot(rows,num_groups,ii+num_groups)
    imshow(perm_img)
    title('Corrupted')

    subplot(rows,num_groups,ii+2*num_groups)
    imshow(reco_img)
    title('FLSRF')
end
% run 'FigureFace'
