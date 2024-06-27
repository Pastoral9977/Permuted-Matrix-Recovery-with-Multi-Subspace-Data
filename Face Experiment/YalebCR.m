close all
clear
%% ��������
% ����·��
addpath(genpath('..\Copy_SSC_ADMM_v1.1'));
addpath(genpath('..\Functions'));
addpath(genpath(pwd));
rng(21)
% ��������
fprintf('Data Loading...')
tic
load('YaleB_cell'); % ���غ��YaleB_cell����Ϊ1*38��cell�������ͣ���38���ˣ�ÿ��cellΪͬһ�����Ķ��źڰ���Ƭ������Ϊ(ÿ������ͼ������һά����ά��*��Ƭ����)
h = 192; % ͼƬ�ĸ�
w = 168; % ͼƬ�Ŀ�
fprintf('Over, costing %0.2fs\n\n', toc)

%% ��ʼ������
% �趨ʵ�����ݵĻ�������
num_groups = 10; % һ��num_groups��������num_groups���ӿռ䣩
start_person_id = 1; % �ӵڼ����˿�ʼѡ
n = [];
for i = 1:num_groups
    [~, n(i)] = size(YaleB_cell{start_person_id-1+i});
end
N = sum(n); % ͼƬ���� 
outliers_num = 16; % number of outliers out of 64 points

% ѹ��������Ĳ���
hh = 48; % ѹ����ͼƬ�ĸ�
ww = 42; % ѹ����ͼƬ�Ŀ�
v_patches = 48; % v_patches: number of vertical blocks in one 'big column', #rows
h_patches = 42; % h_patches: number of horizontal blocks in one 'big row', #cols
perm_flag = 1; % 0 for fully shuffling; 1 for partially shuffling.
X_gt = zeros(hh*ww, N); % ��ʼ��ѹ��������ǰ����������
X_tilde = zeros(hh*ww, N); % ��ʼ��ѹ������������������

%% ���죨ѹ��+��������
tic;
fprintf('Data Construction(downsampling & corrpution)...')
nn = [0, n];
outliers_ID = [];
for i = 1:num_groups
    X_gt(:, sum(nn(1:i))+1:sum(nn(1:i+1))) = DownsamplePiciture(YaleB_cell{1, start_person_id-1+i}, h, w, hh, ww);
    
    outliers_id = randperm(n(i), outliers_num);
%     if length(setdiff(outliers_id, 46)) == outliers_num
%         outliers_id(1) = 46;
%     end
%     if length(setdiff(outliers_id, 23)) == outliers_num
%         outliers_id(2) = 23;
%     end
%     if length(setdiff(outliers_id, 54)) == outliers_num
%         outliers_id(3) = 54;
%     end
%     if length(setdiff(outliers_id, 12)) == outliers_num
%         outliers_id(4) = 12;
%     end

    % ȷ��ĳЩ�ض�ID������outliers_id��
    specific_ids = [46, 23, 54, 12];
    for j = 1:length(specific_ids)
        if length(setdiff(outliers_id, specific_ids(j))) == outliers_num
            outliers_id(j) = specific_ids(j);
        end
    end

    [X] = permute_face_randomly(X_gt(:, sum(nn(1:i))+1:sum(nn(1:i+1))), outliers_id, hh, ww, v_patches, h_patches, perm_flag);
    X_tilde(:, sum(nn(1:i))+1:sum(nn(1:i+1))) = X;
    outliers_ID = [outliers_ID, sort(outliers_id)+sum(nn(1:i))];
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
% �趨ʵ���㷨�Ļ�������
rrank = 8; % low-rank�㷨����
fprintf('\n\trank of each estimated space = %d', rrank)
[Inliers_id, Outliers_id] = Outlier_Detect_v3(X_tilde);
remark(Inliers_id, Outliers_id, inliers_ID, outliers_ID);

% Parameters for SSC
allpoints_class_gt = zeros(N, 1);
for i = 1:num_groups
    allpoints_class_gt(sum(nn(1:i))+1:sum(nn(1:i+1))) = i;
end

r = 0; affine = false; outlier = true; rho = 1; alpha = 60;
fprintf('\t-->Now run SSC algorithm<--\n')
[missrate_in, ggrps] = SSC(X_tilde(:, Inliers_id), r, affine, alpha, outlier, rho, allpoints_class_gt(Inliers_id));
fprintf(['\t-->missrate_in = ' num2str(missrate_in) '<--\n']);

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
for i = 1:size(X_tilde, 2)
    % Inliers_id��inlier_class_solved������ͬ��ǰ����id�������Ƕ�Ӧ��class��
    % Inliers_id == i�ҳ�i����Inliers_id��λ�ã��ٴ���inlier_class_solved��
    % �ҳ�i����Ӧ��class������Outliers_id��outlier_class_solvedͬ��ǰ��
    % ��SSC����õ�class�����������з����������
    if (ismember(i, Inliers_id) == 1)
        allpoints_class_solved(i) = inlier_class_solved(Inliers_id == i); 
    elseif (ismember(i, Outliers_id) == 1)
        allpoints_class_solved(i) = outlier_class_solved(Outliers_id == i);
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

close all;
run 'FigureFace'
