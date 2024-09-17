close all
clear
addpath(genpath('../../functions'))
addpath(genpath(pwd))

D = 20;
V = 50;
num_settings = 5;
num_groups_set = [2, 2, 3, 5, 5];
parts_set      = [2, 3, 3, 3, 5];
shuffled_ratio = 0.5;
outlier_ratio = 0.6;
involve_MRUC = true;
involve_L_PMSDR_MRUC = false;
involve_R_PMSDR_MRUC = false;
involve_PMSDR_MRUC = true;
use_diary = true;

missing_percent=0; % Percentage of missing entries
eps_init=0.1;
eps_decay=0.5;
max_out_iter=50000;
max_in_iter=10000;
verbose=true;
lambda_list= 0.55;
rr = 2;
snr=0;

if use_diary
    diary('DataRecovery_with_MRUC.txt')
end
for setting_id = 1:5
num_groups = num_groups_set(setting_id);
parts = parts_set(setting_id);

gt_num_groups = num_groups;
gt_rr = rr;
% gt_num_groups = 1;
% gt_rr = 6;
% rr = gt_rr/num_groups;
% V = 2*V;

fprintf([repmat('==', 1, 53), '\n'])
fprintf([repmat('==', 1, 53), '\n'])
%% data
% seed = randi(9999999);
seed = 1068839;
fprintf('seed = %d\ngt_num_groups = %d\ngt_rr = %d\n', seed, gt_num_groups, gt_rr)
rng(seed)
[B1s, A1s, Ans, A1_row_inds, An_row_inds, oracle_images, permutation_matrices, test_inds, test_labels, norm_constants, rawimages, ...
          B1_TOTAL, A1_TOTAL, An_TOTAL, A1_row_ind_TOTAL, An_row_ind_TOTAL, permutation_matrix, norm_constant_TOTAL] ...
    = generate_data(D, V, gt_rr, gt_num_groups, shuffled_ratio, outlier_ratio, parts, missing_percent, snr, seed);
shuffle_columns = fix(V * outlier_ratio) * num_groups;

%% L_PMSDR MRUC
if involve_L_PMSDR_MRUC
    L_num_groups = num_groups-1;
    [B1_SUB_solved,A1_SUB_solved,A1_row_ind_SUB_solved,An_SUB_solved,An_row_ind_SUB_solved] = PMSDR_data_handling(B1_TOTAL, A1_TOTAL, An_TOTAL, A1_row_ind_TOTAL, An_row_ind_TOTAL, L_num_groups, rr);
    trials = 4;
    results=zeros(2,trials*L_num_groups);
    seedss=zeros(1,trials*L_num_groups);

    fprintf([repmat('==', 1, 54), '\n'])
    fprintf('L_PMSDR_MRUC: D = %d, V = %d, L_num_groups = %d, shuffle_columns_per_group = %d, parts = %d, shuffled_ratio = %.2f\n',D,V,L_num_groups,shuffle_columns/num_groups,parts,shuffled_ratio)
    fprintf([repmat('==', 1, 54), '\n'])
for t_id = 1:trials
    for jj = 1:L_num_groups
        ii = (t_id-1)*L_num_groups+jj;
        seedd=randi(100000000);
        rng(seedd)
        seedss(t_id)=seedd;
        fprintf('L_PMSDR_MRUC --- Trial: %d, seed:%d\n', ii, seedd)
        shuffle_columns_per_group = fix(V * outlier_ratio);
        shuffle_columns_per_group_per_parts = [];
        for kk = 1:parts
            shuffle_columns_per_group_per_parts = [shuffle_columns_per_group_per_parts,length(An_SUB_solved{jj}{1,kk})];
        end
        Bn=generate_Bn_v2(parts,D,shuffle_columns_per_group_per_parts,norm_constant_TOTAL,seedd);
        [best_P,best_B,final_cert,history,result]...
            =CD_complete_v2(B1_SUB_solved{jj},Bn,A1_SUB_solved{jj},An_SUB_solved{jj},A1_row_ind_SUB_solved{jj},An_row_ind_SUB_solved{jj},lambda_list,...
                         eps_init,eps_decay,verbose,max_out_iter,max_in_iter,100,1.5,0.000001,permutation_matrix,0,1);
        results(:,ii)=[history(1,end); history(5,end)*100/(parts*D)];
        fprintf('obj:%.5f, Perr_ratio:%.2f, seed:%d\n------\n------\n',results(1,ii),results(2,ii), seedd)
    end
end
Perr_ratios_for_PMSDR_MRUC = results(2,:);
mean_for_PMSDR_MRUC = mean(Perr_ratios_for_PMSDR_MRUC);
min_for_PMSDR_MRUC  = min(Perr_ratios_for_PMSDR_MRUC);
max_for_PMSDR_MRUC  = max(Perr_ratios_for_PMSDR_MRUC);
std_for_PMSDR_MRUC  = std(Perr_ratios_for_PMSDR_MRUC);
fprintf('L_PMSDR_MRUC: mean Perr = %.4f, std Perr = %.4f, min Perr = %.4f, max Perr = %.4f\n--------\n--------\n',mean_for_PMSDR_MRUC,std_for_PMSDR_MRUC,min_for_PMSDR_MRUC,max_for_PMSDR_MRUC)
end

%% PMSDR MRUC
if involve_PMSDR_MRUC
[B1_SUB_solved,A1_SUB_solved,A1_row_ind_SUB_solved,An_SUB_solved,An_row_ind_SUB_solved] = PMSDR_data_handling(B1_TOTAL, A1_TOTAL, An_TOTAL, A1_row_ind_TOTAL, An_row_ind_TOTAL, num_groups, rr);
trials = 4;
results=zeros(2,trials*num_groups);
seedss=zeros(1,trials*num_groups);

fprintf([repmat('==', 1, 54), '\n'])
fprintf('PMSDR_MRUC: D = %d, V = %d, num_groups = %d, shuffle_columns_per_group = %d, parts = %d, shuffled_ratio = %.2f\n',D,V,num_groups,shuffle_columns/num_groups,parts,shuffled_ratio)
fprintf([repmat('==', 1, 54), '\n'])
for t_id = 1:trials
    for jj = 1:num_groups
        ii = (t_id-1)*num_groups+jj;
        seedd=randi(100000000);
        rng(seedd)
        seedss(t_id)=seedd;
        fprintf('PMSDR_MRUC --- Trial: %d, seed:%d\n', ii, seedd)
        shuffle_columns_per_group = fix(V * outlier_ratio);
        shuffle_columns_per_group_per_parts = [];
        for kk = 1:parts
            shuffle_columns_per_group_per_parts = [shuffle_columns_per_group_per_parts,length(An_SUB_solved{jj}{1,kk})];
        end
        Bn=generate_Bn(parts,D,shuffle_columns_per_group_per_parts,norm_constant_TOTAL,seedd);
        [best_P,best_B,final_cert,history,result]...
            =CD_complete(B1_SUB_solved{jj},Bn,A1_SUB_solved{jj},An_SUB_solved{jj},A1_row_ind_SUB_solved{jj},An_row_ind_SUB_solved{jj},lambda_list,...
                         eps_init,eps_decay,verbose,max_out_iter,max_in_iter,100,1.5,0.000001,permutation_matrix,0,1);
        results(:,ii)=[history(1,end); history(5,end)*100/(parts*D)];
        fprintf('obj:%.5f, Perr_ratio:%.2f, seed:%d\n------\n------\n',results(1,ii),results(2,ii), seedd)
    end
end
Perr_ratios_for_PMSDR_MRUC = results(2,:);
mean_for_PMSDR_MRUC = mean(Perr_ratios_for_PMSDR_MRUC);
min_for_PMSDR_MRUC  = min(Perr_ratios_for_PMSDR_MRUC);
max_for_PMSDR_MRUC  = max(Perr_ratios_for_PMSDR_MRUC);
std_for_PMSDR_MRUC  = std(Perr_ratios_for_PMSDR_MRUC);
fprintf('PMSDR_MRUC: mean Perr = %.4f, std Perr = %.4f, min Perr = %.4f, max Perr = %.4f\n--------\n--------\n',mean_for_PMSDR_MRUC,std_for_PMSDR_MRUC,min_for_PMSDR_MRUC,max_for_PMSDR_MRUC)
end

%% R_PMSDR MRUC
if involve_R_PMSDR_MRUC
    R_num_groups = num_groups+1;
fprintf([repmat('==', 1, 54), '\n'])
fprintf('R_PMSDR_MRUC: D = %d, V = %d, R_num_groups = %d, shuffle_columns_per_group = %d, parts = %d, shuffled_ratio = %.2f\n',D,V,R_num_groups,shuffle_columns/num_groups,parts,shuffled_ratio)
fprintf([repmat('==', 1, 54), '\n'])
[B1_SUB_solved,A1_SUB_solved,A1_row_ind_SUB_solved,An_SUB_solved,An_row_ind_SUB_solved] = PMSDR_data_handling(B1_TOTAL, A1_TOTAL, An_TOTAL, A1_row_ind_TOTAL, An_row_ind_TOTAL, R_num_groups, rr);
trials = 4;
results=zeros(2,trials*R_num_groups);
seedss=zeros(1,trials*R_num_groups);

for t_id = 1:trials
    for jj = 1:R_num_groups
        ii = (t_id-1)*R_num_groups+jj;
        seedd=randi(100000000);
        rng(seedd)
        seedss(t_id)=seedd;
        fprintf('R_PMSDR_MRUC --- Trial: %d, seed:%d\n', ii, seedd)
        shuffle_columns_per_group = fix(V * outlier_ratio);
        shuffle_columns_per_group_per_parts = [];
        for kk = 1:parts
            shuffle_columns_per_group_per_parts = [shuffle_columns_per_group_per_parts,length(An_SUB_solved{jj}{1,kk})];
        end
        Bn=generate_Bn_v2(parts,D,shuffle_columns_per_group_per_parts,norm_constant_TOTAL,seedd);
        [best_P,best_B,final_cert,history,result]...
            =CD_complete_v2(B1_SUB_solved{jj},Bn,A1_SUB_solved{jj},An_SUB_solved{jj},A1_row_ind_SUB_solved{jj},An_row_ind_SUB_solved{jj},lambda_list,...
                         eps_init,eps_decay,verbose,max_out_iter,max_in_iter,100,1.5,0.000001,permutation_matrix,0,1);
        results(:,ii)=[history(1,end); history(5,end)*100/(parts*D)];
        fprintf('obj:%.5f, Perr_ratio:%.2f, seed:%d\n------\n------\n',results(1,ii),results(2,ii), seedd)
    end
end
Perr_ratios_for_PMSDR_MRUC = results(2,:);
mean_for_PMSDR_MRUC = mean(Perr_ratios_for_PMSDR_MRUC);
min_for_PMSDR_MRUC  = min(Perr_ratios_for_PMSDR_MRUC);
max_for_PMSDR_MRUC  = max(Perr_ratios_for_PMSDR_MRUC);
std_for_PMSDR_MRUC  = std(Perr_ratios_for_PMSDR_MRUC);
fprintf('R_PMSDR_MRUC: mean Perr = %.4f, std Perr = %.4f, min Perr = %.4f, max Perr = %.4f\n--------\n--------\n',mean_for_PMSDR_MRUC,std_for_PMSDR_MRUC,min_for_PMSDR_MRUC,max_for_PMSDR_MRUC)
end

%% MRUC
if involve_MRUC
trials = 10;
results=zeros(2,trials);
seeds=zeros(1,trials);
rng(seed)

fprintf([repmat('==', 1, 52), '\n'])
fprintf('MRUC: D = %d, V = %d, num_groups = %d, shuffle_columns_per_group = %d, parts = %d, shuffled_ratio = %.2f\n',D,V,num_groups,shuffle_columns/num_groups,parts,shuffled_ratio)
fprintf([repmat('==', 1, 52), '\n'])
for t_id = 1:trials
    seedd=randi(100000000);
    rng(seedd)
    seeds(t_id)=seedd;
    fprintf('MRUC --- Trial: %d, seed:%d\n', t_id, seedd)
    Bn_TOTAL=generate_Bn_(parts,D,shuffle_columns/parts,norm_constant_TOTAL,seedd);
    [best_P,best_B,final_cert,history,result]...
        =CD_complete(B1_TOTAL,Bn_TOTAL,A1_TOTAL,An_TOTAL,A1_row_ind_TOTAL,An_row_ind_TOTAL,lambda_list,...
        eps_init,eps_decay,verbose,max_out_iter,max_in_iter,100,1.5,0.000001,permutation_matrix,0,1);
    results(:,t_id)=[history(1,end);history(5,end)*100/(parts*D)];
    fprintf('obj:%.5f, Perr_ratio:%.2f, seed:%d\n------\n------\n', results(1,t_id), results(2,t_id), seedd)
end
Perr_ratios_for_MRUC = results(2,:);
mean_for_MRUC = mean(Perr_ratios_for_MRUC);
min_for_MRUC  = min(Perr_ratios_for_MRUC);
max_for_MRUC  = max(Perr_ratios_for_MRUC);
std_for_MRUC  = std(Perr_ratios_for_MRUC);
fprintf('MRUC: mean Perr = %.4f, std Perr = %.4f, min Perr = %.4f, max Perr = %.4f\n--------\n--------\n',mean_for_MRUC,std_for_MRUC,min_for_MRUC,max_for_MRUC)
end

end
%% utils
if use_diary
    diary off
end