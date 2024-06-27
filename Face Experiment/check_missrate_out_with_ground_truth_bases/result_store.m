% 数据导入
data = [
    2, 1, 0.000, 0.000, 0.031, 0.010526;
    2, 6, 0.000, 0.000, 0.094, 0.030303;
    2, 11, 0.000, 0.000, 0.156, 0;
    2, 16, 0.000, 0.000, 0.125, 0.041667;
    2, 21, 0.000, 0.000, 0.125, 0;
    2, 26, 0.000, 0.000, 0.094, 0;
    3, 1, 0.000, 0.021, 0.021, 0.027972;
    3, 6, 0.000, 0.022, 0.063, 0.020833;
    3, 11, 0.000, 0.023, 0.083, 0.0074074;
    3, 16, 0.000, 0.000, 0.063, 0;
    3, 21, 0.000, 0.000, 0.021, 0.0070922;
    3, 26, 0.000, 0.000, 0.042, 0;
    5, 1, 0.000, 0.026, 0.025, 0.039301;
    5, 6, 0.013, 0.064, 0.025, 0.16034;
    5, 11, 0.025, 0.051, 0.025, 0.049327;
    5, 16, 0.025, 0.013, 0.050, 0.016949;
    5, 21, 0.025, 0.039, 0.050, 0.021552;
    5, 26, 0.013, 0.013, 0.025, 0.0086207;
    8, 1, 0.023, 0.079, 0.008, 0.12121;
    8, 6, 0.016, 0.097, 0.031, 0.23055;
    8, 11, 0.039, 0.063, 0.000, 0.10086;
    8, 16, 0.023, 0.031, 0.008, 0.013889;
    8, 21, 0.023, 0.024, 0.031, 0.0083799;
    8, 26, 0.023, 0.063, 0.016, 0.10056;
    10, 1, 0.019, 0.088, 0.006, 0.12896;
    10, 6, 0.013, 0.210, 0.019, 0.26128;
    10, 11, 0.031, 0.081, 0.000, 0.18357;
    10, 16, 0.025, 0.044, 0.013, 0.015695;
    10, 21, 0.037, 0.032, 0.013, 0.013793;
    10, 26, 0.025, 0.146, 0.013, 0.19818;
    12, 1, 0.016, 0.089, 0.005, 0.17761;
    12, 6, 0.021, 0.117, 0.021, 0.2008;
    12, 11, 0.036, 0.130, 0.000, 0.21667;
    12, 16, 0.026, 0.042, 0.010, 0.017375;
    12, 21, 0.031, 0.042, 0.010, 0.017682;
    12, 26, 0.031, 0.195, 0.010, 0.21538;
];

% 转换为表格
T = array2table(data, 'VariableNames', {'num_groups', 'start_person_id', 'missrate_out_gt', 'missrate_out_reconstructed', 'undetected_outlier_num_ratio', 'missrate_in'});

% 生成透视矩阵
missrate_out_gt_pivot = unstack(T(:, {'num_groups', 'start_person_id', 'missrate_out_gt'}), 'missrate_out_gt', 'start_person_id');
missrate_out_reconstructed_pivot = unstack(T(:, {'num_groups', 'start_person_id', 'missrate_out_reconstructed'}), 'missrate_out_reconstructed', 'start_person_id');
undetected_outlier_num_ratio_pivot = unstack(T(:, {'num_groups', 'start_person_id', 'undetected_outlier_num_ratio'}), 'undetected_outlier_num_ratio', 'start_person_id');
missrate_in_pivot = unstack(T(:, {'num_groups', 'start_person_id', 'missrate_in'}), 'missrate_in', 'start_person_id');

% 提取矩阵并计算均值和中位数
missrate_out_gt_matrix = table2array(missrate_out_gt_pivot(:, 2:end));
missrate_out_reconstructed_matrix = table2array(missrate_out_reconstructed_pivot(:, 2:end));
undetected_outlier_num_ratio_matrix = table2array(undetected_outlier_num_ratio_pivot(:, 2:end));
missrate_in_matrix = table2array(missrate_in_pivot(:, 2:end));

% 计算均值和中位数
missrate_out_gt_stats = [mean(missrate_out_gt_matrix, 2), median(missrate_out_gt_matrix, 2)];
missrate_out_reconstructed_stats = [mean(missrate_out_reconstructed_matrix, 2), median(missrate_out_reconstructed_matrix, 2)];
undetected_outlier_num_ratio_stats = [mean(undetected_outlier_num_ratio_matrix, 2), median(undetected_outlier_num_ratio_matrix, 2)];
missrate_in_stats = [mean(missrate_in_matrix, 2), median(missrate_in_matrix, 2)];

% 显示结果
disp('Missrate Out GT Matrix:');
disp(missrate_out_gt_matrix);

disp('Missrate Out Reconstructed Matrix:');
disp(missrate_out_reconstructed_matrix);

disp('Undetected Outlier Num Ratio Matrix:');
disp(undetected_outlier_num_ratio_matrix);

disp('missrate In Matrix:');
disp(missrate_in_matrix);

disp('Missrate Out GT Stats (Mean, Median):');
disp(missrate_out_gt_stats);

disp('Missrate Out Reconstructed Stats (Mean, Median):');
disp(missrate_out_reconstructed_stats);

disp('Undetected Outlier Num Ratio Stats (Mean, Median):');
disp(undetected_outlier_num_ratio_stats);

disp('Missrate In (Mean, Median):');
disp(missrate_in_stats);
