clear
close all
% 
load('RESULT_Others')
addpath(genpath(pwd))
% Define ranges and settings for plots
folderPath = 'Experiment_Images/Synthetic';
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end
rank_ratio_set = 1:1:25;
shuffled_ratio_set = 0.1:0.1:0.6;
colorbar_limits = [0, 0.57];

left = 160;
bottom = 500;
width = 260; 
height = 230;
offset = 280;

%% Plot and save figure 1: RKPCA
figure('Position', [left, bottom, width, height]);
imagesc(shuffled_ratio_set, rank_ratio_set, abs(recovery_error_RKPCA_mat));
colormap(flipud(gray));
colorbar;
caxis(colorbar_limits); 
set(gca, 'YDir', 'normal', 'FontSize', 17, 'FontName', 'Times New Roman');
set(gca, 'XTick', shuffled_ratio_set, 'XTickLabel', cellfun(@(x) strrep(num2str(x, '%.1f'), '0.', '.'), num2cell(shuffled_ratio_set), 'UniformOutput', false));
set(gca, 'YTick', 1:6:length(rank_ratio_set), 'YTickLabel', rank_ratio_set(1:6:length(rank_ratio_set)));
h = ylabel('r');
set(h, 'Rotation', pi/2);
xlabel('Shuffled Ratio');
title('RKPCA', 'Interpreter', 'latex')
saveas(gcf, 'Experiment_Images/Synthetic/RKPCA_synthetic.png')

%% Plot and save figure 2: RKPCA PMSDR
figure('Position', [left+offset, bottom, width, height]);
imagesc(shuffled_ratio_set, rank_ratio_set, abs(recovery_error_PMSDR_RKPCA_mat));
colormap(flipud(gray));
colorbar;
caxis(colorbar_limits);  
set(gca, 'YDir', 'normal', 'FontSize', 17, 'FontName', 'Times New Roman');
set(gca, 'XTick', shuffled_ratio_set, 'XTickLabel', cellfun(@(x) strrep(num2str(x, '%.1f'), '0.', '.'), num2cell(shuffled_ratio_set), 'UniformOutput', false));
set(gca, 'YTick', 1:6:length(rank_ratio_set), 'YTickLabel', rank_ratio_set(1:6:length(rank_ratio_set)));
h = ylabel('r');
set(h, 'Rotation', pi/2);
xlabel('Shuffled Ratio');
title('RKPCA-S', 'Interpreter', 'latex')
saveas(gcf, 'Experiment_Images/Synthetic/RKPCA_PMSDR_synthetic.png')

%% Plot and save figure 3: RPCA
figure('Position', [left+offset*2, bottom, width, height]);
imagesc(shuffled_ratio_set, rank_ratio_set, abs(recovery_error_RPCA_mat));
colormap(flipud(gray));
colorbar;
caxis(colorbar_limits); 
set(gca, 'YDir', 'normal', 'FontSize', 17, 'FontName', 'Times New Roman');
set(gca, 'XTick', shuffled_ratio_set, 'XTickLabel', cellfun(@(x) strrep(num2str(x, '%.1f'), '0.', '.'), num2cell(shuffled_ratio_set), 'UniformOutput', false));
set(gca, 'YTick', 1:6:length(rank_ratio_set), 'YTickLabel', rank_ratio_set(1:6:length(rank_ratio_set)));
h = ylabel('r');
set(h, 'Rotation', pi/2);
xlabel('Shuffled Ratio');
title('RPCA', 'Interpreter', 'latex')
saveas(gcf, 'Experiment_Images/Synthetic/RPCA_synthetic.png')

%% Plot and save figure 4: PMSDR RPCA
figure('Position', [left+offset*3, bottom, width, height]);
imagesc(shuffled_ratio_set, rank_ratio_set, abs(recovery_error_PMSDR_RPCA_mat));
colormap(flipud(gray));
colorbar;
caxis(colorbar_limits);  
set(gca, 'YDir', 'normal', 'FontSize', 17, 'FontName', 'Times New Roman');
set(gca, 'XTick', shuffled_ratio_set, 'XTickLabel', cellfun(@(x) strrep(num2str(x, '%.1f'), '0.', '.'), num2cell(shuffled_ratio_set), 'UniformOutput', false));
set(gca, 'YTick', 1:6:length(rank_ratio_set), 'YTickLabel', rank_ratio_set(1:6:length(rank_ratio_set)));
h = ylabel('r');
set(h, 'Rotation', pi/2);
xlabel('Shuffled Ratio');
title('RPCA-S', 'Interpreter', 'latex')
saveas(gcf, 'Experiment_Images/Synthetic/RPCA_PMSDR_synthetic.png')

%% Plot and save figure 5: Subspace Clustering Error Missrate In
figure('Position', [left, bottom-offset*1.15, width, height]);
imagesc(shuffled_ratio_set, rank_ratio_set, abs(recovery_error_SSC_mat));
colormap(flipud(gray));
colorbar;
caxis(colorbar_limits);  
set(gca, 'YDir', 'normal', 'FontSize', 17, 'FontName', 'Times New Roman');
set(gca, 'XTick', shuffled_ratio_set, 'XTickLabel', cellfun(@(x) strrep(num2str(x, '%.1f'), '0.', '.'), num2cell(shuffled_ratio_set), 'UniformOutput', false));
set(gca, 'YTick', 1:6:length(rank_ratio_set), 'YTickLabel', rank_ratio_set(1:6:length(rank_ratio_set)));
h = ylabel('r');
set(h, 'Rotation', pi/2);
xlabel('Shuffled Ratio');
title('SSC', 'Interpreter', 'latex')
saveas(gcf, 'Experiment_Images/Synthetic/SSC_synthetic.png')

%% Plot and save figure 6: Undetected Outlier Ratio
figure('Position', [left+offset, bottom-offset*1.15, width, height]);
imagesc(shuffled_ratio_set, rank_ratio_set, abs(recovery_error_PMSDR_SSC_mat));
colormap(flipud(gray));
colorbar;
caxis(colorbar_limits);  
set(gca, 'YDir', 'normal', 'FontSize', 17, 'FontName', 'Times New Roman');
set(gca, 'XTick', shuffled_ratio_set, 'XTickLabel', cellfun(@(x) strrep(num2str(x, '%.1f'), '0.', '.'), num2cell(shuffled_ratio_set), 'UniformOutput', false));
set(gca, 'YTick', 1:6:length(rank_ratio_set), 'YTickLabel', rank_ratio_set(1:6:length(rank_ratio_set)));
h = ylabel('r');
set(h, 'Rotation', pi/2);
xlabel('Shuffled Ratio');
title('SSC-S', 'Interpreter', 'latex')
saveas(gcf, 'Experiment_Images/Synthetic/SSC_PMSDR_synthetic.png')
