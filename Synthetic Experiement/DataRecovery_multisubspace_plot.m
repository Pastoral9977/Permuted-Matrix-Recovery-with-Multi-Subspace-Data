clear
close all

% Load data from a MATLAB file
load('RESULT_Group.mat')

% Get field names from loaded data
fields = fieldnames(RESULT);

% Iterate through each field
for i = 1:numel(fields)
    field_name = fields{i};
    tables_cell = RESULT.(field_name);
    
    % Convert tables to matrices and compute median
    max_matrix = tables2maxmatrices(tables_cell);
    
    % Assign matrices to dynamically named variables in the workspace
    assignin('caller', field_name, max_matrix);
end

% Define ranges and settings for plots
rank_ratio_set = 1:1:25;
shuffled_ratio_set = 0.1:0.1:0.6;
colorbar_limits = [0, 0.3];

left = 160;
bottom = 500;
width = 260; 
height = 230;
offset = 280;

%% Plot and save figure 1: Missrate Out with GT Information
figure('Position', [left, bottom, width, height]);
imagesc(shuffled_ratio_set, rank_ratio_set, abs(missrate_out_gtbases_allgroup));
colormap(flipud(gray));
colorbar;
caxis(colorbar_limits); 
set(gca, 'YDir', 'normal', 'FontSize', 17, 'FontName', 'Times New Roman');
set(gca, 'XTick', shuffled_ratio_set, 'XTickLabel', cellfun(@(x) strrep(num2str(x, '%.1f'), '0.', '.'), num2cell(shuffled_ratio_set), 'UniformOutput', false));
set(gca, 'YTick', 1:6:length(rank_ratio_set), 'YTickLabel', rank_ratio_set(1:6:length(rank_ratio_set)));
h = ylabel('r');
set(h, 'Rotation', pi/2);
xlabel('Shuffled Ratio');
title('$\mathit{\mathbf{CE_{gt}}}$', 'Interpreter', 'latex')
saveas(gcf, 'Missrate_Out_with_GT_Information.png')

%% Plot and save figure 2: Missrate Out with Recon Information
figure('Position', [left+offset, bottom, width, height]);
imagesc(shuffled_ratio_set, rank_ratio_set, abs(missrate_out_allgroup));
colormap(flipud(gray));
colorbar;
caxis(colorbar_limits);  
set(gca, 'YDir', 'normal', 'FontSize', 17, 'FontName', 'Times New Roman');
set(gca, 'XTick', shuffled_ratio_set, 'XTickLabel', cellfun(@(x) strrep(num2str(x, '%.1f'), '0.', '.'), num2cell(shuffled_ratio_set), 'UniformOutput', false));
set(gca, 'YTick', 1:6:length(rank_ratio_set), 'YTickLabel', rank_ratio_set(1:6:length(rank_ratio_set)));
h = ylabel('r');
set(h, 'Rotation', pi/2);
xlabel('Shuffled Ratio');
title('$\mathit{\mathbf{CE_{recon}}}$', 'Interpreter', 'latex')
saveas(gcf, 'Missrate_Out_with_Recon_Information.png')

%% Plot and save figure 3: Recovery Error with GT Information
figure('Position', [left+offset*2, bottom, width, height]);
imagesc(shuffled_ratio_set, rank_ratio_set, abs(recover_error_gtbases_allgroup));
colormap(flipud(gray));
colorbar;
caxis(colorbar_limits); 
set(gca, 'YDir', 'normal', 'FontSize', 17, 'FontName', 'Times New Roman');
set(gca, 'XTick', shuffled_ratio_set, 'XTickLabel', cellfun(@(x) strrep(num2str(x, '%.1f'), '0.', '.'), num2cell(shuffled_ratio_set), 'UniformOutput', false));
set(gca, 'YTick', 1:6:length(rank_ratio_set), 'YTickLabel', rank_ratio_set(1:6:length(rank_ratio_set)));
h = ylabel('r');
set(h, 'Rotation', pi/2);
xlabel('Shuffled Ratio');
title('$\mathit{\mathbf{RE_{gt}}}$', 'Interpreter', 'latex')
saveas(gcf, 'Recovery_Error_with_GT_Information.png')

%% Plot and save figure 4: Recovery Error with Recon Information
figure('Position', [left+offset*3, bottom, width, height]);
imagesc(shuffled_ratio_set, rank_ratio_set, abs(recover_error_allgroup));
colormap(flipud(gray));
colorbar;
caxis(colorbar_limits);  
set(gca, 'YDir', 'normal', 'FontSize', 17, 'FontName', 'Times New Roman');
set(gca, 'XTick', shuffled_ratio_set, 'XTickLabel', cellfun(@(x) strrep(num2str(x, '%.1f'), '0.', '.'), num2cell(shuffled_ratio_set), 'UniformOutput', false));
set(gca, 'YTick', 1:6:length(rank_ratio_set), 'YTickLabel', rank_ratio_set(1:6:length(rank_ratio_set)));
h = ylabel('r');
set(h, 'Rotation', pi/2);
xlabel('Shuffled Ratio');
title('$\mathit{\mathbf{RE_{recon}}}$', 'Interpreter', 'latex')
saveas(gcf, 'Recovery_Error_with_Recon_Information.png')

%% Plot and save figure 5: Subspace Clustering Error Missrate In
figure('Position', [left, bottom-offset*1.15, width, height]);
imagesc(shuffled_ratio_set, rank_ratio_set, abs(missrate_in_allgroup));
colormap(flipud(gray));
colorbar;
caxis(colorbar_limits/3);  
set(gca, 'YDir', 'normal', 'FontSize', 17, 'FontName', 'Times New Roman');
set(gca, 'XTick', shuffled_ratio_set, 'XTickLabel', cellfun(@(x) strrep(num2str(x, '%.1f'), '0.', '.'), num2cell(shuffled_ratio_set), 'UniformOutput', false));
set(gca, 'YTick', 1:6:length(rank_ratio_set), 'YTickLabel', rank_ratio_set(1:6:length(rank_ratio_set)));
h = ylabel('r');
set(h, 'Rotation', pi/2);
xlabel('Shuffled Ratio');
title('$\mathit{\mathbf{SCerr}}$', 'Interpreter', 'latex')
saveas(gcf, 'Missrate_In.png')

%% Plot and save figure 6: Undetected Outlier Ratio
figure('Position', [left+offset, bottom-offset*1.15, width, height]);
imagesc(shuffled_ratio_set, rank_ratio_set, abs(undetected_outlier_ratio_allgroup));
colormap(flipud(gray));
colorbar;
caxis(colorbar_limits);  
set(gca, 'YDir', 'normal', 'FontSize', 17, 'FontName', 'Times New Roman');
set(gca, 'XTick', shuffled_ratio_set, 'XTickLabel', cellfun(@(x) strrep(num2str(x, '%.1f'), '0.', '.'), num2cell(shuffled_ratio_set), 'UniformOutput', false));
set(gca, 'YTick', 1:6:length(rank_ratio_set), 'YTickLabel', rank_ratio_set(1:6:length(rank_ratio_set)));
h = ylabel('r');
set(h, 'Rotation', pi/2);
xlabel('Shuffled Ratio');
title('$\mathit{\mathbf{UOratio}}$', 'Interpreter', 'latex')
saveas(gcf, 'Undetected_Outlier_Ratio.png')

%% Function to convert tables to matrices and compute median
function max_matrix = tables2maxmatrices(tables_cell)
    matrices_cell = cell(size(tables_cell));
    for i = 1:numel(tables_cell)
        matrices_cell{i} = table2array(tables_cell{i});
    end
%     max_matrix = median(cat(3, matrices_cell{:}), 3);
%     max_matrix = mean(cat(3, matrices_cell{:}), 3);
    max_matrix = max(cat(3, matrices_cell{:}),[], 3);
%     max_matrix = min(cat(3, matrices_cell{:}),[], 3);
end
