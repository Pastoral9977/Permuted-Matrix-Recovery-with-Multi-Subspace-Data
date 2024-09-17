clear
close all

% Load data from a MATLAB file
load('RESULT_single_subspace.mat')

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

%% Plot and save figure 1: Recovery Error For Single Subspace
figure('Position', [left+offset, bottom-offset*0.5, width, height]);
imagesc(shuffled_ratio_set, rank_ratio_set, abs(recover_error_singlegroup));
colormap(flipud(gray));
colorbar;
caxis(colorbar_limits); 
set(gca, 'YDir', 'normal', 'FontSize', 17, 'FontName', 'Times New Roman');
set(gca, 'XTick', shuffled_ratio_set, 'XTickLabel', cellfun(@(x) strrep(num2str(x, '%.1f'), '0.', '.'), num2cell(shuffled_ratio_set), 'UniformOutput', false));
set(gca, 'YTick', 1:6:length(rank_ratio_set), 'YTickLabel', rank_ratio_set(1:6:length(rank_ratio_set)));
h = ylabel('r');
set(h, 'Rotation', pi/2);
xlabel('Shuffled ratio');
title('$\mathit{\mathbf{SingleRE_{gt}}}$', 'Interpreter', 'latex')
saveas(gcf, 'Recovery_Error_For_Single_Subspace.png')



%% Function to convert tables to matrices and compute median
function max_matrix = tables2maxmatrices(tables_cell)
    matrices_cell = cell(size(tables_cell));
    for i = 1:numel(tables_cell)
        matrices_cell{i} = table2array(tables_cell{i});
    end
    max_matrix = median(cat(3, matrices_cell{:}), 3);
end
