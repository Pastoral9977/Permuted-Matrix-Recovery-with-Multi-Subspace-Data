% 加载保存的结果数据
load('RESULT_Group.mat', 'RESULT');

% 获取不同的num_groups取值
num_groups_set = [2, 3, 5, 8, 10];

% 初始化汇报内容
report = '';

% 遍历每个num_groups取值，并生成汇报
for g_idx = 1:length(num_groups_set)
    num_groups = num_groups_set(g_idx);
    
    % 获取对应的结果表
    missrate_out_table = RESULT.missrate_out_allgroup{g_idx};
    missrate_out_gtbases_table = RESULT.missrate_out_gtbases_allgroup{g_idx};
    recover_error_table = RESULT.recover_error_allgroup{g_idx};
    recover_error_gtbases_table = RESULT.recover_error_gtbases_allgroup{g_idx};
    missrate_in_table = RESULT.missrate_in_allgroup{g_idx};
    undetected_outlier_ratio_table = RESULT.undetected_outlier_ratio_allgroup{g_idx};
    
    % 生成汇报内容
    report = [report, sprintf('========================================\n')];
    report = [report, sprintf(' Experiement Settings: num_groups = %d\n', num_groups)];
    report = [report, sprintf('========================================\n')];
    
    % 添加表格内容
    report = [report, 'Missrate Out (Reconstructed Bases):\n'];
    report = [report, tableToString(missrate_out_table), '\n'];
    
    report = [report, 'Missrate Out (GT Bases):\n'];
    report = [report, tableToString(missrate_out_gtbases_table), '\n'];
    
    report = [report, 'Recover Error (Reconstructed Bases):\n'];
    report = [report, tableToString(recover_error_table), '\n'];
    
    report = [report, 'Recover Error (Reconstructed Bases with Detected True Outliers):\n'];
    report = [report, tableToString(recover_error_gtbases_table), '\n'];
    
    report = [report, 'Recover Error (GT Bases):\n'];
    report = [report, tableToString(recover_error_gtbases_table), '\n'];
    
    report = [report, 'Missrate In (SSC):\n'];
    report = [report, tableToString(missrate_in_table), '\n'];
    
    report = [report, 'Undetected Outlier Ratio:\n'];
    report = [report, tableToString(undetected_outlier_ratio_table), '\n'];
    
    report = [report, '\n'];
end

% 保存汇报到文本文件
fid = fopen('Report.txt', 'wt');
fprintf(fid, report);
fclose(fid);

% 显示汇报内容
fprintf(report);

function str = tableToString(T)
    % 将表转换为字符串并右对齐
    str = '';
    % 获取列名
    colnames = T.Properties.VariableNames;
    % 获取行名
    rownames = T.Properties.RowNames;
    % 获取数据
    data = table2cell(T);
    
    % 获取列宽，初始为列名的长度
    col_width = cellfun(@length, colnames, 'UniformOutput', true);
    row_width = max(cellfun(@length, rownames));
    
    % 更新列宽，确保能容纳所有数据
    for i = 1:width(T)
        col_width(i) = max([col_width(i), max(cellfun(@(x) length(sprintf('%.4f', x)), data(:, i)))]);
    end
    
    % 添加表头
    u = 1;
    v = 5;
    header = [' ' * ones(1, row_width + v)];
    for i = 1:length(colnames)
        header = [header, sprintf(['%-', num2str(col_width(i) + u), 's'], colnames{i})];
    end
    str = [str, '', header, '\n'];
    
    % 添加数据
    for i = 1:height(T)
        row = sprintf(['%-', num2str(row_width + v), 's'], rownames{i});
        for j = 1:width(T)
            row = [row, sprintf(['%-', num2str(col_width(j) + u), '.4f'], data{i, j})];
        end
        str = [str, row, '\n'];
    end
end

