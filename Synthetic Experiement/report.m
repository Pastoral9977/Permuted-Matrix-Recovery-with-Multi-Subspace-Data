% ���ر���Ľ������
load('RESULT_Group.mat', 'RESULT');

% ��ȡ��ͬ��num_groupsȡֵ
num_groups_set = [2, 3, 5, 8, 10];

% ��ʼ���㱨����
report = '';

% ����ÿ��num_groupsȡֵ�������ɻ㱨
for g_idx = 1:length(num_groups_set)
    num_groups = num_groups_set(g_idx);
    
    % ��ȡ��Ӧ�Ľ����
    missrate_out_table = RESULT.missrate_out_allgroup{g_idx};
    missrate_out_gtbases_table = RESULT.missrate_out_gtbases_allgroup{g_idx};
    recover_error_table = RESULT.recover_error_allgroup{g_idx};
    recover_error_gtbases_table = RESULT.recover_error_gtbases_allgroup{g_idx};
    missrate_in_table = RESULT.missrate_in_allgroup{g_idx};
    undetected_outlier_ratio_table = RESULT.undetected_outlier_ratio_allgroup{g_idx};
    
    % ���ɻ㱨����
    report = [report, sprintf('========================================\n')];
    report = [report, sprintf(' Experiement Settings: num_groups = %d\n', num_groups)];
    report = [report, sprintf('========================================\n')];
    
    % ��ӱ������
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

% ����㱨���ı��ļ�
fid = fopen('Report.txt', 'wt');
fprintf(fid, report);
fclose(fid);

% ��ʾ�㱨����
fprintf(report);

function str = tableToString(T)
    % ����ת��Ϊ�ַ������Ҷ���
    str = '';
    % ��ȡ����
    colnames = T.Properties.VariableNames;
    % ��ȡ����
    rownames = T.Properties.RowNames;
    % ��ȡ����
    data = table2cell(T);
    
    % ��ȡ�п���ʼΪ�����ĳ���
    col_width = cellfun(@length, colnames, 'UniformOutput', true);
    row_width = max(cellfun(@length, rownames));
    
    % �����п�ȷ����������������
    for i = 1:width(T)
        col_width(i) = max([col_width(i), max(cellfun(@(x) length(sprintf('%.4f', x)), data(:, i)))]);
    end
    
    % ��ӱ�ͷ
    u = 1;
    v = 5;
    header = [' ' * ones(1, row_width + v)];
    for i = 1:length(colnames)
        header = [header, sprintf(['%-', num2str(col_width(i) + u), 's'], colnames{i})];
    end
    str = [str, '', header, '\n'];
    
    % �������
    for i = 1:height(T)
        row = sprintf(['%-', num2str(row_width + v), 's'], rownames{i});
        for j = 1:width(T)
            row = [row, sprintf(['%-', num2str(col_width(j) + u), '.4f'], data{i, j})];
        end
        str = [str, row, '\n'];
    end
end

