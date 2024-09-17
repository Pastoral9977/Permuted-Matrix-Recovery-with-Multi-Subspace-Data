function update_progress(iteration, total_iterations)
    % 更新进度条
    progress_bar_length = 50; % 进度条长度
    percent_complete = iteration / total_iterations;
    progress_chars = round(percent_complete * progress_bar_length);
    
%     display_progress();
    % 回退至进度条的起始位置
    fprintf(repmat('\b', 1, progress_bar_length + 7));
    
    % 打印更新的进度条
    fprintf('[');
    fprintf(repmat('#', 1, progress_chars));
    fprintf(repmat(' ', 1, progress_bar_length - progress_chars));
    fprintf('] %3.0f%%', percent_complete * 100);

    % 在最后一轮时换行
    if iteration == total_iterations
        fprintf('\n');
    end
end

% function display_progress()
%     % 初始化变量
%     progress_bar_length = 50; % 进度条长度
% 
%     % 打印初始的进度条
%     fprintf('Progress: [');
%     fprintf(repmat(' ', 1, progress_bar_length));
%     fprintf('] 0%%\n');  % 在末尾添加换行符以避免覆盖
% end