function display_progress()
    % 初始化变量
    progress_bar_length = 50; % 进度条长度

    % 打印初始的进度条
    fprintf('\tProgress: [');
    fprintf(repmat(' ', 1, progress_bar_length));
    fprintf('] 0%%\n');  % 在末尾添加换行符以避免覆盖
end