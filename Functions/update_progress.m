function update_progress(iteration, total_iterations)
    % ���½�����
    progress_bar_length = 50; % ����������
    percent_complete = iteration / total_iterations;
    progress_chars = round(percent_complete * progress_bar_length);
    
%     display_progress();
    % ����������������ʼλ��
    fprintf(repmat('\b', 1, progress_bar_length + 7));
    
    % ��ӡ���µĽ�����
    fprintf('[');
    fprintf(repmat('#', 1, progress_chars));
    fprintf(repmat(' ', 1, progress_bar_length - progress_chars));
    fprintf('] %3.0f%%', percent_complete * 100);

    % �����һ��ʱ����
    if iteration == total_iterations
        fprintf('\n');
    end
end

% function display_progress()
%     % ��ʼ������
%     progress_bar_length = 50; % ����������
% 
%     % ��ӡ��ʼ�Ľ�����
%     fprintf('Progress: [');
%     fprintf(repmat(' ', 1, progress_bar_length));
%     fprintf('] 0%%\n');  % ��ĩβ��ӻ��з��Ա��⸲��
% end