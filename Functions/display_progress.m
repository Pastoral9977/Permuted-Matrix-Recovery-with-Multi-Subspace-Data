function display_progress()
    % ��ʼ������
    progress_bar_length = 50; % ����������

    % ��ӡ��ʼ�Ľ�����
    fprintf('\tProgress: [');
    fprintf(repmat(' ', 1, progress_bar_length));
    fprintf('] 0%%\n');  % ��ĩβ��ӻ��з��Ա��⸲��
end