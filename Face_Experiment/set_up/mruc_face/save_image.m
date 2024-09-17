function save_image(xj, f_id, img_id, h, w, img_title, base_folder)
    % Reshape and display the image
    Imagej = reshape(xj, [h, w]);
    f = figure;
    imshow(Imagej);
    f.Units = 'centimeters';
    f.Position = [1+3.5*fix(f_id/3)+9*mod(f_id,3) 15-4*fix(f_id/3) 3.6 4.8];
    parts = strsplit(img_title, '_');
    title(parts(1));
    set(gca, 'Units', 'normalized', 'Position', [0.125 0.125 0.75 0.75]);  % ½« axes ÌîÂúÕû¸ö figure
    set(gca, 'FontSize', 12, 'Fontname', 'Times New Roman');
    xticks([]);
    yticks([]);
    yticklabels([]);
    xticklabels([]);
    
    % Create the folder path if it does not exist
    folder_path = fullfile(base_folder, img_title);
    if ~exist(folder_path, 'dir')
        mkdir(folder_path);
    end
    
    % Determine the next available sequence number
    file_pattern = fullfile(folder_path, [img_title, '_', num2str(img_id), '_*.png']);
    files = dir(file_pattern);
    max_seq = 0;
    for k = 1:length(files)
        [~, name, ~] = fileparts(files(k).name);
        tokens = regexp(name, ['_', num2str(img_id), '_(\d+)$'], 'tokens');
        if ~isempty(tokens)
            seq_num = str2double(tokens{1}{1});
            if seq_num > max_seq
                max_seq = seq_num;
            end
        end
    end
    next_seq = max_seq + 1;
    
    % Construct the full file path
    file_name = [img_title, '_', num2str(img_id), '_', num2str(next_seq), '.png'];
    invalid_chars = ['\', '/', ':', '*', '?', '"', '<', '>', '|'];
    for i = 1:length(invalid_chars)
        file_name = strrep(file_name, invalid_chars(i), '_');
    end
    full_file_path = fullfile(folder_path, file_name);
    
    % Save the figure
    try
        saveas(f, full_file_path);
    catch ME
        warning('Failed to save the image: %s', ME.message);
    end
    close(f);
end