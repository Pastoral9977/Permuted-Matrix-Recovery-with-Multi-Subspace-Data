base_folder = 'Experiment_Images';
if ~exist(base_folder, 'dir')
    mkdir(base_folder);
end

for j = [54]
    f_id = 0;
    for i = 1:num_groups

        img_title = 'Ground Truth';
        folder_path = fullfile(base_folder, img_title);
        if ~exist(folder_path, 'dir')
            mkdir(folder_path);
        end
        save_image(X_gt(:, sum(nn(1:i))+j), f_id, j, hh, ww, img_title, folder_path);
        f_id = f_id + 1;

        img_title = 'Outlier';
        folder_path = fullfile(base_folder, img_title);
        if ~exist(folder_path, 'dir')
            mkdir(folder_path);
        end
        save_image(X_tilde(:, sum(nn(1:i))+j), f_id, j, hh, ww, img_title, folder_path);
        f_id = f_id + 1;

        img_title = 'Recovered';
        folder_path = fullfile(base_folder, img_title);
        if ~exist(folder_path, 'dir')
            mkdir(folder_path);
        end
        save_image(X_solved(:, sum(nn(1:i))+j), f_id, j, hh, ww, img_title, folder_path);
        f_id = f_id + 1;

    end
end

function save_image(xj, f_id, img_id, h, w, img_title, folder_path)
    Imagej = reshape(xj, [h, w]);
    f = figure;
    image(Imagej); 
    colormap(gray);
    f.Units = 'centimeters';
    f.Position = [1+3.6*fix(f_id/3)+10*mod(f_id,3) 15-4*fix(f_id/3) 3.6 4.8];
    title(img_title);
    set(gca, 'FontSize', 12, 'Fontname','times new Roman'); 
    xticks([]);
    yticks([]); 
    yticklabels([]);
    xticklabels([]);

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

    % Save the figure
    saveas(f, fullfile(folder_path, [img_title, '_', num2str(img_id), '_', num2str(next_seq), '.png']));
    close(f);
end

