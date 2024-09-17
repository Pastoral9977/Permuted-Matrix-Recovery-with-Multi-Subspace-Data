if ~exist('X_RPCA', 'var')
    load('figure_plot_results')
end
if ~exist('use_PMSDR', 'var')
    use_PMSDR = true;
end

%% Plot figures
close all
show_full_image = true;
if show_full_image
    g = figure;
    g.Units = 'centimeters';
    shrink = 0.9;
    if ~use_PMSDR
        g.Position = [4 1.1 3.5*num_groups*shrink 13*shrink];
        rows = 3;
    else
        g.Position = [8 1.1 3.5*num_groups*shrink 20*shrink];
        rows = 5;
    end
    sample_id = 23;
    for ii  = 1:num_groups
        id = sample_id+(ii-1)*64;
        img = reshape(X_gt(:, id),n,m);
        perm_img = reshape(X_tilde(:, id),n,m);
        reco_img = reshape(X_RPCA(:,id),n,m);
        img = img/max(img, [], 'all');
        perm_img = perm_img/max(perm_img, [], 'all');
        reco_img = reco_img/max(reco_img, [], 'all');

        subplot(rows,num_groups,ii)
        imshow(img)
        title('Original')

        subplot(rows,num_groups,ii+num_groups)
        imshow(perm_img)
        title('Corrupted')

        subplot(rows,num_groups,ii+2*num_groups)
        imshow(reco_img)
        title('RPCA')

        if use_PMSDR
            reco_img_PMSDR = reshape(X_PMSDR(:,id),n,m);
            reco_img_PMSDR_RPCA = reshape(X_PMSDR_RPCA(:,id),n,m);
            reco_img_PMSDR = reco_img_PMSDR/max(reco_img_PMSDR, [], 'all');
            reco_img_PMSDR_RPCA = reco_img_PMSDR_RPCA/max(reco_img_PMSDR_RPCA, [], 'all');

            subplot(rows,num_groups,ii+3*num_groups)
            imshow(reco_img_PMSDR_RPCA)
            title('RPCA-S')

            subplot(rows,num_groups,ii+4*num_groups)
            imshow(reco_img_PMSDR)
            title('PMSDR')
        end
    end
end

restoredefaultpath
base_folder = 'Experiment_Images';
if ~exist(base_folder, 'dir')
    mkdir(base_folder);
end

for sample_id = [23]
    f_id = 0;
    for ii = 1:num_groups
        
        img = reshape(X_gt(:, sample_id+(ii-1)*64),n, m);
        perm_img = reshape(X_tilde(:, sample_id+(ii-1)*64),n, m);
        reco_img = reshape(X_RPCA(:,sample_id+(ii-1)*64),n, m);
        img = img/max(img, [], 'all');
        perm_img = perm_img/max(perm_img, [], 'all');
        reco_img = reco_img/max(reco_img, [], 'all');
        
        img_title = 'Original';
        folder_path = fullfile(base_folder, img_title);
        if ~exist(folder_path, 'dir')
            mkdir(folder_path);
        end
        save_image(img, f_id, sample_id, img_title, folder_path);
        f_id = f_id + 1;

        img_title = 'Corrupted';
        folder_path = fullfile(base_folder, img_title);
        if ~exist(folder_path, 'dir')
            mkdir(folder_path);
        end
        save_image(perm_img, f_id, sample_id, img_title, folder_path);
        f_id = f_id + 1;

        img_title = 'RPCA';
        folder_path = fullfile(base_folder, img_title);
        if ~exist(folder_path, 'dir')
            mkdir(folder_path);
        end
        save_image(reco_img, f_id, sample_id, img_title, folder_path);
        f_id = f_id + 1;
        
        if use_PMSDR
            reco_img_PMSDR = reshape(X_PMSDR(:,sample_id+(ii-1)*64),n,m);
            reco_img_PMSDR_RPCA = reshape(X_PMSDR_RPCA(:,sample_id+(ii-1)*64),n,m);
            reco_img_PMSDR = reco_img_PMSDR/max(reco_img_PMSDR, [], 'all');
            reco_img_PMSDR_RPCA = reco_img_PMSDR_RPCA/max(reco_img_PMSDR_RPCA, [], 'all');

            img_title = 'RPCA-S';
            folder_path = fullfile(base_folder, img_title);
            if ~exist(folder_path, 'dir')
                mkdir(folder_path);
            end
            save_image(reco_img_PMSDR_RPCA, f_id, sample_id, img_title, folder_path);
            f_id = f_id + 1;

            img_title = 'PMSDR';
            folder_path = fullfile(base_folder, img_title);
            if ~exist(folder_path, 'dir')
                mkdir(folder_path);
            end
            save_image(reco_img_PMSDR, f_id, sample_id, img_title, folder_path);
            f_id = f_id + 1;
        end
    end
end

save('figure_plot_results_rpca', 'g', 'X_gt', 'X_tilde', 'X_RPCA', 'num_groups', 'n', 'm', 'use_PMSDR');
if use_PMSDR
    save('figure_plot_results_rpca', 'g', 'X_gt', 'X_tilde', 'X_RPCA', 'X_PMSDR_RPCA','X_PMSDR', 'num_groups', 'n', 'm', 'use_PMSDR');
end

function save_image(img, f_id, img_id, img_title, folder_path)
    f = figure;
    imshow(img)
    f.Units = 'centimeters';
    f.Position = [1+3.6*fix(f_id/3)+10*mod(f_id,3) 15-4*fix(f_id/3) 3.6 4.8];
    title(img_title);
    set(gca, 'Units', 'normalized', 'Position', [0.125 0.125 0.75 0.75]); 
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