% if ~exist('X_RKPCA', 'var')
%     load('figure_plot_results_rkpca')
% end
if ~exist('use_PMSDR', 'var')
    use_PMSDR = true;
end
save_subfigs = false;

%% Plot figures
close all
restoredefaultpath
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

sub_id = 23;
for ii  = 1:num_groups
    id = sub_id+(ii-1)*64;
    img = reshape(X_gt(:, id),n,m);
    perm_img = reshape(X_tilde(:, id),n,m);
    reco_img = reshape(X_RKPCA(:,id),n,m);
    img = img/max(img, [], 'all');
    perm_img = perm_img/max(perm_img, [], 'all');
    reco_img = reco_img/max(reco_img, [], 'all');

    subplot(rows,num_groups,ii)
    imshow(img)
    title('Original')

    subplot(rows,num_groups,ii+num_groups)
    imshow(perm_img)
    title('Corrupted')
    colormap(gray);

    subplot(rows,num_groups,ii+2*num_groups)
    imshow(reco_img)
    title('RKPCA')
    colormap(gray);
    
    if use_PMSDR
        reco_img_PMSDR = reshape(X_PMSDR(:,id),n,m);
        reco_img_PMSDR_RKPCA = reshape(X_PMSDR_RKPCA(:,id),n,m);
        reco_img_PMSDR = reco_img_PMSDR/max(reco_img_PMSDR, [], 'all');
        reco_img_PMSDR_RKPCA = reco_img_PMSDR_RKPCA/max(reco_img_PMSDR_RKPCA, [], 'all');
        
        subplot(rows,num_groups,ii+3*num_groups)
        imshow(reco_img_PMSDR_RKPCA)
        title('PMSDR RKPCA')
        
        subplot(rows,num_groups,ii+4*num_groups)
        imshow(reco_img_PMSDR)
        title('PMSDR')
    end
end
sgtitle(['$\lambda_{whole} = ' num2str(lambda) ', \lambda_{sub} = ' num2str(lambda) ', c = ' num2str(c) '$'], 'interpreter', 'latex')


%% Save figures
base_folder = 'Experiment_Images';
if ~exist(base_folder, 'dir')
    mkdir(base_folder);
end
saveas(g, fullfile(base_folder, 'A.png'));

for ii  = 1:num_groups
    id = sub_id+(ii-1)*64;
    img = reshape(X_gt(:, id),n,m);
    perm_img = reshape(X_tilde(:, id),n,m);
    reco_img = reshape(X_RKPCA(:,id),n,m);
    img = img/max(img, [], 'all');
    perm_img = perm_img/max(perm_img, [], 'all');
    reco_img = reco_img/max(reco_img, [], 'all');

    save_image(img, ii, sub_id, 'Original', base_folder);
    save_image(perm_img, ii, sub_id, 'Corrupted', base_folder);
    save_image(reco_img, ii, sub_id,'RKPCA', base_folder);
    
    if use_PMSDR
        reco_img_PMSDR_RKPCA = reshape(X_PMSDR_RKPCA(:,id),n,m);
        reco_img_PMSDR = reshape(X_PMSDR(:,id),n,m);
        reco_img_PMSDR_RKPCA = reco_img_PMSDR_RKPCA/max(reco_img_PMSDR_RKPCA, [], 'all');
        reco_img_PMSDR = reco_img_PMSDR/max(reco_img_PMSDR, [], 'all');

        save_image(reco_img_PMSDR_RKPCA, ii, sub_id, 'RKPCA-S', base_folder);
        save_image(reco_img_PMSDR, ii, sub_id, 'PMSDR', base_folder);
    end
end

save('figure_plot_results_rkpca', 'g', 'X_gt', 'X_tilde', 'X_RKPCA', 'num_groups', 'n', 'm', 'use_PMSDR','lambda','c');
if use_PMSDR
    save('figure_plot_results_rkpca', 'g', 'X_gt', 'X_tilde', 'X_RKPCA', 'X_PMSDR_RKPCA','X_PMSDR', 'num_groups', 'n', 'm', 'use_PMSDR','lambda','c');
end

function save_image(img, f_id, img_id, img_title, base_folder)
    f = figure;
    imshow(img)
    f.Units = 'centimeters';
    f.Position = [1+3.6*fix(f_id/3)+10*mod(f_id,3) 15-4*fix(f_id/3) 3.6 4.8];
    title(img_title);
    set(gca, 'FontSize', 12, 'Fontname','times new Roman'); 
    set(gca, 'Units', 'normalized', 'Position', [0.125 0.125 0.75 0.75]); 
    xticks([]);
    yticks([]); 
    yticklabels([]);
    xticklabels([]);
    
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

    % Save the figure
    saveas(f, fullfile(folder_path, [img_title, '_', num2str(img_id), '_', num2str(next_seq), '.png']));
    close(f);
end