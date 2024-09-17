%% figure
save_subfigs = false;
g = figure;
g.Units = 'centimeters';
shrink = 0.9;
g.Position = [4 3.1 3.5*num_groups*shrink 13*shrink];
rows = 3;
for ii = 1:num_groups
    subplot(rows, num_groups, ii)
    img = imgs{ii};
    img = img/max(img, [], 'all');
    imshow(img)
    title('Original')

    subplot(rows, num_groups, num_groups+ii)
    perm_img = perm_imgs{ii};
    perm_img = perm_img/max(perm_img, [], 'all');
    imshow(perm_img)
    title('Corrupted')
    
    if reco
        subplot(rows, num_groups, 2*num_groups+ii)
        reco_img = reco_imgs{ii};
        reco_img = reco_img/max(reco_img, [], 'all');
        imshow(reco_img)
        title('MRUC')
    end
end

if save_subfigs
    restoredefaultpath
    base_folder = 'Experiment_Images';
    if ~exist(base_folder, 'dir')
        mkdir(base_folder);
    end

    for ii  = 1:num_groups
        id = sample_id+(ii-1)*64;
        img = imgs{ii};
        perm_img = perm_imgs{ii};
        reco_img = reco_imgs{ii};
        img = img/max(img, [], 'all');
        perm_img = perm_img/max(perm_img, [], 'all');
        reco_img = reco_img/max(reco_img, [], 'all');

        [n, m] = size(img);
        save_image(img, ii, sample_id, n, m, 'Original_NoDSP', base_folder);
        save_image(perm_img, ii, sample_id, n, m, 'Corrupted_NoDSP', base_folder);
        save_image(reco_img, ii, sample_id, n, m, 'MRUC', base_folder);
    end
end