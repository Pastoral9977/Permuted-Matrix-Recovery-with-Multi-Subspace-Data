close all
clear
addpath(genpath('../Functions'));
addpath(genpath(pwd));
%% load and corrupt data
use_downsampled_version = true;
if ~use_downsampled_version
    suffix = '';
    n = 192;
    m = 168;
    patch = 6;
else
    suffix = '_light';
    n = 48;
    m = 42;
    np = 8;
    mp = 7;
    patch = 1;
    missing_ratio = 0.0;
    shuffled_ratio = 0.4; 
    seed = 6218348;
end
rng(seed)
reco = true;
fprintf('seed = %d\n', seed)
sample_id = 23;
outliers_num = 16;
ids = 1:10; 
num_groups = length(ids);
imgs = cell(1, num_groups);
perm_imgs = cell(1, num_groups);
perm_indices = cell(1, num_groups);
for ii = 1:num_groups
    id = ids(ii);
    load('YaleB_cell.mat')
    images = YaleB_cell{id};
    images = DownsamplePicture(images, 192, 168, n, m);
    imgs{ii} = reshape(images(:, sample_id),n,m);
    outliers_id = randperm(size(images,2), outliers_num);
    specific_ids = [46, 23, 54, 12];
    for j = 1:length(specific_ids)
        if length(setdiff(outliers_id, specific_ids(j))) == outliers_num
            outliers_id(j) = specific_ids(j);
        end
    end
    [perm_images, perm_indices{ii}] = permute_face_randomly_MRUC_for_perm_indices(images, outliers_id, n, m, np, mp, sample_id, missing_ratio, shuffled_ratio, seed*ii);
    perm_imgs{ii} = reshape(perm_images(:, sample_id), n, m);
end

%% recover by MRUC
% algorithm set up
max_in_iter=10000;
eps_init=100;
max_out_iter=200000;
verbose=false;
use_rank=4;
lr1=0.1;
lr2=1;

init_imgs = cell(1, num_groups);
reco_imgs = cell(1, num_groups);
sub_observes = cell(1, num_groups);
for ii = 1:num_groups
    fprintf('initialization progress: %d/%d\n',ii,num_groups)
    [init_imgs{ii}, sub_observes{ii}] = initialize_image(perm_imgs{ii}, np, mp);
end
fprintf('\n')

if reco
    for ii = 1:num_groups
        tic
        fprintf('recovery progress: %d/%d',ii,num_groups)
        reco_imgs{ii} = MRUC_FaceRecovery(max_in_iter,eps_init,max_out_iter,verbose,use_rank,init_imgs{ii},sub_observes{ii},lr1,lr2,np,mp,perm_indices{ii});
        fprintf('...%.2fs\n',toc)
    end
end

run Figure_Face_MRUC





