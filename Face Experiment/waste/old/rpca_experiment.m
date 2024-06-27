close all
clear
addpath(genpath('C:\Users\Pastoral\Desktop\Copy_SSC_ADMM_v1.1'));
addpath(genpath('C:\Users\Pastoral\Desktop\Unlabeled_PCA_Code'));
addpath(genpath('C:\Users\Pastoral\Desktop\Pursuit\Functions'));

load('YaleB_cell');
num_groups = 4;

h = 192; w = 168; 
hh = 48; ww = 42;
[~,n] = size(YaleB_cell{1,1});
perm_flag = 1;
outliers_num = 2; % number of outliers out of 64 points
v_patches = 48;
h_patches = 42;


X = zeros(hh*ww,n*num_groups);

for i = 1:num_groups
    X(:,(i-1)*n+1:i*n) = DSP(YaleB_cell{1,i},h,w,hh,ww);
        shuffled_ids = randperm(n);
    outliers_id = shuffled_ids(1:outliers_num);
    if length(setdiff(outliers_id, 46))== outliers_num
        outliers_id(1) = 46;
    end
    if length(setdiff(outliers_id, 23))== outliers_num
        outliers_id(2) = 23;
    end

    [Y] = permute_face(X(:,(i-1)*n+1:i*n), outliers_id, hh, ww, v_patches, h_patches, perm_flag);
    X(:,(i-1)*n+1:i*n) = Y;
    
end

lambda = 0.03;
[A,E] = rpca(X,lambda);
f_id = 0;

for i = 1:num_groups
for j = 45

img_title = 'raw';
image_face(X(:, (i-1)*n+j), f_id, hh, ww, img_title);
f_id = f_id + 1;


img_title = 'rpca';
image_face(A(:, (i-1)*n+j), f_id, hh, ww, img_title);
f_id = f_id + 1;

img_title = 'noise';
image_face(E(:, (i-1)*n+j), f_id, hh, ww, img_title);
f_id = f_id + 1;

end
end
