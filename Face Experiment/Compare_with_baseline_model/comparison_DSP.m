close all
clear 
addpath(genpath('..\Copy_SSC_ADMM_v1.1'));
addpath(genpath('..\Functions'));

load('YaleB_cell');
rng(21)

h = 192; w = 168;
hh = 48; ww = 42;

r = 10;

%% Single subspace comparison
num_groups = 1;
f_id = 0;
n = [];
sp = 16;
for i = 1:num_groups
    [~,n(i)] = size(YaleB_cell{1,sp-1+i});
end
N = sum(n); 

outliers_num = 16; % number of outliers out of 64 points
v_patches = 12;
h_patches = 14;
% perm_pattern = 0 for fully shuffling; 1 for partially shuffling.
perm_flag = 0;
detect_flag =2;

X_raw = zeros(hh*ww,N);
X_tilde = zeros(hh*ww,N);
outliers_ID = [];

nn = [0 n];

for i = 1:num_groups
    X_raw(:,sum(nn(1:i))+1:sum(nn(1:i+1))) = DownsamplePiciture(YaleB_cell{1,sp-1+i},h,w,hh,ww);
    
    shuffled_ids = randperm(n(i));
    outliers_id = shuffled_ids(1:outliers_num);
    if length(setdiff(outliers_id, 46))== outliers_num
        outliers_id(1) = 46;
    end
    if length(setdiff(outliers_id, 23))== outliers_num
        outliers_id(2) = 23;
    end
    if length(setdiff(outliers_id, 54))== outliers_num
        outliers_id(3) = 54;
    end
    if length(setdiff(outliers_id, 12))== outliers_num
        outliers_id(4) = 12;
    end
    [X] = permute_face_randomly(X_raw(:,sum(nn(1:i))+1:sum(nn(1:i+1))), outliers_id, hh, ww, v_patches, h_patches, perm_flag);
    X_tilde(:,sum(nn(1:i))+1:sum(nn(1:i+1))) = X;
    outliers_ID = [outliers_ID, sort(outliers_id)+sum(nn(1:i))];
end
inliers_ID = setdiff(1:N,outliers_ID);

% parameters for DPCP
c = n-r;
budget = 10000;
epsilon_J = 1e-9;
mu_min = 1e-15;
maxIter = 1000;
delta = 1e-9;
tic
[B_solved, ~] = DPCP_IRLS_proj(X_tilde,c,delta,maxIter,epsilon_J,budget);
time_upca_basis = toc;

tic
alpha = 1e10;
[Inliers_id, Outliers_id] = Outlier_Detect_v2(X_tilde, detect_flag, alpha);
[Basis_estimated] = sub_esti(X_tilde(:,Inliers_id),num_groups,r);
time_se_basis = toc;

if (num_groups  == 1)
    for j = [12,23,46,54]

    y = X_tilde(:, j);
    img_title = 'outlier';
    image_face(y, f_id, hh, ww, img_title);
    f_id = f_id + 1;

    c_hat = LSR(B_solved,y); 
    Xj_solved = B_solved* c_hat;
    img_title = 'upca - recovered';
    image_face(Xj_solved, f_id, hh, ww, img_title);
    f_id = f_id + 1;

    c_hat = LSR(Basis_estimated,y); 
    Xj_solved = Basis_estimated* c_hat;
    img_title = 'se - recovered';
    image_face(Xj_solved, f_id, hh, ww, img_title);
    f_id = f_id + 1;

    end
else
    for j = 23
        f_id = 0;
        for i = 1:num_groups

        y = X_tilde(:, sum(nn(1:i))+j);
        img_title = 'outlier';
        image_face(y, f_id, hh, ww, img_title);
        f_id = f_id + 1;


        c_hat = LSR(B_solved,y); 
        Xj_solved = B_solved* c_hat;
        img_title = 'upca - recovered';
        image_face(Xj_solved, f_id, hh, ww, img_title);
        f_id = f_id + 1;


        [X_solved] = Solve_partial_LSR(X_tilde, sum(nn(1:i))+j, Basis_estimated);
        Xj_solved = X_solved(:,sum(nn(1:i))+j);
        img_title = 'se - recovered';
        image_face(Xj_solved, f_id, hh, ww, img_title);
        f_id = f_id + 1;

        end
    end

end




