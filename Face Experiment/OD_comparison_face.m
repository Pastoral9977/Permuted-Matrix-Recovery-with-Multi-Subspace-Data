close all
clear 
addpath(genpath('..\Copy_SSC_ADMM_v1.1'));
addpath(genpath('..\Functions'));

load('YaleB_cell');
num_groups = 5;

h = 192; w = 168;
hh = 48; ww = 42;
n = [];
sp = 1;
for i = 1:num_groups
    [~,n(i)] = size(YaleB_cell{1,sp-1+i});
end
N = sum(n); 

outliers_num = 16; % number of outliers out of 64 points
rrank = 18;
v_patches = 48;
h_patches = 42;
% perm_pattern = 0 for fully shuffling; 1 for partially shuffling.
perm_flag = 1;
detect_flag = 0;


X_raw = zeros(hh*ww,N);
X_tilde = zeros(hh*ww,N);
outliers_ID = [];

nn = [0 n];
for i = 1:num_groups
    X_raw(:,sum(nn(1:i))+1:sum(nn(1:i+1))) = DSP(YaleB_cell{1,sp-1+i},h,w,hh,ww);
    
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

%% SE_outliers_detection
affine = false ; alpha = 1e10;
tic
[Inliers_id, Outliers_id] = Outlier_Detect(X_tilde, detect_flag, alpha);

err_in = length(setdiff(Inliers_id,inliers_ID))/length(Inliers_id);
err_out = length(setdiff(Outliers_id,outliers_ID))/length(Outliers_id);
str = ['For SE: ', '(In_out, Out_in) = (', num2str(length(setdiff(Inliers_id,inliers_ID))), ',',...
    num2str(length(setdiff(Outliers_id,outliers_ID))), ')',','...
    ' err_in = ' num2str(err_in) ', err_out = ' num2str(err_out), ... 
    ', running time = ', num2str(toc)];
disp(str);


%% OP_outliers_detection
lambda = 1.00;
tic
[L, C] = rpca_OP(X_tilde,w*8,lambda);
Outliers_id = find(C(1,:));
Inliers_id = setdiff(1:N, Outliers_id);

err_in = length(setdiff(Inliers_id,inliers_ID))/length(Inliers_id);
err_out = length(setdiff(Outliers_id,outliers_ID))/length(Outliers_id);

str = ['For OP: ', '(In_num, Out_num) = (', num2str(length(Inliers_id)), ...
    ',', num2str(length(Outliers_id)), ')', ', err_in = ' num2str(err_in) ,...
    ', err_out = ' num2str(err_out), ', running time = ', num2str(toc)];
disp(str);





