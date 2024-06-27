close all
clear 
addpath(genpath('C:\Users\Pastoral\Desktop\Copy_SSC_ADMM_v1.1'));
addpath(genpath('C:\Users\Pastoral\Desktop\Unlabeled_PCA_Code'));
addpath(genpath('C:\Users\Pastoral\Desktop\Pursuit\Functions'));

load('YaleB_cell');
num_groups = 8;

h = 192; w = 168; 
hh = 48; ww = 42;
[m,n] = size(YaleB_cell{1,1});


rank = 4;
k = 20; % number of outliers out of 64 points
v_patches = 24;
h_patches = 21;
% perm_pattern = 0 for fully shuffling; 1 for partially shuffling.
perm_flag = 1;

X_raw = zeros(hh*ww,n*num_groups);
X_tilde = zeros(hh*ww,n*num_groups);
outliers_ID = [];

for i = 1:num_groups
    X_raw(:,(i-1)*n+1:i*n) = DSP(YaleB_cell{1,i},h,w,hh,ww);
    
    shuffled_ids = randperm(n);
    outliers_id = shuffled_ids(1:k);
    if length(setdiff(outliers_id, 1))== k
        outliers_id(1) = 1;
    end
    if length(setdiff(outliers_id, 23))==k
        outliers_id(2) = 23;
    end
    [X] = permute_face(X_raw(:,(i-1)*n+1:i*n), outliers_id, hh, ww, v_patches, h_patches, perm_flag);
    X_tilde(:,(i-1)*n+1:i*n) = X;
    outliers_ID = [outliers_ID, sort(outliers_id)+(i-1)*n];

end

inliers_ID = setdiff(1:num_groups*n,outliers_ID);


%SSC
affine = false ;alpha = 10^10;

%Round 1
Z = admmLasso_mat_func(X_tilde, affine, alpha);
%想办法解释不标准化的效果比标准化要好（SSC的代码似乎也没有标准化）
value = vecnorm(Z,1);
T = kmeans(value',2);
I = find(T == 1); J = find(T == 2);
I_mean = mean(value(:,I));J_mean = mean(value(:,J));

if (I_mean < J_mean)
    Inliers_id = J;Outliers_id = I;
else
    Inliers_id = I;Outliers_id = J;
end

err_in = length(setdiff(Inliers_id,inliers_ID))/length(Inliers_id);
err_out = length(setdiff(Outliers_id,outliers_ID))/length(Outliers_id);

%Round 2
Z1 = admmLasso_mat_func(X_tilde(:,Inliers_id), affine, alpha);
value1 = vecnorm(Z1,1);
Q = kmeans(value1',2);
I1 = find(Q == 1); J1 = find(Q == 2);
I1_mean = mean(value1(:,I1));J1_mean = mean(value1(:,J1));

if (I1_mean < J1_mean)
    Inliers1_id = J1;Outliers1_id = I1;
else
    Inliers1_id = I1;Outliers1_id = J1;
end

Inliers_id = Inliers_id(setdiff(1:length(Inliers_id),Outliers1_id));
Outliers_id = setdiff(1:size(X_tilde,2),Inliers_id);

err1_in = length(setdiff(Inliers_id,inliers_ID))/length(Inliers_id);
err1_out = length(setdiff(Outliers_id,outliers_ID))/length(Outliers_id);

%Round 3
Z2 = admmLasso_mat_func(X_tilde(:,Outliers_id), affine, alpha);
value2 = vecnorm(Z2,1);
R = kmeans(value2',2);
I2 = find(R == 1); J2 = find(R == 2);
I2_mean = mean(value2(:,I2));J2_mean = mean(value2(:,J2));

if (I2_mean < J2_mean)
    Inliers2_id = J2;Outliers2_id = I2;
else
    Inliers2_id = I2;Outliers2_id = J2;
end

Outliers_id = Outliers_id(setdiff(1:length(Outliers_id),Inliers2_id));
Inliers_id = setdiff(1:size(X_tilde,2),Outliers_id);

err2_in = length(setdiff(Inliers_id,inliers_ID))/length(Inliers_id);
err2_out = length(setdiff(Outliers_id,outliers_ID))/length(Outliers_id);










img_title = 'outlier';
image_face(X_tilde(:, 23), 1, hh, ww, img_title);
img_title = 'inlier';
image_face(X_raw(:, 23), 10, hh, ww, img_title);

Z_s = Z(:,[inliers_ID,outliers_ID]);
value = vecnorm(Z_s,1);
figure
plot(value);grid on; hold on;
line([length(inliers_ID) length(inliers_ID)],[0,2],'linestyle','--', 'Color','r', 'LineWidth', 1);
title('1-norm of the consequential coefficients for clustering')
text(0.4*length(inliers_ID),0.5,'inliers','color','red','FontSize',16);
text(n*num_groups-0.6*length(outliers_ID),0.5,'outliers','color','red','FontSize',16);










