close all
clear
addpath(genpath('C:\Users\Pastoral\Desktop\Copy_SSC_ADMM_v1.1'));
addpath(genpath('C:\Users\Pastoral\Desktop\Unlabeled_PCA_Code'));
addpath(genpath('C:\Users\Pastoral\Desktop\Pursuit\Functions'));

load('YaleB_cell');
load('YaleBCrop025');
num_groups = 10;

h = 192; w = 168; 
hh = 48; ww = 42;
[~,n] = size(YaleB_cell{1,1}); n1 = 40;
lambda = 0.05; outliers_num = 1;
v_patches = 48;
h_patches = 42;
perm_flag = 1;

X1 = zeros(hh*ww,n*num_groups);
X2 = zeros(hh*ww,n1*num_groups);
X3 = zeros(hh*ww,n1*num_groups);

for i = 1:num_groups
    X1(:,(i-1)*n+1:i*n) = DSP(YaleB_cell{1,i},h,w,hh,ww);
    p = sort(randperm(n,n1));
    X2(:,(i-1)*n1+1:i*n1) = DSP(YaleB_cell{1,i}(:,p),h,w,hh,ww);
    outliers_id = randperm(n1,outliers_num);
    X3(:,(i-1)*n1+1:i*n1) = permute_face(X2(:,(i-1)*n1+1:i*n1), ...
        outliers_id, hh, ww, v_patches, h_patches, perm_flag);
end

Y1 = rpca(X1,lambda);
Y2 = rpca(X2,lambda);
Y3 = rpca(X3,lambda);

t1 = zeros(num_groups*n,1);
t2 = zeros(num_groups*n1,1);
for i = 1:num_groups
    t1((i-1)*n+1:i*n,1) = i;
    t2((i-1)*n1+1:i*n1,1) = i;
end
t3 = t2;

r = 0; affine = false; outlier = true; rho = 1;alpha = 20;
[missrate1_in,~,~,~] = SSC(Y1,r,affine,alpha,outlier,rho,t1);
[missrate2_in,~,~,~] = SSC(Y2,r,affine,alpha,outlier,rho,t2);
[missrate3_in,~,~,~] = SSC(Y3,r,affine,alpha,outlier,rho,t3);
[missrate4_in,~] = SSC(X1,r,affine,alpha,outlier,rho,t1);
[missrate5_in,~] = SSC(X2,r,affine,alpha,outlier,rho,t2);
[missrate6_in,~] = SSC(X3,r,affine,alpha,outlier,rho,t3);

Z = [];
for i = 1:num_groups
    t((i-1)*n+1:i*n,1) = i;
    Z = [Z Y(:,:,i)];
end

[missrate,~] = SSC(Z,r,affine,alpha,outlier,rho,t);

close all
W = X1-Z;
id = 54;
image_face(Y1(:, id), 3, hh, ww, 'mine');
image_face(Z(:, id), 4, hh, ww, 'their');
image_face(W(:, id), 5, hh, ww, 'difference');







