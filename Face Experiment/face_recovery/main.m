% close all
% clear 

num_groups = 3; % 一共num_groups张人脸（num_groups个子空间）
start_person_id = 1; % 从第几个人开始选
outliers_num = 16; % number of outliers out of 64 points
perm_flag = 1; % 0 for fully shuffling; 1 for partially shuffling.
rrank = 9; % low-rank算法参数

% run Data_Loading.m
% run Initialization.m
% run Data_Construction.m
% run Inliers_Detection.m
% run Basis_Reconstrution.m
run Matrix_Reconstrution.m
run Final_Assesment.m
% run FigureFace.m