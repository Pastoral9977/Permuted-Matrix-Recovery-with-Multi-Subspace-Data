% close all
% clear 

num_groups = 3; % һ��num_groups��������num_groups���ӿռ䣩
start_person_id = 1; % �ӵڼ����˿�ʼѡ
outliers_num = 16; % number of outliers out of 64 points
perm_flag = 1; % 0 for fully shuffling; 1 for partially shuffling.
rrank = 9; % low-rank�㷨����

% run Data_Loading.m
% run Initialization.m
% run Data_Construction.m
% run Inliers_Detection.m
% run Basis_Reconstrution.m
run Matrix_Reconstrution.m
run Final_Assesment.m
% run FigureFace.m