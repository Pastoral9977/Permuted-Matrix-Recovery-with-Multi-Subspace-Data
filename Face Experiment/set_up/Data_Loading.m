%% Data Loading
% 加载路径
global DATA PIC
addpath(genpath('..\Copy_SSC_ADMM_v1.1'));
addpath(genpath('..\Functions'));

% 加载数据
fprintf('Data Loading...')
tic
load('YaleB_cell'); % 加载后的YaleB_cell数据为1*38的cell数据类型，即38个人，每个cell为同一人脸的多张黑白照片，具体为(每张人脸图被拉成一维向量维度*照片数量)
h = 192; % 图片的高
w = 168; % 图片的宽

DATA.YaleB_cell = YaleB_cell;
PIC.ori_height = h;
PIC.ori_width = w;

fprintf('Over, costing %0.2fs\n\n', toc)