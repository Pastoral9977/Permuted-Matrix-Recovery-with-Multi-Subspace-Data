%% Data Loading
% ����·��
global DATA PIC
addpath(genpath('..\Copy_SSC_ADMM_v1.1'));
addpath(genpath('..\Functions'));

% ��������
fprintf('Data Loading...')
tic
load('YaleB_cell'); % ���غ��YaleB_cell����Ϊ1*38��cell�������ͣ���38���ˣ�ÿ��cellΪͬһ�����Ķ��źڰ���Ƭ������Ϊ(ÿ������ͼ������һά����ά��*��Ƭ����)
h = 192; % ͼƬ�ĸ�
w = 168; % ͼƬ�Ŀ�

DATA.YaleB_cell = YaleB_cell;
PIC.ori_height = h;
PIC.ori_width = w;

fprintf('Over, costing %0.2fs\n\n', toc)