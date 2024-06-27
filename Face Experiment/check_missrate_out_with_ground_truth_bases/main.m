%--------------------------------------------------------------------------
% This function is designed to check the performance of outlier
% classification, etc missrate_out, with repect to the ground truth bases
% and solved bases.
%--------------------------------------------------------------------------
% Copyright @ XLQ
%--------------------------------------------------------------------------

%% Set the paths
close all
clear 
addpath(genpath('..'));
addpath(genpath('..\..\Copy_SSC_ADMM_v1.1'));
addpath(genpath('..\..\Functions'));


%% Set global values
global Involve_reconstructed_bases Involve_gt_bases DATA LABEL BASES RESULT INIT PIC

DATA = struct();
LABEL = struct();
BASES = struct();
RESULT = struct();
Involve_reconstructed_bases = true;
Involve_gt_bases = false;
INIT = struct();
PIC = struct();


%% Set initial parameters
rng(42)
num_groups = 12; % һ��num_groups��������num_groups���ӿռ䣩
start_person_id = 21; % �ӵڼ����˿�ʼѡ
outliers_num_each = 16; % number of outliers out of 64 points
perm_flag = 1; % 0 for fully shuffling; 1 for partially shuffling.
rrank = 8; % low-rank�㷨����

INIT.num_groups = num_groups; % һ��num_groups��������num_groups���ӿռ䣩
INIT.start_person_id = start_person_id; % �ӵڼ����˿�ʼѡ
INIT.outliers_num_each = outliers_num_each; % number of outliers out of 64 points
INIT.perm_flag = perm_flag; % 0 for fully shuffling; 1 for partially shuffling.
INIT.rrank = rrank; % low-rank�㷨����
DATA.num_groups = num_groups;


%% main
prepare_data_and_label();
Construct_Bases();
Outlier_Classification();
Assessment();







