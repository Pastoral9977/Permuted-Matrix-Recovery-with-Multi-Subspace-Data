%% Initialization
% 设定实验数据的基本参数
global DATA INIT PIC
DATA.num_eachgroup = [];
for i = 1:INIT.num_groups
    [~, DATA.num_eachgroup(i)] = size(DATA.YaleB_cell{INIT.start_person_id-1+i});
end
N = sum(DATA.num_eachgroup); % 图片总数
DATA.N = N;

% 压缩与乱序的参数
PIC.height = 48; % 压缩后图片的高
PIC.width = 42; % 压缩后图片的宽
PIC.v_patches = 48; % v_patches: number of vertical blocks in one 'big column', #rows
PIC.h_patches = 42; % h_patches: number of horizontal blocks in one 'big row', #cols

DATA.X_gt = zeros(PIC.height*PIC.width, DATA.N); % 初始化压缩后乱序前的样本矩阵
DATA.X_tilde = zeros(PIC.height*PIC.width, DATA.N); % 初始化压缩后乱序后的样本矩阵