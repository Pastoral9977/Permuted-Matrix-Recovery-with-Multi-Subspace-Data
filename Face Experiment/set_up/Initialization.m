%% Initialization
% �趨ʵ�����ݵĻ�������
global DATA INIT PIC
DATA.num_eachgroup = [];
for i = 1:INIT.num_groups
    [~, DATA.num_eachgroup(i)] = size(DATA.YaleB_cell{INIT.start_person_id-1+i});
end
N = sum(DATA.num_eachgroup); % ͼƬ����
DATA.N = N;

% ѹ��������Ĳ���
PIC.height = 48; % ѹ����ͼƬ�ĸ�
PIC.width = 42; % ѹ����ͼƬ�Ŀ�
PIC.v_patches = 48; % v_patches: number of vertical blocks in one 'big column', #rows
PIC.h_patches = 42; % h_patches: number of horizontal blocks in one 'big row', #cols

DATA.X_gt = zeros(PIC.height*PIC.width, DATA.N); % ��ʼ��ѹ��������ǰ����������
DATA.X_tilde = zeros(PIC.height*PIC.width, DATA.N); % ��ʼ��ѹ������������������