%% Data Loading
global DATA PIC
fprintf('Data Loading...')
tic
load('YaleB_cell'); 
h = 192; 
w = 168; 

DATA.YaleB_cell = YaleB_cell;
PIC.ori_height = h;
PIC.ori_width = w;

fprintf('Over, costing %0.2fs\n\n', toc)