%% Initialization
global DATA INIT PIC
DATA.num_eachgroup = [];
for i = 1:INIT.num_groups
    [~, DATA.num_eachgroup(i)] = size(DATA.YaleB_cell{INIT.start_person_id-1+i});
end
N = sum(DATA.num_eachgroup);
DATA.N = N;

PIC.height = 48; 
PIC.width = 42; 
PIC.v_patches = 8; % v_patches: number of vertical blocks in one 'big column', #rows
PIC.h_patches = 7; % h_patches: number of horizontal blocks in one 'big row', #cols

DATA.X_gt = zeros(PIC.height*PIC.width, DATA.N); 
DATA.X_tilde = zeros(PIC.height*PIC.width, DATA.N); 