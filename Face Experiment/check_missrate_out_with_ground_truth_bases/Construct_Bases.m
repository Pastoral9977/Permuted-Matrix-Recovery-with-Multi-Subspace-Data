%% Bases construction
function Construct_Bases()
    global Involve_reconstructed_bases INIT BASES LABEL DATA
    U_gt = ground_truth_bases_construction(DATA.X_gt, LABEL.allpoints_class_gt, INIT.rrank);
    BASES.U_gt = U_gt;
    
    if Involve_reconstructed_bases
        run Inliers_Detection.m
        run Basis_Reconstrution.m
    end
    
end