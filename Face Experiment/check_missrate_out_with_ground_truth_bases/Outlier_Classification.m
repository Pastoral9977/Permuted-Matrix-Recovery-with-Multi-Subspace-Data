%% Outlier classification   
function Outlier_Classification()
    global DATA LABEL BASES Involve_reconstructed_bases Involve_gt_bases
    
    tic
    fprintf('Outlier classification...\n')
    
    if Involve_gt_bases
        % With ground truth bases
        fprintf('\twith ground truth bases...\n\t')
        [X_solved_gt, outlier_class_solved_gtbases] = Solve_partial_LSR(DATA.X_tilde, LABEL.outliers_ID, BASES.U_gt);
        DATA.X_solved_gt = X_solved_gt;
        LABEL.outlier_class_solved_gtbases = outlier_class_solved_gtbases;
    end
    
    % With reconstructed bases
    if Involve_reconstructed_bases
        fprintf('\twith reconstruced bases...\n\t')
        [X_solved_reconstructed, outlier_class_solved_reconstructedbases] = Solve_partial_LSR(DATA.X_tilde, LABEL.Outliers_id, BASES.U_esti);
        DATA.X_solved_reconstructed = X_solved_reconstructed;
        LABEL.outlier_class_solved_reconstructedbases = outlier_class_solved_reconstructedbases;
    end
    
    fprintf('\tOver, costing %0.2fs\n\n', toc)
end