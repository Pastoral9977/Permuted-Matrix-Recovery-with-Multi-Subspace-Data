%% Matrix reconstrution
global DATA LABEL INIT BASES
tic
fprintf('Matrix reconstrution...\n')
if (INIT.perm_flag == 0)
    if size(BASES.U_esti,2) > 4
        BASES.U_esti = BASES.U_esti(:,1:4,:);
    end
    [DATA.X_solved, ~, LABEL.outlier_class_solved] = Solve_all(DATA.X_tilde, LABEL.Outliers_id, BASES.U_esti);
else
    [DATA.X_solved, ~, LABEL.outlier_class_solved] = Solve_partial_LSR(DATA.X_tilde, LABEL.Outliers_id, BASES.U_esti);
end
fprintf('\tOver, costing %0.2fs\n', toc)

