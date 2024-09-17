function [err_ratio] = Evaluate(X_solved, X_gt)
    err_refine = norm(X_solved - X_gt, 'fro');
    err_ratio = err_refine/norm(X_gt, 'fro');
end