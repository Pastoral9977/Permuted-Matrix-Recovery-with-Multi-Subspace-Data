function [err_ratio] = EvaluateRefined_v2(X_solved, X_gt, num_groups, rrank)

[D,N,num_or] = size(X_solved);
V = N/num_groups;
M = zeros(D, V, num_groups);
X_solved_refine = X_solved;
err_ratio = zeros(1,num_or);

for p = 1:num_or
    for i = 1:num_groups
        M(:,:,i) = X_solved(:,(i-1)*V+1:i*V, p);
        [U,~,~] = svd(M(:,:,i),'econ');
        B = U(:,1:rrank);
        P = B*B';
        D = P * M(:,:,i);
        X_solved_refine(:,(i-1)*V+1:i*V,p) = D ./ vecnorm(D);
    end
err_refine = norm(X_solved_refine(:,:,p) - X_gt, 'fro');
err_ratio(p) = err_refine / norm(X_gt, 'fro');
end

end
