function [X_solved, time, IDX_solved] = Solve_partial_L1RR(X_tilde, idOutliers, U_solved)

tic
X_solved = X_tilde;
[a,b,c] = size(U_solved);
b = fix(b/2);
X_hat = zeros(size(X_tilde,1),size(U_solved,3));
X_out = zeros(size(X_tilde,1),size(U_solved,3));
L = zeros(1,size(U_solved,3)); 
D = zeros(1,size(U_solved,3)); 
IDX_solved = zeros(length(idOutliers),1);
UU_solved = zeros(a,(c-1)*b,c);M = [];
for i = 1:c
    for k = setdiff(1:c,i)
        M = [M U_solved(:,1:b,k)];
    end
    UU_solved(:,:,i) = M;
    M = [];
end


for i = 1:length(idOutliers)
    %Since we may not be clear about whether idOutliers is a column of row,
    %we use function max() to eliminate the discrepency.
    j = idOutliers(i);
    y = X_tilde(:, j);
    for k = 1:size(U_solved,3)
        [c_hat] = L1_RR(U_solved(:,:,k), y);
        [c_out] = L1_RR(UU_solved(:,:,k), y);
        X_hat(:,k) = U_solved(:,:,k) * c_hat;
        X_out(:,k) = UU_solved(:,:,k) * c_out;
        L(k) = vecnorm(sort(y) - sort(X_hat(:,k)));
        D(k) = vecnorm(sort(y) - sort(X_out(:,k)));
    end
    es = L./D;
    [~,p] = min(es);
    X_solved(:,j) = X_hat(:,p);
    IDX_solved(i) = p;
    
    str = ['Solving process has been finished by ',num2str(i*100/length(idOutliers)), ...
        '%'];
%   disp(str);
    
end
time = toc;
end