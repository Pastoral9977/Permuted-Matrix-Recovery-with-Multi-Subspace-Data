function [X_solved, outlier_class_solved] = Solve_partial_LSR(X_tilde, idOutliers, U_solved, broadcast)
    %--------------------------------------------------------------------------
    % This is the function to reconstructed the original matrix based on the
    % estimated bases U_solved and estimated outliers id idOutliers.
    % 
    %--------------------------------------------------------------------------
    % Copyright @ ...
    %--------------------------------------------------------------------------
    if nargin < 4
        broadcast = true;
    end
    
    X_solved = X_tilde;
%     if class(U_solved) == cell
        
    [a,b,c] = size(U_solved); % b = rrank, c = num_groups
    L = zeros(1,c); % To store the 'in' distance
    D = zeros(1,c); % To store the 'out' distance
    outlier_class_solved = zeros(length(idOutliers),1); % To store the outliers' classified labels

    n = length(idOutliers);
    if broadcast
        display_progress()
    end
    for i = 1:n
        j = idOutliers(i);
        y = X_tilde(:, j);
        for k = 1:c
%             Method 8 for yaleb face data
%             [~, res] = LSR_v3(U_solved(:,:,k), y, 0.5, 20); % for Synthetic data
            [~, res] = LSR_v3(U_solved(:,:,k), y, 0.5, 40); % for Hopkins data
%             [~, res] = LSR_v3(U_solved(:,:,k), y, 0.7, 20); % for yaleb face data
            L(k) = res.subspace_cos_dist;
            D(k) = 1;

        end
        es = L./D;
        [~,p] = min(es);
        outlier_class_solved(i) = p;
        c_hat = LSR(U_solved(:,:,p), y);
        X_solved(:,j) = U_solved(:,:,p) * c_hat;
        
        if broadcast
            update_progress(i, n);
        end
    end
end




