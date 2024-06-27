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
    L = zeros(1,size(U_solved,3)); % To store the 'in' distance
    D = zeros(1,size(U_solved,3)); % To store the 'out' distance
    outlier_class_solved = zeros(length(idOutliers),1); % To store the outliers' classified labels

    n = length(idOutliers);
    if broadcast
        display_progress()
    end
    for i = 1:n
        j = idOutliers(i);
        y = X_tilde(:, j);
        for k = 1:c

    % %         Method 0 for breast cancer data and scores data
    %         [c_hat] = LSR(U_solved(:,:,k), y);
    %         X_hat(:,k) = U_solved(:,:,k) * c_hat;
    %         subspace_cos_sim1 = cosine_subspace_angle(y, U_solved(:,:,k));
    %         L(k) = (1-subspace_cos_sim1);
    %         D(k) = 1;

% %             Method 1 for Hopkins data
%             [c_hat] = LSR(U_solved(:,:,k), y);
%             X_hat(:,k) = U_solved(:,:,k) * c_hat;
%             v1 = retain_smallest_part(sort(y) - sort(X_hat(:,k)));
%             L(k) = norm(v1);
%             D(k) = 1;
%             [c_out] = LSR(UU_solved(:,:,k), y);
%             X_out(:,k) = UU_solved(:,:,k) * c_out;
%             v2 = retain_smallest_part(sort(y) - sort(X_out(:,k)));
%             D(k) = norm(v2);

% %             Method 2 for Hopkins data
%             percent = 1;
%             q = 1;
%             [c_hat] = LSR(U_solved(:,:,k), y);
%             X_hat(:,k) = U_solved(:,:,k) * c_hat;
%             [c_out] = LSR(UU_solved(:,:,k), y);
%             X_out(:,k) = UU_solved(:,:,k) * c_out;
%             subspace_cos_sim1 = cosine_subspace_angle(y, U_solved(:,:,k));
%             v1 = retain_smallest_part(sort(y) - sort(X_hat(:,k)), percent);
%             v2 = retain_smallest_part(sort(y) - sort(X_out(:,k)), percent);
%             u1 = retain_smallest_part(X_hat(:,k), q, y);
%             cos_sim1 = calculate_cosine_similarity(y, u1);
%             L(k) = vecnorm(v1)*(1-subspace_cos_sim1)*(1-cos_sim1);
%             D(k) = vecnorm(v2);

%     %         Method 3
%             [c_hat] = LSR(U_solved(:,:,k), y);
%             X_hat(:,k) = U_solved(:,:,k) * c_hat;
%             [c_out] = LSR(UU_solved(:,:,k), y);
%             X_out(:,k) = UU_solved(:,:,k) * c_out;
%             L(k) = calculate_cosine_distance(sort(y), sort(X_hat(:,k)));
%             D(k) = calculate_cosine_distance(sort(y), sort(X_out(:,k)));

% %           Method 4
%             [c_hat] = LSR(U_solved(:,:,k), y);
%             X_hat(:,k) = U_solved(:,:,k) * c_hat;
%             [c_out] = LSR(UU_solved(:,:,k), y);
%             X_out(:,k) = UU_solved(:,:,k) * c_out;
% %             L(k) = vecnorm(sort(y) - sort(X_hat(:,k)))*(1-calculate_cosine_similarity(sort(y), sort(X_hat(:,k))));
% %             D(k) = vecnorm(sort(y) - sort(X_out(:,k)))*(1-calculate_cosine_similarity(sort(y), sort(X_out(:,k))));
%             L(k) = vecnorm(sort(y) - sort(X_hat(:,k)));
%             D(k) = vecnorm(sort(y) - sort(X_out(:,k)));

    %         Method 5 for yaleb face data
%             [c_hat] = LSR(U_solved(:,:,k), y);
%             X_hat(:,k) = U_solved(:,:,k) * c_hat;
%             subspace_cos_sim1 = cosine_subspace_angle(y, U_solved(:,:,k));
%             L(k) = (1-subspace_cos_sim1);
%             D(k) = 1;
    %         subspace_cos_sim2 = cosine_subspace_angle(y, UU_solved(:,:,k));
    %         D(k) = subspace_cos_sim2;

    % %         Method 6 for yaleb face data
    %         q = 1;
    %         [c_hat] = LSR(U_solved(:,:,k), y);
    %         X_hat(:,k) = U_solved(:,:,k) * c_hat;
    %         subspace_cos_sim1 = cosine_subspace_angle(y, U_solved(:,:,k));
    %         u1 = retain_smallest_part(X_hat(:,k), q, y);
    %         cos_sim1 = calculate_cosine_similarity(y, u1);
    %         L(k) = (1-subspace_cos_sim1)*(1-cos_sim1);
    %         D(k) = 1;

    % %         Method 7 for yaleb face data
    %         percent = 1; q = 1;
    %         [c_hat] = LSR(U_solved(:,:,k), y);
    %         X_hat(:,k) = U_solved(:,:,k) * c_hat;
    %         subspace_cos_sim1 = cosine_subspace_angle(y, U_solved(:,:,k));
    %         v1 = retain_smallest_part(sort(y) - sort(X_hat(:,k)), percent);
    %         u1 = retain_smallest_part(X_hat(:,k), q, y);
    %         cos_sim1 = calculate_cosine_similarity(y, u1);
    %         a = vecnorm((v1).^0.1);
    %         L(k) = a*(1-subspace_cos_sim1)*(1-cos_sim1);
    %         D(k) = 1;
    % %         L(k) = vecnorm(v1)*(1-subspace_cos_sim1);
    % %         L(k) = vecnorm(v1);
    % %         L(k) = vecnorm(v1)*(1-cos_sim1);
    
%             Method 8 for yaleb face data
%             [c_hat, res] = LSR_v3(U_solved(:,:,k), y, 0.5, 20); % for Synthetic data
%             [~, res] = LSR_v3(U_solved(:,:,k), y, 0.5, 40); % for Hopkins data
            [~, res] = LSR_v3(U_solved(:,:,k), y, 0.7, 20); % for yaleb face data
%             X_hat(:,k) = U_solved(:,1:r_refine,k) * c_hat;
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




