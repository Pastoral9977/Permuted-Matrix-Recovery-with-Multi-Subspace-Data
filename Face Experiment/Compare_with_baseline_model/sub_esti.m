function [Basis_estimated] = sub_esti(X,n,r)
    %% Subspace Clustering
    [grps] = SSC_with_numgroups_only(X,n);
    D = size(X,1);
    O_id = 1:size(X,2);
    num_groups = max(grps);
    
    %% Subscript Matching
    X_selected_inlier_groups = cell(1, num_groups);
    for k = 1:num_groups
        kth_group_ids = O_id(grps == k);
        X_selected_inlier_groups{k} = X(:, kth_group_ids); 
    end

    %% Subspace Estimation
    Basis_estimated = zeros(D, r, n);
    for k = 1:n
        [U,~,~] = svd(X_selected_inlier_groups{k});
        Basis_estimated(:,:,k) = U(:,1:r) ;
    end

end








