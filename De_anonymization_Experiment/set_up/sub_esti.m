function [U_esti] = sub_esti(X,n,r)
    %% Subspace Clustering
    [grps] = SSC_with_numgroups_only(X,n);
    D = size(X,1);
    all_ids = 1:size(X,2);
    %% Subscript Matching
    label_hat = zeros(size(X,2), 2);
    X_selected_inlier_groups = cell(1, n);
    index = 1;
    for k = 1:n
        kth_group_ids = all_ids(grps == k);
        X_selected_inlier_groups{k} = X(:,all_ids(grps == k)); 
        
        num_kth_group_ids = length(kth_group_ids);
        label_hat(index:index+num_kth_group_ids-1, :) = [kth_group_ids', k*ones(num_kth_group_ids, 1)];
        index = index + num_kth_group_ids;
    end

    %% Subspace Estimation
    U_esti = zeros(D, r, n);
    for k = 1:n
        [U,~,~] = svd(X_selected_inlier_groups{k});
        U_esti(:,:,k) = U(:,1:r) ;
    end
end








