function [X_tilde, target_perm_index] = permute_face_randomly_MRUC_for_perm_indices(X, outliers_ids, n, m, np, mp, sample_id, missing_ratio, shuffled_ratio, seed)
    % This function generate the observed data corrupted by unknown permutations on some columns.
    %   v_patches: number of vertical blocks in one 'big column', #rows
    %   h_patches: number of horizontal blocks in one 'big row', #cols
    %   X_tilde: X_gt with shuffling
    %   perm_flag = 0 when patches are fully shuffled; 
    %   perm_flag = 1 when patches are partially shuffled with shuffled ratio = 0.4;
    
    if nargin < 10
        seed = randi(100000000);
    end
    if nargin < 9
        shuffled_ratio = 0.4;
    end
    if nargin < 8
        missing_ratio = 0.0;
    end
    if nargin < 7
        error('Please use permute_face_randomly_v2 since you did not input sample_id')
    end
    [M, N] = size(X);
    X_tilde = X;
    for ii = 1:length(outliers_ids)
        id = outliers_ids(ii);
        x = X(:,id);
        img = reshape(x, [n, m]);
        rng(seed*ii)
        perm_num = round(shuffled_ratio*np*mp);
        perm_index = sort(randperm(np*mp, perm_num));
        if id == sample_id
            target_perm_index = perm_index;
        end
        perm_img = permute_image_with_given_index(img, np, mp, perm_index, missing_ratio, seed);
        X_tilde(:, id) = perm_img(:);
    end
end