function [X_tilde] = permute_face_randomly(X, outliers_ids, h, w, v_patches, h_patches, perm_flag)
    % This function generate the observed data corrupted by unknown permutations on some columns.
    %   v_patches: number of vertical blocks in one 'big column', #rows
    %   h_patches: number of horizontal blocks in one 'big row', #cols
    %   X_tilde: X_gt with shuffling
    %   perm_flag = 0 when patches are fully shuffled; 
    %   perm_flag = 1 when patches are partially shuffled;

    [m, n] = size(X); 
    edge = fix(sqrt(v_patches*h_patches*0.4));
    num_fixed = v_patches*h_patches - edge^2;
    hw = h*w; 
    if (m == hw) && (mod(h, v_patches)==0) && (mod(w, h_patches)==0)
        num_patches = v_patches * h_patches;
        vb_size = round(h/v_patches); % #pixels in one column of each block
        hb_size = round(w/h_patches); % #pixels in one row of each block
        row_start_id = ones(v_patches, h_patches);
        row_end_id = ones(v_patches, h_patches);
        col_start_id = ones(v_patches, h_patches);
        col_end_id = ones(v_patches, h_patches);
        part_id = cell(1, num_patches);
        kl = 0;
        for k = 1: v_patches
            for l = 1: h_patches
                kl = kl + 1;
                row_start_id(k, l) = (k-1) * vb_size + 1;
                row_end_id(k, l) = k * vb_size;
                col_start_id(k, l) = (l-1) * hb_size + 1;
                col_end_id(k, l) = l * hb_size;
                part_id{kl} = [k, l];
            end
        end

    %     num_outliers = fix(outlier_ratio * n); % # permuted columns
        X_tilde = X;

        for j = outliers_ids

            if perm_flag==0 && h_patches==2 && v_patches ==2
                perm_pattern = [4 3 1 2]; %[2 4 1 3]
            elseif perm_flag==0 && h_patches==4 && v_patches ==4
                perm_pattern = randperm(v_patches * h_patches);%[8 4 6 2 1 9 5 3 7]; %
            elseif perm_flag==0
                perm_pattern = randperm(v_patches*h_patches);
            elseif perm_flag==1
            perm_start_row = randperm(v_patches-edge+1,1);
            perm_start_col = randperm(h_patches-edge+1,1);
            perm_end_row = perm_start_row+edge-1;
            perm_end_col = perm_start_col+edge-1;

            perm_pattern = 1:num_patches;
            matrix_process = reshape(perm_pattern,[h_patches, v_patches]);
            matrix_process = matrix_process';
            matrix_perm = matrix_process(perm_start_row:perm_end_row, ...
                perm_start_col: perm_end_col);
            vec = reshape(matrix_perm, [1,numel(matrix_perm)]);
            vec = vec(randperm(length(vec)));
            matrix_perm = reshape(vec,[edge, edge]);
            matrix_process(perm_start_row: perm_end_row, perm_start_col:perm_end_col)...
                        = matrix_perm;
            perm_pattern = reshape(matrix_process', [1,num_patches]);
            end

            Xj = X(:, j);
            Xj_mat = reshape(Xj, [h, w]);
            Xj_shuffled = Xj_mat;
            kl = 0;
            for k = 1: v_patches
                for l = 1: h_patches
                    kl = kl + 1;
                    kl_tuple = part_id{perm_pattern(kl)};
                    k_prime = kl_tuple(1);
                    l_prime = kl_tuple(2);
                    Xj_shuffled( row_start_id(k, l): row_end_id(k, l), col_start_id(k, l): col_end_id(k, l)) ...
                        = Xj_mat( row_start_id(k_prime, l_prime): row_end_id(k_prime, l_prime), ...
                        col_start_id(k_prime, l_prime): col_end_id(k_prime, l_prime));
                end
            end
            X_tilde(:, j) = reshape(Xj_shuffled, [hw, 1]);
        end
    %     idOutliers_binary = [ones(1, num_outliers), zeros(1, n - num_outliers)]; % n-dimensional logical vector to indicate outliers
    %     idOutliers = find(idOutliers_binary .* (1:1:n));
    else
        error('m is not hw, or the number of patches is not suitable to h*w');
    end

end