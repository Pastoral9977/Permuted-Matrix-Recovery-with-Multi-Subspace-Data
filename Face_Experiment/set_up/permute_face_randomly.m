function [X_tilde] = permute_face_randomly(X, outliers_ids, n, m, np, mp, missing_ratio, shuffled_ratio, seed)
    % This function generate the observed data corrupted by unknown permutations on some columns.
    %   v_patches: number of vertical blocks in one 'big column', #rows
    %   h_patches: number of horizontal blocks in one 'big row', #cols
    %   X_tilde: X_gt with shuffling
    %   perm_flag = 0 when patches are fully shuffled; 
    %   perm_flag = 1 when patches are partially shuffled with shuffled ratio = 0.4;
    
    if nargin < 9
        seed = randi(100000000);
    end
    if nargin < 8
        shuffled_ratio = 0.4;
    end
    if nargin < 7
        missing_ratio = 0.0;
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
        perm_img = permute_image_with_given_index(img, np, mp, perm_index, missing_ratio, seed);
        X_tilde(:, id) = perm_img(:);
    end
end

function permuted_img = permute_image_with_given_index(img, np, mp, perm_index, missing_ratio, seed)
    % Input:
    % img - the original image with dimension n*m
    % np - number of partitions along the height
    % mp - number of partitions along the width
    % perm_index - a vector containing the linear indices of the blocks that need to be permuted
    % missing_ratio - ratio of NaN (or 0) entries.
    % seed - random seed
    
    if nargin < 6
        seed = randi(10000000000);
    end
    rng(seed)
    if nargin < 5
        missing_ratio = 0;
    end
    
    % Get the size of the image
    [n, m] = size(img);
    
    % Calculate the size of each block
    block_height = floor(n / np);
    block_width = floor(m / mp);
    
    % Initialize cell array to store blocks
    blocks = cell(np, mp);
    
    % Extract blocks
    for i = 1:np
        for j = 1:mp
            row_start = (i-1)*block_height + 1;
            row_end = i*block_height;
            col_start = (j-1)*block_width + 1;
            col_end = j*block_width;
            blocks{i, j} = img(row_start:row_end, col_start:col_end);
        end
    end
    
    % Convert blocks to a linear array in column-major order
    blocks_linear = reshape(blocks, [], 1);
    
    % Ensure that all perm_index positions are permuted and not self-mapped
    permuted_indices = perm_index;
    while true
        permuted_indices = perm_index(randperm(length(perm_index)));
        if all(permuted_indices ~= perm_index)
            break;
        end
    end
    
    % Permute the blocks according to perm_index and permuted_indices
    permuted_blocks_linear = blocks_linear;
    permuted_blocks_linear(perm_index) = blocks_linear(permuted_indices);

    
    % Reshape permuted_blocks_linear back to a 2D cell array
    permuted_blocks = reshape(permuted_blocks_linear, np, mp);
    
    % Initialize the permuted image
    permuted_img = zeros(n, m);
    
    % Combine the permuted blocks into the final image
    for i = 1:np
        for j = 1:mp
            row_start = (i-1)*block_height + 1;
            row_end = i*block_height;
            col_start = (j-1)*block_width + 1;
            col_end = j*block_width;
            permuted_img(row_start:row_end, col_start:col_end) = permuted_blocks{i, j};
        end
    end
    permuted_img = add_missing_values(permuted_img, missing_ratio);    
end

function Am = add_missing_values(A, missing_ratio)
[n, m] = size(A);
Am = A(:);
total_elements = numel(Am);
num_missing = round(total_elements * missing_ratio);
missing_indices = randperm(total_elements, num_missing);
Am(missing_indices) = 0;
Am = reshape(Am, [n, m]);
end
