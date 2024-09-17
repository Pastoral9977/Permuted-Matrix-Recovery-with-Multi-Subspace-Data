function [B1s, A1s, Ans, A1_row_inds, An_row_inds, oracle_images, permutation_matrices, test_inds, test_labels, norm_constants, raw_images, ...
          B1_TOTAL, A1_TOTAL, An_TOTAL, A1_row_ind_TOTAL, An_row_ind_TOTAL, permutation_matrix, norm_constant_TOTAL] ...
    = generate_data(D, V, ranks, num_groups, shuffled_ratio, outlier_ratio, parts, missing_percent, snr, seed)

% Set default values for optional parameters
if nargin < 10, seed = randi(10000000); end
if nargin < 9, snr = 0; end
if nargin < 8, missing_percent = 0; end
if nargin < 7, parts = 1; end
if nargin < 6, outlier_ratio = 0.5; end
if nargin < 5, shuffled_ratio = 0.5; end

% Calculate the number of columns to permute and the number of shuffled columns per group
perm_num = fix(D * shuffled_ratio);
shuffle_columns_per_group = fix(V * outlier_ratio);

assert(mod(shuffle_columns_per_group, parts) == 0, 'Error: shuffle_columns_per_group is not divisible by parts.');

% Initialize cell arrays to hold the results
B1s = cell(1, num_groups);
A1s = cell(1, num_groups);
Ans = cell(1, num_groups);
A1_row_inds = cell(1, num_groups);
An_row_inds = cell(1, num_groups);
oracle_images = cell(1, num_groups);
permutation_matrices = cell(1, num_groups);
test_inds = cell(1, num_groups);
test_labels = cell(1, num_groups);
norm_constants = cell(1, num_groups);
raw_images = cell(1, num_groups);
outlier_perm_belongings = cell(1, num_groups);

% Generate data for each group
for jj = 1:num_groups
    [B1s{jj}, A1s{jj}, Ans{jj}, A1_row_inds{jj}, An_row_inds{jj}, oracle_images{jj},...
        permutation_matrices{jj}, test_inds{jj}, test_labels{jj}, norm_constants{jj}, raw_images{jj}, outlier_perm_belongings{jj}] ...
        = generate_matrix(D, ranks, parts, missing_percent, shuffle_columns_per_group, perm_num, V, snr, seed*jj, seed);
end

% Generate B1 for whole matrix
B1_TOTAL = zeros(D, num_groups*size(B1s{1},2));
for jj = 1:num_groups
    B1_TOTAL(:,(jj-1)*size(B1s{jj},2)+1:jj*size(B1s{jj},2)) = B1s{jj};
end

% Generate permutation matrix for whole matrix
if num_groups>1
    assert(norm(permutation_matrices{1}{1}-permutation_matrices{2}{1})==0, 'Something wrong with permutation matrix!')
end
permutation_matrix = permutation_matrices{1};

% Generate A1 for whole matrix
A1_row_ind_TOTAL = cell(D, 1);
A1_TOTAL = cell(D, 1);
for ii = 1:D
    for jj = 1:num_groups
        A1_TOTAL{ii} = [A1_TOTAL{ii}, A1s{jj}{ii}];
        A1_row_ind_TOTAL{ii} = [A1_row_ind_TOTAL{ii}, A1_row_inds{jj}{ii}+maxWithDefault(A1_row_ind_TOTAL{ii})];
    end
end

% Generate An for whole matrix
An_row_ind_TOTAL = cell(D, parts);
An_TOTAL = cell(D, parts);
for kk = 1:parts
    for ii = 1:D
        for jj = 1:num_groups
            An_TOTAL{ii, kk} = [An_TOTAL{ii, kk}, Ans{jj}{ii, kk}];
            An_row_ind_TOTAL{ii, kk} = [An_row_ind_TOTAL{ii, kk}, An_row_inds{jj}{ii, kk}+maxWithDefault(An_row_ind_TOTAL{ii, kk})];
        end
    end
end

% Generate norm_constant for whole matrix
cons = 0;
for jj = 1:num_groups
    cons = cons + norm_constants{jj};
end
norm_constant_TOTAL = cons/num_groups;

end

function result = maxWithDefault(array)
    if isempty(array)
        result = 0;
    else
        result = max(array);
    end
end