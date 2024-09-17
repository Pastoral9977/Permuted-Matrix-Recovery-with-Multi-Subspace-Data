function [perm_data, outliers_ID_in_parts, perm_labels] = permute_data(data, sr, or, s, parts, seed)
[D, V] = size(data);
if nargin < 6
    seed = randi(1000000);
end
if nargin < 5
    parts = 1;
end

% perm in different segment same part is the same.
perm_labels = cell(parts,1);
num_shuffled = fix(D*sr);
for jj = 1:parts
    perm = generate_perm(D,num_shuffled,seed*jj);                
    perm_labels{jj} = perm;
end

% Example indices to split the data, representing n-1 splits
split_indices = find(diff(s));

% Append N+1 to the split indices to handle the last segment
split_indices = reshape(split_indices,[1,length(split_indices)]);
split_indices = [split_indices, V];
num_splits = numel(split_indices);
if num_splits ~= length(unique(s))
    error('Something wrong with data permutation!')
end

% Initialize variables
previous_index = 1;
perm_data = data;

outliers_ID_in_parts = cell(parts,1);
for ii = 1:num_splits
    current_index = split_indices(ii);
    segment_range = previous_index:current_index;
    num_shuffled_segment = fix(length(segment_range)*or);
    if num_shuffled_segment > 0
        outliers_ID_segment = segment_range(sort(randperm(length(segment_range),num_shuffled_segment)));
        outliers_ID_segment_parts = split_array(outliers_ID_segment, parts);
        for jj = 1:parts
            perm = perm_labels{jj};                                     % perm in different segment same part is the same.
            ids = outliers_ID_segment_parts{jj};
            perm_data(:,ids) = perm_data(perm,ids);
            outliers_ID_in_parts{jj} = [outliers_ID_in_parts{jj}, ids];
        end
    end
    previous_index = current_index;
end

end

function perm = generate_perm(D,num_shuffled,seed)
if nargin < 2
    seed = randi(1000000);
end
rng(seed)
rp = randperm(D);
shuffled_index = rp(1:num_shuffled);
shuffled_how = randperm(num_shuffled);
shuffled_pattern= shuffled_index(shuffled_how);
perm = 1:1:D;
perm(shuffled_index) = shuffled_pattern;
end

function arrays = split_array(array, parts)
arrays = cell(1,parts);
n = length(array);
if n < parts
    error('n should be greater than k!')
end
sub = fix(n/parts);
for i = 1:parts-1
    arrays{i} = array((i-1)*sub+1:i*sub);
end
arrays{parts} = array((parts-1)*sub+1:end);
end