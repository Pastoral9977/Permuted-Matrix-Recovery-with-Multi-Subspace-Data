function [perm_data, outliers_ids] = permute_data(data, sr, or, s)
[D, V] = size(data);
if nargin < 4
    outliers_ids = sort(randperm(V, fix(or*V)));
    outliers_ids = reshape(outliers_ids, [1, length(outliers_ids)]);
else
    split_indices = find(diff(s));  % Example indices to split the data, representing n-1 splits

    % Append N+1 to the split indices to handle the last segment
    split_indices = reshape(split_indices, [1, length(split_indices)]);
    
    split_indices = [split_indices, V];
    num_splits = numel(split_indices);

    % Initialize variables
    outliers_ids = [];
    previous_index = 1;
    
    for i = 1:num_splits
        current_index = split_indices(i);
        segment_range = previous_index:current_index;

        num_shuffled_segment = fix(length(segment_range) * or);
        if num_shuffled_segment > 0
            outliers_ID_segment = segment_range(randperm(length(segment_range), num_shuffled_segment));
            outliers_ids = [outliers_ids, outliers_ID_segment];
        end

        previous_index = current_index;
    end
end

num_shuffled = fix(D*sr);
perm_data = data;
for j = 1:length(outliers_ids)
    x = data(:,outliers_ids(j));
    perm = 1:D;
    shuffled_index = randperm(D, num_shuffled); % Randomly choose the shuffled indices
    perm(shuffled_index) = shuffled_index(randperm(num_shuffled)); % Randomly shuffle the choosed indices
    x = x(perm); 
    perm_data(:,outliers_ids(j)) = x;
end
end