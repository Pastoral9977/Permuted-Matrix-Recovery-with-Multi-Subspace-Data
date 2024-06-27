function [outcome] = shuffle_h(inputtensor, shuffled_ratio, outliers_ID)
N = size(inputtensor,2);
outcome = inputtensor;
for i = 1:length(outliers_ID)
    j = outliers_ID(i);
    n = size(inputtensor,3);
    indexnum_shuffled = fix(n*shuffled_ratio);
    perm = 1:n;
    shuffled_index = randperm(n, indexnum_shuffled);
    perm(shuffled_index) = shuffled_index(randperm(indexnum_shuffled));
    outcome(:,j,:) = inputtensor(:,j,perm); 
end
end