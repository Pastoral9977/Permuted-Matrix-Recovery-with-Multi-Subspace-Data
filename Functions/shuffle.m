function [outcome] = shuffle(inputvector, shuffled_ratio)
    n = length(inputvector);
    indexnum_shuffled = fix(n*shuffled_ratio);
    perm = 1:n;
    shuffled_index = randperm(n, indexnum_shuffled);
    perm(shuffled_index) = shuffled_index(randperm(indexnum_shuffled));
    outcome = inputvector(perm); 
end