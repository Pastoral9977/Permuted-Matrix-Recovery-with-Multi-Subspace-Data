function [X_tilde] = permute_corruption(X,id_out,shuffled_ratio,seed)
rng(seed);
n = size(X,1);
X_tilde = X;
num_shuffled = fix(shuffled_ratio * n);
for ii = 1:length(id_out)
    j = id_out(ii);
    rp = randperm(n);
    shuffled_index = rp(1:num_shuffled);
    shuffled_how = randperm(num_shuffled);
    shuffled_pattern=shuffled_index(shuffled_how);
    perm = 1:1:n;
    perm(shuffled_index) = shuffled_pattern;
    X_tilde(:, j) = X(perm, j);
end
end