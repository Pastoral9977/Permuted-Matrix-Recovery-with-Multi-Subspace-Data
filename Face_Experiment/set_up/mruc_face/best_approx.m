function [lambda,mat]=best_approx(B,ranks)
    [U,S,V] = svd(B);
    mat=U(:,1:ranks)*S(1:ranks,1:ranks)*V(:,1:ranks)';
    lambda=S(ranks,ranks);  % change lambda with respect to ranks
end