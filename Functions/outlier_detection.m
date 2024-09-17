function [Inliers_id, Outliers_id] = outlier_detection(X_tilde, alpha, lambda, thres)

    if nargin < 4
        thres = 5*1e-4;
    end
    if nargin < 3
        lambda = 0.95;
    end
    if nargin < 2
        alpha = 5;
    end
    N = size(X_tilde,2);

    %% step 1: compute representation R from data
    gamma = @(X, y, lambda, alpha)  alpha*lambda/max(abs(X'*y));
    EN_solver =  @(X, y) rfss( full(X), full(y), lambda / gamma(X, y, lambda, alpha), (1-lambda) / gamma(X, y, lambda, alpha) );
    
    R = selfRepresentation(cnormalize(X_tilde), EN_solver);
    
    %% step 2: compute transition P from R (line 2 of Alg. 1)
    P = cnormalize(abs(R), 1)';
    
    %% step 3: compute \pi from P (line 3 - 7 of Alg. 1)
    T = 1000;
    pi = ones(1, N) / N;
    pi_bar = zeros(1, N);
    for ii = 1:T
        pi = pi * P;
        pi_bar = pi_bar + pi;
    end
    pi_bar = pi_bar / T;
    ids = 1:N;
    Inliers_id = ids(pi_bar>=thres);
    Outliers_id = ids(pi_bar<thres);
    
end