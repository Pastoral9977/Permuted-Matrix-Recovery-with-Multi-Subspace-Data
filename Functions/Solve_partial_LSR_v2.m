function [X_solved, outlier_class_solved] = Solve_partial_LSR_v2(X_tilde, idOutliers, inlier_class_solved , U_solved)
    %--------------------------------------------------------------------------
    % This is the function to reconstructed the original matrix based on the
    % estimated bases U_solved and estimated outliers id idOutliers.
    % 
    %--------------------------------------------------------------------------
    % Copyright @ XLQ
    %--------------------------------------------------------------------------

    tic
    X_solved = X_tilde;
    [a,b,c] = size(U_solved); % b = rrank, c = num_groups
    b = fix(b/2); % for speeding up!
    U = reshape(U_solved(:,1:b,:), [a, b*c]);
    N = size(X_tilde,2);
    n = length(idOutliers);
    idInliers = inlier_class_solved(:,1);
    inlier_class = inlier_class_solved(:,2);
    [U2,~,~] = svd(X_tilde(:,idInliers));
    U = U2(:,1:(b*c));
    display_progress()
    
%   step1: preliminarily recovery
    for i = 1:n
        j = idOutliers(i);
        y = X_tilde(:, j);
        [c_hat] = LSR_v2(U, y, 50);
        X_tilde(:, j) = U * c_hat;
        update_progress(i, n);
    end

%   step2: calculate transition matrix
    lambda = 0.95;
    alpha = 5;
    gamma = @(X, y, lambda, alpha)  alpha*lambda/max(abs(X'*y));
    EN_solver =  @(X, y) rfss( full(X), full(y), lambda / gamma(X, y, lambda, alpha), (1-lambda) / gamma(X, y, lambda, alpha) );
    R = selfRepresentation(cnormalize(X_tilde), EN_solver);
    P = cnormalize(abs(R), 1)';
        
%   step3: initial probability
    row_indices = 1:n;
    col_indices = idOutliers;
    values = ones(1, n)*(n-1);
    f = sparse(row_indices, col_indices, values, n, N);
    pi = ones(n, N)+f;
    pi(:,idInliers) = 0;
    row_norms = vecnorm(pi, 1, 2);
    pi = pi./row_norms;
    pi_bar = sparse(n, N);
    
%   step4: transition process
    T = 100;
    for ii = 1:T
        pi = pi * P;
        pi_bar = pi_bar + pi;
    end
    pi_bar = pi_bar / T;
    
%   step5: final process
    scoresMat = zeros(n, c);
    for ii = 1:c
        idx = idInliers(inlier_class==ii);
        score = mean(pi_bar(:,idx),2);
        scoresMat(:,ii) = score;
    end
    
    [~,outlier_class_solved] = max(scoresMat, [], 2);
    outlier_class_solved = reshape(outlier_class_solved, [n,1]);
    
end





