function [alpha, beta] = beta_parameters(mu, variance)
    % Check for valid mean and variance values
    if mu <= 0 || mu >= 1
        error('Mean must be between 0 and 1 (exclusive).');
    end
    
    if variance <= 0
        error('Variance must be positive.');
    end

    % Compute alpha and beta
    alpha = ((1 - mu) * mu^2 / variance) - mu;
    beta = alpha * (1 - mu) / mu;
end