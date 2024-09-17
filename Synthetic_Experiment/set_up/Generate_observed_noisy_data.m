function [X_observed, outliers_ID, inliers_ID] = Generate_observed_noisy_data(M_gt, outlier_ratio, shuffled_ratio, noise_SNR)
    if nargin < 4
        noise_SNR = inf;
    end
    
    [D, V, num_groups] = size(M_gt);
    q = sort(randperm(V, fix(outlier_ratio*V)));
    num_shuffled = fix(D*shuffled_ratio);
    outliers_ID = zeros(1,num_groups*length(q));
    noise_level = 0.1^(noise_SNR/20);

    XX = reshape(M_gt, [D, V*num_groups]);

    noise = randn(D, num_groups*V);
    noise = noise_level * (noise ./ vecnorm(noise));

    for i = 1:num_groups
        outliers_ID((i-1)*length(q)+1:i*length(q)) = q+(i-1)*V;
    end
    
    inliers_ID = setdiff(1:num_groups*V, outliers_ID);
    for i =1:num_groups
         for j = 1:length(q)
            x = M_gt(:,q(j),i);
            perm = 1:D;
            shuffled_index = randperm(D, num_shuffled); % Randomly choose the shuffled indices
            perm(shuffled_index) = shuffled_index(randperm(num_shuffled)); % Randomly shuffle the choosed indices
            x = x(perm); 
            M_gt(:,q(j),i) = x;
         end
    end
    XX = reshape(M_gt, [D, V*num_groups]);

    if (noise_level == 0 )
        X_observed = XX;
    else
        X_observed = XX + noise;
    end
end


