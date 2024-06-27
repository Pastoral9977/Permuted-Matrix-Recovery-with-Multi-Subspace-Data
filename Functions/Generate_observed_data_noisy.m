function [X_observed, X_gt_1, order_inliers] = Generate_observed_data_noisy(M_gt, outlier_ratio, shuffled_ratio, noise_SNR)
    [D, V, w] = size(M_gt);
    XX = zeros(D,w*V);
    q = sort(randperm(V, fix(outlier_ratio*V)));
    num_shuffled = fix(D*shuffled_ratio);
    order_outliers = zeros(1,w*length(q));
    noise_level = 0.1^(noise_SNR/20);

    for i = 1:w
        XX(:,(i-1)*V+1:i*V) = M_gt(:,:,i);
    end

    if (noise_level == 0 )
        X_gt_1 = XX;
    else
        noise = randn(D, w*V);
        noise = noise_level * (noise ./ vecnorm(noise));
        X_gt_1 = XX + noise;
    end

    for i = 1:w
        order_outliers((i-1)*length(q)+1:i*length(q)) = q+(i-1)*V;
    end

    order_inliers = setdiff(1:w*V, order_outliers);

    for i =1:w
         for j = 1:length(q)
            x = M_gt(:,q(j),i);
            perm = 1:D;
            shuffled_index = randperm(D, num_shuffled);
            perm(shuffled_index) = shuffled_index(randperm(num_shuffled));
            x = x(perm); 
            M_gt(:,q(j),i) = x;
        end
        XX(:,(i-1)*V+1:i*V) = M_gt(:,:,i);
    end

    if (noise_level == 0 )
        X_observed = XX;
    else
        X_observed = XX + noise;
    end
end


