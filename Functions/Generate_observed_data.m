function [X_observed, order_clean] = Generate_observed_data(M_gt, outlier_ratio, shuffled_ratio, noise_SNR)

    [D, V, w] = size(M_gt);
    XX = zeros(D,w*V);
    q = sort(randperm(V, fix(outlier_ratio*V)));
    num_shuffled = fix(D*shuffled_ratio);
    qq = zeros(1,w*length(q));
    
    for i = 1:w
        qq((i-1)*length(q)+1:i*length(q)) = q+(i-1)*V;
    end
    pp = 1:w*V;
    [~,pos] = ismember(qq,pp);
    pp(pos) = [];
    order_clean = pp;

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

    if (noise_SNR == Inf)
        X_observed = XX;
    else
        X_observed = awgn(XX, noise_SNR);
    end
end





