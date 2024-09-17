% Algorithm 2 in the submitted paper "Unlabeled Principal Component Analysis"
function [x_hat, res] = LSR_v3(A,y, retain_ratio, num_iterations)
    [m,r] = size(A);

    if nargin < 4
        num_iterations = int64(m-r);
    end
    
    if nargin < 3
        retain_ratio = 0.7;
    end
    
    num_iterations = min(num_iterations, int64(m-r));
    rm = r+int64(m-r)*retain_ratio;
    
    diffA = decelerate_to_end(m, rm, num_iterations);
    
    for i=1:num_iterations
        
        k = diffA(i);
        
        x_hat = A \ y;
        [~,idx] = sort(abs(y-A*x_hat), 'descend');
        idx = idx(1:min(k, m - r));
        
        y(idx) = [];
        A(idx, :) = [];
        
        [m, ~] = size(A);
        if m<=r
            break;
        end
        
    end
    x_hat = A \ y;
        
    res = struct();
    y_hat = A*x_hat;
    res.norm_diff = norm(y-y_hat);
    res.vec_cos_dist = 1-calculate_cosine_similarity(y, y_hat);
    res.subspace_cos_dist =1-cosine_subspace_angle(y, A);

    
end

function result = decelerate_to_end(n, m, k)
    result = zeros(1, k + 1);
    result(1) = n;
    
    decay_factor = logspace(0, 1, k + 1) - 1;
    decay_factor = decay_factor(end:-1:1);
    decay_factor = decay_factor / max(decay_factor);
    
    for i = 2:k+1
        result(i) = m + int64(decay_factor(i) * (n - m));
    end
    
    result = -diff(result);
    
end



function cosine_similarity = calculate_cosine_similarity(vector1, vector2)
    norm_vector1 = norm(vector1);
    norm_vector2 = norm(vector2);
    
    cosine_similarity = dot(vector1, vector2) / (norm_vector1 * norm_vector2+eps);
end

function cos_theta = cosine_subspace_angle(v, U)
    [Q, ~] = qr(U, 0);
    v_proj = Q * (Q' * v);
    cos_theta = dot(v, v_proj) / (norm(v) * norm(v_proj)+eps);
end