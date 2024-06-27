% Algorithm 2 in the submitted paper "Unlabeled Principal Component Analysis"
function [x_hat, res] = LSR_v2(A,y,num_iterations)
    [m,r] = size(A);
    
    if nargin < 3
        num_iterations = int64(m-r);
    
        for i=1:num_iterations
            x_hat = A \ y;
            [~,idx] = max(abs(y-A*x_hat));
            y(idx) = [];
            A(idx, :) = [];
        end
        
    else
        for i = 1:num_iterations
            x_hat = A \ y;
            [~, idx] = sort(abs(y - A * x_hat), 'descend');

%             k = max(int64((m - r) / (0.5*num_iterations+i)), 1);
            k = int64((m - r) / (0.5*num_iterations+i));

            idx = idx(1:min(k, m - r));
            y(idx) = [];
            A(idx, :) = [];
            [m, ~] = size(A);

            if m<=r
                break;
            end
        end
        
        x_hat = A \ y;
        
    end
    res = struct();
    y_hat = A*x_hat;
    res.norm_diff = norm(y-y_hat);
    res.vec_cos_dist = 1-calculate_cosine_similarity(y, y_hat);
    res.subspace_cos_dist =1-cosine_subspace_angle(y, A);

    
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