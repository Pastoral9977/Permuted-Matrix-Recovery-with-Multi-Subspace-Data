% This is the function for solving the SLR of unselected points.
% INPUT
% X_tilde: observed data matrix.
% idOutliers: indices of unselected points.
% U_solved: 3-d tensor of solved bases of subspaces


function [X_solved, outlier_class_solved] = Solve_all(X_tilde, idOutliers, U_solved)

X_solved = X_tilde;
X_hat = zeros(size(X_tilde,1),size(U_solved,3));
es = zeros(1,size(U_solved,3)); 
outlier_class_solved = zeros(1, length(idOutliers));
[a,b,c] = size(U_solved);
outer_rank = min(4,(c-1)*b);
UU_solved = zeros(a,outer_rank,c);

for i = 1:c
    M = zeros(a, b*(c-1));
    for k = 1:(c-1)
        kk = setdiff(1:c,i);
        ink = kk(k);
        M(:,(k-1)*b+1:k*b) = U_solved(:,1:b,ink);
    end
    UU_solved(:,:,i) = M;
end

for i = 1:length(idOutliers)
    j = idOutliers(i);
    y = X_tilde(:, j);
    for k = 1:size(U_solved,3)
        percent = 1;
        [c_hat] = AIEM(U_solved(:,:,k), y);
        X_hat(:,k) = U_solved(:,:,k) * c_hat;
        v1 = retain_smallest_part(sort(y)/norm(y) - sort(X_hat(:,k))/norm(X_hat(:,k)), percent);
        a = vecnorm((v1).^1,1);
%         u1 = retain_smallest_part(X_hat(:,k), 1, y);
%         cos_sim1 = calculate_cosine_similarity(y, u1);
%         b = abs(norm(y,1)-norm(X_hat(:,k),1))^0.5;
        
        
        es(k) = a;
    end
    [~,p] = min(es);
    X_solved(:, j) = X_hat(:,p);
    outlier_class_solved(i) = p;
end

end
            
function cosine_similarity = calculate_cosine_similarity(vector1, vector2)
    % ����������ģ
    norm_vector1 = norm(vector1)+eps;
    norm_vector2 = norm(vector2)+eps;
    
    % �����������ƶ�
    cosine_similarity = dot(vector1, vector2) / (norm_vector1 * norm_vector2)+eps;
    
    % �������Ҿ��루�������ƶȵĲ�����
%     cosine_distance = 1 - cosine_similarity;
end

function [a_modified, retain_indices] = retain_smallest_part(a, percent, b)
    if nargin < 3
        b = zeros(size(a));
    end

    % ��ȡ�����ĳ���
    n = length(a);
    
    % ������Ҫ������Ԫ������
    num_to_retain = ceil(percent * n);
    
    % �������������򣬲���ȡ����������
    [~, sorted_indices] = sort(a);
    
    % ������С��percent������Ԫ�ص�����
    retain_indices = sorted_indices(1:num_to_retain);
    
    % ��ʼ��������������Ԫ����Ϊ0
    a_modified = b;
    
    % ���������б�����Сpercent������Ԫ��
    a_modified(retain_indices) = a(retain_indices);
end
        


















