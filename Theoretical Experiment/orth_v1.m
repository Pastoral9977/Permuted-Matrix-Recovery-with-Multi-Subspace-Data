close all
clear

n = 500;
r = 20;
sigma = 1/sqrt(n);
B = randn(n)*sigma;
B = B(:, 1:r);

[Q, ~] = qr(B, 0);

W = B*B';
G = Q*Q';
subplot(1,3,1)
imagesc(W)
colorbar
colormap('hot')
title('Normal Matrices Dot Product W')
subplot(1,3,2)
imagesc(G)
colorbar
colormap('hot')
title('Orthononal Matrices Dot Product G')
subplot(1,3,3)
T = (W-G)./abs(G);
% T = (W-G);
imagesc(T)
colorbar
colormap('hot')
caxis([-0.2,0.2]);
title('Differences (W-G)./abs(G)')

% 发现用W模拟G不太合适，还是得用orth_v7的正态参数模拟具体每个entry



































% hFig = figure;
% subplot(1,3,1)
% imagesc(W);
% colorbar;
% colormap('hot');
% % caxis([0, 0.12]);
% title('Heatmap of Matrix Generated From Normal Distribution');
% xlabel('Column Index');
% ylabel('Row Index');
% 
% subplot(1,3,2)
% imagesc(G);
% colorbar;
% colormap('hot');
% % caxis([0, 0.12]);
% title('Heatmap of Matrix Generated From Orthonogal Basis');
% xlabel('Column Index');
% ylabel('Row Index');
% 
% subplot(1,3,3)
% imagesc(abs(G-W));
% colorbar;
% colormap('hot');
% caxis([0, 1]);
% title('Heatmap of abs(G-W)');
% xlabel('Column Index');
% ylabel('Row Index');
% 
% % 设置图形窗口的位置和大小
% % 格式：[left, bottom, width, height]
% left = 680; % 左边缘到屏幕左边缘的距离
% bottom = 120; % 下边缘到屏幕下边缘的距离
% width = 800; % 窗口的宽度
% height = 300; % 窗口的高度
% 
% % 设置图形窗口的位置和大小
% set(hFig, 'Position', [left, bottom, width, height]);
% 
% 
% disp(['Normal Distribution: ', num2str(diag_non_diag_ratio(W))]);
% disp(['Orthonogal Basis: ', num2str(diag_non_diag_ratio(G))]);
% 
% 
% mask = ~eye(size(W));
% non_diag_elements = W(mask);
% [max_abs_value, linear_idx] = max(abs(non_diag_elements));
% median_abs_value = median(abs(non_diag_elements));
% 
% [row, col] = find(mask);
% row = row(linear_idx);
% col = col(linear_idx);
% 
% disp('Maximum absolute value of non-diagonal elements:');
% disp(max_abs_value);
% disp('Position of the maximum absolute value:');
% disp(['Row: ', num2str(row), ', Column: ', num2str(col)]);
% disp('Median absolute value of non-diagonal elements:');
% disp(median_abs_value);









function ratio = diag_non_diag_ratio(matrix)
    % 获取矩阵的大小
    [rows, cols] = size(matrix);

    % 确保输入是一个方阵
    if rows ~= cols
        error('输入必须是一个方阵');
    end

    % 获取对角元素的绝对值
    diag_elements = abs(diag(matrix));

    % 计算对角元素的最小绝对值
    min_diag = min(diag_elements);

    % 创建非对角元素掩码
    non_diag_mask = ~eye(rows);

    % 获取非对角元素的绝对值
    non_diag_elements = abs(matrix(non_diag_mask));

    max_non_diag = max(non_diag_elements);

    % 计算比值
    ratio = min_diag / max_non_diag;
end



% [norms1, norms2] = orth_v1(n,r,d);
% norm(norms1-norms2, 'inf')

function [norms1, norms2] = orth_V1(n, r, d)

    [U,~,~] = svd(randn(n));
    U = U(:,1:r);
    u2 = U(d+1:end,:);
    
    norms1 = vecnorm(U,2,2);
    figure;
    subplot(1,3,1)
    plot(norms1)
    ylim([0,1])
    grid on

    A = U*u2';
    
    norms2 = vecnorm(A,2,2);
    subplot(1,3,2)
    plot(norms2)
    ylim([0,1])
    grid on
    
    norms3 = norms1-norms2;
    subplot(1,3,3)
    plot(norms3)
    ylim([0,1])
    grid on

end















