
% 验证大自由度的卡方分布与正态分布的差别
clear;
close all;

% 设置自由度
k_values = [2, 5, 10, 50, 100, 500, 1000];

% 创建新的 figure
figure;

% 绘制PDF
for i = 1:length(k_values)
    k = k_values(i);
    
    % 智能定义 x 范围
    x_range = 0:0.1:(k + 5 * sqrt(2 * k));
    
    % 计算卡方分布PDF
    chi2_pdf = chi2pdf(x_range, k);
    
    % 计算相应的正态分布参数
    mu = k;
    sigma = sqrt(2 * k);
    
    % 计算正态分布PDF
    normal_pdf = normpdf(x_range, mu, sigma);
    
    % 绘制PDF
    subplot(length(k_values), 1, i);
    hold on;
    plot(x_range, chi2_pdf, 'b-', 'DisplayName', '\chi^2(k) PDF');
    plot(x_range, normal_pdf, 'r--', 'DisplayName', 'N(k, 2k) PDF');
    
    % 计算并标注3-sigma位置
    upper_3sigma = mu + 3 * sigma;
    y_max = max([chi2_pdf, normal_pdf]); % 最大 y 值，用于绘制线条高度
    plot([upper_3sigma, upper_3sigma], [0, y_max*1.1], 'p--', 'DisplayName', '3-\sigma Bound');
    
    % 计算卡方分布的0.01水平置信上限
    chi2_001 = chi2inv(0.99, k); % 0.01 水平对应的概率是 0.99
    
    % 标注卡方分布的0.01水平置信上限
    plot([chi2_001, chi2_001], [0, y_max*1.1], '-.', 'color', [0, 160/255, 0], 'DisplayName', '\chi^2_{0.99}(k) Confidence Bound');
    
    title(sprintf('PDF Comparison for degrees of freedom k = %d', k));
    lgd = legend('show', 'Location', 'bestoutside');
    lgd.FontSize = 8.5;
    hold off;
end

% 调整图像布局
sgtitle('Comparison of Chi-Squared and Normal Distributions for Different Degrees of Freedom');



