
% ��֤�����ɶȵĿ����ֲ�����̬�ֲ��Ĳ��
clear;
close all;

% �������ɶ�
k_values = [2, 5, 10, 50, 100, 500, 1000];

% �����µ� figure
figure;

% ����PDF
for i = 1:length(k_values)
    k = k_values(i);
    
    % ���ܶ��� x ��Χ
    x_range = 0:0.1:(k + 5 * sqrt(2 * k));
    
    % ���㿨���ֲ�PDF
    chi2_pdf = chi2pdf(x_range, k);
    
    % ������Ӧ����̬�ֲ�����
    mu = k;
    sigma = sqrt(2 * k);
    
    % ������̬�ֲ�PDF
    normal_pdf = normpdf(x_range, mu, sigma);
    
    % ����PDF
    subplot(length(k_values), 1, i);
    hold on;
    plot(x_range, chi2_pdf, 'b-', 'DisplayName', '\chi^2(k) PDF');
    plot(x_range, normal_pdf, 'r--', 'DisplayName', 'N(k, 2k) PDF');
    
    % ���㲢��ע3-sigmaλ��
    upper_3sigma = mu + 3 * sigma;
    y_max = max([chi2_pdf, normal_pdf]); % ��� y ֵ�����ڻ��������߶�
    plot([upper_3sigma, upper_3sigma], [0, y_max*1.1], 'p--', 'DisplayName', '3-\sigma Bound');
    
    % ���㿨���ֲ���0.01ˮƽ��������
    chi2_001 = chi2inv(0.99, k); % 0.01 ˮƽ��Ӧ�ĸ����� 0.99
    
    % ��ע�����ֲ���0.01ˮƽ��������
    plot([chi2_001, chi2_001], [0, y_max*1.1], '-.', 'color', [0, 160/255, 0], 'DisplayName', '\chi^2_{0.99}(k) Confidence Bound');
    
    title(sprintf('PDF Comparison for degrees of freedom k = %d', k));
    lgd = legend('show', 'Location', 'bestoutside');
    lgd.FontSize = 8.5;
    hold off;
end

% ����ͼ�񲼾�
sgtitle('Comparison of Chi-Squared and Normal Distributions for Different Degrees of Freedom');



