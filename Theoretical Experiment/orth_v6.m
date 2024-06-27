close all
n = 1000;
r = 20;
ratios = [0.1, 0.2, 0.5, 0.8, 1.0]; % ���Ը�����Ҫ����
num_trials = 100; % ÿ��ratio�µ�ʵ�����

% ��ʼ�����ڱ������ľ���
dotps_all = cell(length(ratios), 1);
dotps_all_np = cell(length(ratios), 1);

for i = 1:length(ratios)
    ratio = ratios(i);
    m = fix(n * ratio);
    dotps_trials = zeros(r, num_trials);
    dotps_trials_np = zeros(r, num_trials);
    
    for trial = 1:num_trials
        % �����������A��������QR�ֽ�
        A = randn(n, r) / sqrt(n);
        [U, ~] = qr(A, 0);
        
        % �������beta��������y
        beta = randn(r, 1);
        beta = beta / norm(beta);
        y = U * beta;
        
        % ��ȡ���m��Ԫ��
        y2 = y(end-m+1:end);
        U2 = U(end-m+1:end, :);
        
        % ��������û�����P
        P = eye(m);
        perms = randperm(m);
        P = P(perms, :);
        
        % �����ڻ�
        dotps = U2' * P * y2;
        dotps_np = U2'* y2;
        
        % ����ÿ��ʵ��Ľ��
        dotps_trials(:, trial) = dotps;
        dotps_trials_np(:, trial) = dotps_np;
    end
    
    dotps_all{i} = dotps_trials;
    dotps_all_np{i} = dotps_trials_np;
end

% ���ƽ��ͼ
figure;
colors = lines(length(ratios)); % Ԥ����ɫ
colors_grey = linspace(0.1, 0.9, length(ratios)); % �Ҷ�ֵ��Χ

subplot(1,2,1)
hold on
legend_entries = cell(2*length(ratios), 1);
for i = 1:length(ratios)
    ratio = ratios(i);
    m = fix(n * ratio);
    dotps_mean = mean(dotps_all{i}, 2);
    dotps_std = std(dotps_all{i}, 0, 2);
    
    % ����ƽ��ֵ�ͷ�������
    fill([1:r, r:-1:1], [dotps_mean - dotps_std; flipud(dotps_mean + dotps_std)], ...
        [colors_grey(i), colors_grey(i), colors_grey(i)], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(1:r, dotps_mean, 'Color', colors(i, :), 'LineWidth', 2);
    legend_entries{2*i-1} = sprintf('Ratio = %.1f (Standard Error)', ratio);
    legend_entries{2*i} = sprintf('Ratio = %.1f (Mean)', ratio);
end
grid on
xlabel('Dimention Index');
ylabel('Dot Product Value');
ylim([-0.5,0.5]);
title("Dot Products U2' * P * y2 for Different Ratios");
legend(legend_entries, 'Location', 'Best');
hold off

subplot(1,2,2)
hold on
for i = 1:length(ratios)
    ratio = ratios(i);
    m = fix(n * ratio);
    dotps_mean_np = mean(dotps_all_np{i}, 2);
    dotps_std_np = std(dotps_all_np{i}, 0, 2);
    
    % ����ƽ��ֵ�ͷ�������
    fill([1:r, r:-1:1], [dotps_mean_np - dotps_std_np; flipud(dotps_mean_np + dotps_std_np)], ...
        [colors_grey(i), colors_grey(i), colors_grey(i)], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(1:r, dotps_mean_np, 'Color', colors(i, :), 'LineWidth', 2);
    legend_entries{2*i-1} = sprintf('Ratio = %.1f (Standard Error)', ratio);
    legend_entries{2*i} = sprintf('Ratio = %.1f (Mean)', ratio);
end
grid on
xlabel('Dimention Index');
ylabel('Dot Product Value');
ylim([-0.5,0.5]);
title("Dot Products U2' * y2 for Different Ratios without Permutation");
% legend(arrayfun(@(x) sprintf('Ratio = %.1f', x), ratios, 'UniformOutput', false));
legend(legend_entries, 'Location', 'Best');
hold off







% % ����ͼ��
% saveas(gcf, 'dot_products_variance.png');


% close all
% n = 1000;
% r = 20;
% ratio = 1;
% m = fix(n*ratio);
% 
% A = randn(n, r)/sqrt(n);
% [U, ~] = qr(A, 0);
% beta = randn(r, 1);
% beta = beta/norm(beta);
% y = U*beta;
% 
% y2 = y(end-m+1
% );
% U2 = U(end-m+1
% ,:);
% 
% P = eye(m);
% perms = randperm(m);
% P = P(perms, :);
% 
% dotps = U2'Py2;
% 
% figure;
% plot(dotps)
% grid on
