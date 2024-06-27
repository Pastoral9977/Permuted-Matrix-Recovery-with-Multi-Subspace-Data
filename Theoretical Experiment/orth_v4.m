close all
n = 1000;
y = randn(n,1);
y = y/norm(y,2);

shuffled_ratios = 0:0.05:1;
m = length(shuffled_ratios);
% m = 20;
% shuffled_ratios = 1/m:1/m:1;
num_samples = 100; % 每个sr取样次数
res_mean = zeros(1,m);
res_std = zeros(1,m);
res_mean2 = zeros(1,m);
res_std2 = zeros(1,m);
for i = 1:m
    sr_results = zeros(1, num_samples);
    sr_results2 = zeros(1, num_samples);
    for j = 1:num_samples
        z = y;
        sr = shuffled_ratios(i);
        k = fix(sr*n);
        last_elements = y(end-k+1:end);
        z(end-k+1:end) = last_elements(randperm(k));
        sr_results(j) = y'*z;
        if i == 1
            continue
        end
        sr_results2(j) = last_elements'*z(end-k+1:end);
    end
    res_mean(i) = mean(sr_results);
    res_std(i) = std(sr_results);
    res_mean2(i) = mean(sr_results2);
    res_std2(i) = std(sr_results2);
end

figure;
x = shuffled_ratios;
subplot(1,2,1)
    y1 = res_mean + res_std;
    y2 = res_mean - res_std;

    % 绘制灰度区域
    fill([x fliplr(x)], [y1 fliplr(y2)], [0.8 0.8 0.8], 'EdgeColor', 'none');
    hold on;

    % 绘制平均值线
    plot(shuffled_ratios, res_mean, 'b', 'LineWidth', 1.5);

    % 绘制 y = 1 - x 虚线
    y = 1-x;
    plot(x, y, 'LineStyle', '--', 'LineWidth', 1.5, 'color', 'red');

    ylim([0,1]);
    grid on;
    xlabel('Shuffled Ratios');
    ylabel('Average Dot Product');
    title("y'*\tilde(y) for Different Shuffled Ratios");

    % 添加图例
    legend('Standard Deviation', 'Average Dot Product', 'y = 1 - x', 'Location', 'Best');
    
subplot(1,2,2)
    y3 = res_mean2 + res_std2;
    y4 = res_mean2 - res_std2;
    
    % 绘制灰度区域
    fill([x fliplr(x)], [y3 fliplr(y4)], [0.8 0.8 0.8], 'EdgeColor', 'none');
    hold on;

    % 绘制平均值线
    plot(shuffled_ratios, res_mean2, 'b', 'LineWidth', 1.5);
    ylim([-0.5, 0.5])
    grid on
    xlabel('Shuffled Ratios');
    ylabel('Average Dot Product');
    title('Average Dot Product for Fully Shuffled Ratios with Different Size in a Norm 1 vector');

    % 添加图例
    legend('Standard Deviation', 'Average Dot Product', 'y = 1 - x', 'Location', 'Best');
    
