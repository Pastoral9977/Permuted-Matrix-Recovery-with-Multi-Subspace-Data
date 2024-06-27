
% 考察列正交阵Q与正态阵U的行之一：模长
% 更精准请看orth_v10

% f1()
% f2()
f3()


function f1()
    close all
    n_values = [1000, 500, 200];  % 不同的n值
    m = 20;

    figure;
    hold on;

    colors = {'b', 'g', 'r'};  % 用于不同n值的颜色
    legends = {};

    for idx = 1:length(n_values)
        n = n_values(idx);
        rs = n/m:n/m:n;

        A = randn(n)/sqrt(n);
        res_mean = zeros(1, m);
        res_std = zeros(1, m);
        res_mean_normal = zeros(1, m);
        res_std_normal = zeros(1, m);

        for i = 1:m
            r = rs(i);
            U = A(:, 1:r); 
            [Q, ~] = qr(U, 0);
            norms = vecnorm(Q, 2, 2);
            norms_normal = vecnorm(U, 2, 2);
            res_mean(i) = mean(norms);
            res_std(i) = std(norms);
            res_mean_normal(i) = mean(norms_normal);
            res_std_normal(i) = std(norms_normal);
        end

        x = rs;
        y1 = res_mean + res_std;
        y2 = res_mean - res_std;
        y3 = res_mean_normal + res_std_normal;
        y4 = res_mean_normal - res_std_normal;

        % 绘制灰度区域
        fill([x fliplr(x)], [y1 fliplr(y2)], colors{idx}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        fill([x fliplr(x)], [y3 fliplr(y4)], colors{idx}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

        % 绘制平均值线
        plot(rs, res_mean, colors{idx}, 'LineWidth', 1.5);
        plot(rs, res_mean_normal, colors{idx}, 'LineStyle', '--', 'LineWidth', 1.5);

        % 添加图例条目
        legends{end+1} = ['Mean (Q), n = ' num2str(n)];
        legends{end+1} = ['Mean (U), n = ' num2str(n)];
    end

    % 绘制参考线 y = sqrt(r/n) 对于 n = 1000 的情况
    r = linspace(0, max(n_values), 100);
    for idx = 1:length(n_values)
        n = n_values(idx);
        rr = r(r<=n);
        plot(rr, sqrt(rr/n), 'color', colors{idx}, 'LineStyle','--', 'LineWidth', 1.5);
        legends{end+1} = ['y = sqrt(r/n), n = ' num2str(n)];
    end

    % 添加图例
    legend(legends, 'Location', 'Best');

    % 设置图表属性
    ylim([0, 1.1]);
    grid on;
    xlabel('r');
    ylabel('Values');
    title('Comparison of Norms for Different n Values with Shaded Standard Deviation Areas');
end


function f2()
    close all
    n = 500;
    m = 20;
    rs = n/m:n/m:n;
    
    A = randn(n)/sqrt(n);
    res_mean = zeros(1, m);
    res_std = zeros(1, m);
    res_mean_normal = zeros(1, m);
    res_std_normal = zeros(1, m);
    
    display_progress()
    for i = 1:m
        r = rs(i);
        U = A(:, 1:r); 
        [Q, ~] = qr(U, 0);
        norms = vecnorm(Q, 2, 2);
        norms_normal = vecnorm(U, 2, 2);
        res_mean(i) = mean(norms);
        res_std(i) = std(norms);
        res_mean_normal(i) = mean(norms_normal);
        res_std_normal(i) = std(norms_normal);
        update_progress(i, m)
    end
    
    x = rs;
    y1 = res_mean + res_std;
    y2 = res_mean - res_std;
    y3 = res_mean_normal + res_std_normal;
    y4 = res_mean_normal - res_std_normal;
    
    % 绘制Q的子图
    figure;
    subplot(3, 1, 1);
    fill([x fliplr(x)], [y1 fliplr(y2)], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    hold on;
    plot(rs, res_mean, 'b', 'LineWidth', 1.5);
    plot(rs, sqrt(rs/n), 'r--', 'LineWidth', 1.5);
    ylim([0, 1.1]);
    grid on;
    xlabel('r');
    ylabel('Values');
    title('Norms of Q with Shaded Standard Deviation Area');
    legend('Std Dev (Q)', 'Mean (Q)', 'y = sqrt(r/n)', 'Location', 'Best');

    % 绘制U的子图
    subplot(3, 1, 2);
    fill([x fliplr(x)], [y3 fliplr(y4)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    hold on;
    plot(rs, res_mean_normal, 'g', 'LineWidth', 1.5);
    plot(rs, sqrt(rs/n), 'r--', 'LineWidth', 1.5);
    ylim([0, 1.1]);
    grid on;
    xlabel('r');
    ylabel('Values');
    title('Norms of U with Shaded Standard Deviation Area');
    legend('Std Dev (U)', 'Mean (U)', 'y = sqrt(r/n)', 'Location', 'Best');
    
    % 绘制U的子图
    subplot(3, 1, 3);
    plot(rs, (res_mean-res_mean_normal)./res_mean)
    grid on;
    ylim([-0.01, 0.01])
    xlabel('r');
    ylabel('Values');
    legend('(Mean (U)-Mean (Q))./Mean (Q)', 'Location', 'Best');
end

function f3()
    
%     rng(42)
    close all
    combinations = [200, 10; 500, 30; 1000, 100]; % r<<n
%     combinations = [200, 100; 500, 250; 1000, 500]; % r=n/2
%     combinations = [200, 150; 500, 450; 1000, 990]; % r接近n
    num_combinations = size(combinations, 1);

    figure;
    for i = 1:num_combinations
        n = combinations(i, 1);
        r = combinations(i, 2);

        A = randn(n)/sqrt(n);

        U = A(:, 1:r); 
%         [Q, ~] = qr(U, 0);
        [Q, ~, P] = svd(U);
        Q = Q(:, 1:r)*P';
        norms = vecnorm(Q, 2, 2);
        norms_normal = vecnorm(U, 2, 2);

        x = 1:n;

        subplot(num_combinations, 1, i);
        yyaxis left
        h1 = plot(x, norms, 'b', 'LineWidth', 1.5);
        hold on;
        h2 = plot(x, norms_normal, 'g', 'LineWidth', 1.5);
        h3 = line([1, n], [sqrt(r/n), sqrt(r/n)], 'LineStyle' , '--', 'color', 'r', 'LineWidth', 1);
        ylim([0, 1.1]);
        xlabel('n');
        ylabel('Norms');
        title(sprintf('Norms of Matrice Rows and Difference (n=%d, r=%d)', n, r));
        grid on;

        yyaxis right
        h4 = plot(x, (norms - norms_normal)./norms, 'LineWidth', 0.5, 'Marker','.');

        % Create filled area for -0.05 to 0.05 range
        hold on;
        h5 = fill([x, fliplr(x)], [0.05*ones(1, n), -0.05*ones(1, n)], [0.4, 0.5, 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        h6 = fill([x, fliplr(x)], [0.1*ones(1, n), -0.1*ones(1, n)], [0.8, 0.7, 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

        ylabel('diff ratio');
        ylim([-0.2, 0.15]*1.8);
        grid on;

        % Create two separate legends for each subplot
        legend([h1, h2, h3], {'Orthogonal Q', 'Normal U', 'y = sqrt(r/n)'}, 'Location', 'NorthWest');
        
        % Create a second invisible axes for the second legend
        ah = axes('Position', get(gca, 'Position'), 'Color', 'none', 'XTick', [], 'YTick', [], 'Box', 'off');
        ah.Visible = 'off';
        legend(ah, [h4, h5, h6], {'Difference ratio','Area Under 5%','Area Under 10%'}, 'Location', 'NorthEast');
    end
end






