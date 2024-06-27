% 考察列正交阵Q与正态阵U的行之二：不同行的内积

clear
close all

% 不同的 n 和 r 组合
n_values = [100, 500, 1000];
rn_ratios = [1/0.9, 2, 10, 25];

% 设置子图布局
num_rows = length(rn_ratios);
num_cols = length(n_values);
rows = num_rows;
cols = num_cols;
num_plots = num_cols * num_rows;

% Margin settings
left_margin = 0.08;
right_margin = 0.08;
top_margin = 0.12; % Increased top margin to make space for the title
bottom_margin = 0.20; % Increased bottom margin to make space for the legend
h_spacing = 0.05;
v_spacing = 0.05;

% Calculate subplot width and height
subplot_width = (1 - left_margin - right_margin - (num_cols - 1) * h_spacing) / num_cols;
subplot_height = (1 - top_margin - bottom_margin - (num_rows - 1) * v_spacing) / num_rows;

figure('Position', [280, 80, 1000, 600]); % Adjust the figure size

% 存储每个 subplot 的 handle
axes_handles = gobjects(rows, cols);

% 循环不同的 n 和 r 组合进行实验
plot_index = 0;
for r_idx = 1:rows
    for n_idx = 1:cols
        n = n_values(n_idx);
        r = fix(n / rn_ratios(end - r_idx + 1));

        plot_index = plot_index + 1;

        fprintf('Running experiment for n = %d, r = %d\n', n, r);

        A = randn(n, r) / sqrt(n);
        [U, ~] = qr(A, 0);
        R = U * U';
        dotps = R(triu(true(size(R)), 1));

        edges = linspace(-0.15, 0.15, 150 + 1);

        % Calculate position for the subplot
        subplot_left = left_margin + (n_idx - 1) * (subplot_width + h_spacing);
        subplot_bottom = 1 - top_margin - r_idx * (subplot_height + v_spacing) + v_spacing;

        % Create subplot
        subplot_handle = subplot('Position', [subplot_left, subplot_bottom, subplot_width, subplot_height]);
        hold on

        h1 = histogram(dotps, 'BinEdges', edges, 'Normalization', 'probability');
        ylim([0, max(h1.Values) * 1.6]);
        
        if r_idx == rows
            xlabel('$u_i^\top u_j$', 'FontSize', 16, 'Interpreter', 'latex');
        end
        if n_idx == 1
            ylabel('Probability', 'FontSize', 13);
        end
        title(['$n = ' num2str(n) ', \ r = ' num2str(r) '$'], 'Interpreter', 'latex');
        grid on;

        % 计算正态分布概率密度曲线
        mean_value = 0;
        std_value_1 = sqrt(r / n^2);
        std_value_2 = sqrt(r * (n - r) / (n * (n - 1) * (n + 2)));
        x_values = linspace(min(edges), max(edges), 1000);
        y_values_1 = normpdf(x_values, mean_value, std_value_1);
        y_values_2 = normpdf(x_values, mean_value, std_value_2);

        hold on;
        h2 = plot(x_values, y_values_1 * (edges(2) - edges(1)), 'g-', 'LineWidth', 1.2);
        h3 = plot(x_values, y_values_2 * (edges(2) - edges(1)), 'r-', 'LineWidth', 1.8);

        hold off;

        summary(dotps);
    end
end

% 添加整体标题
sgtitle('Approximation for $u_i^{\top}u_j (i \neq j)$ with Gaussian Distribution', 'Interpreter', 'Latex', 'FontSize', 18);

% Create a separate axis for the legend
legend_handle = axes('Position', [0.15, 0.05, 0.8, 0.05], 'Visible', 'off');
legend(legend_handle, [h1, h2, h3], {'Empirical Frequency for $u_i^\top u_j$',...
                                    '$\mathcal{N}(0, \sigma_{\alpha}^2)$ with $\sigma^2_{\alpha} = r/n^2$',...
                                    '$\mathcal{N}(0, \sigma_{\alpha}^2)$ with $\sigma^2_{\alpha} = \frac{r(n-r)}{n(n-1)(n+2)}$'},...
                    'Interpreter', 'latex', 'Orientation', 'horizontal', 'FontSize', 15, 'Location', 'south');


function summary(array)
    mean_value = mean(array);
    std_value = std(array);
    min_value = min(array);
    max_value = max(array);
    median_value = median(array);

    fprintf('均值: %.3e\n', mean_value);
    fprintf('标准差: %.3e\n', std_value);
    fprintf('最小值: %.3e\n', min_value);
    fprintf('最大值: %.3e\n', max_value);
    fprintf('中位数: %.3e\n', median_value);
    fprintf('\n');
end
