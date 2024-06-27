
% 考察列正交阵Q的行模长的平方：Q*Q'的对角元

clear
close all

n_values = [100, 500, 1000]; 
rn_ratios = [1/0.9, 2, 10, 25];
shuffled_ratio = 0.4;
num_samples = 500;  

num_rows = length(n_values);
num_cols = length(rn_ratios);

figure;

% Margin settings
left_margin = 0.08;
right_margin = 0.08;
top_margin = 0.12; % Increased top margin to make space for the title
bottom_margin = 0.15; % Increased bottom margin to make space for the legend
h_spacing = 0.05;
v_spacing = 0.05;

% Calculate subplot width and height
subplot_width = (1 - left_margin - right_margin - (num_cols - 1) * h_spacing) / num_cols;
subplot_height = (1 - top_margin - bottom_margin - (num_rows - 1) * v_spacing) / num_rows;

for i = 1:length(n_values)
    for j = 1:length(rn_ratios)
        n = n_values(i);
        r = round(n / rn_ratios(end - j + 1));
        n2 = round(n * shuffled_ratio); 

        % 存储所有采样得到的 f 值
        all_f = zeros(num_samples * n2, 1);

        for k = 1:num_samples
            % 生成随机矩阵并进行QR分解
            A = randn(n, r) / sqrt(n);  % 确保 A 的大小是 n x r
            [Q, ~] = qr(A, 0);
            Q = Q(end-n2+1:end, :);

            % 计算 F 矩阵和相关量
            F = Q * Q';
            f = diag(F);

            % 将当前采样得到的 f 值存储起来
            start_idx = (k-1) * n2 + 1;
            end_idx = k * n2;
            all_f(start_idx:end_idx) = f;
        end
        
        % Calculate position for the subplot
        subplot_left = left_margin + (j - 1) * (subplot_width + h_spacing);
        subplot_bottom = 1 - top_margin - i * (subplot_height + v_spacing) + v_spacing;

        % Create subplot
        subplot_handle = subplot('Position', [subplot_left, subplot_bottom, subplot_width, subplot_height]);
        hold on

        % 设置直方图的边界和显示
        edges = linspace(0, 1, 150 + 1);
        h1 = histogram(all_f, 'BinEdges', edges, 'Normalization', 'probability');
        
        % 设置 Beta 分布的参数
        alpha = r/2;  % 形状参数 alpha
        beta = (n - r)/2;   % 形状参数 beta

        % 生成 x 轴的值（从 0 到 1）
        x = linspace(0, 1, 1000);

        % 计算 Beta 分布的概率密度函数值
        pdf_values = betapdf(x, alpha, beta);

        % 绘制 Beta 分布的概率密度函数图像
        h2 = plot(x, pdf_values * (edges(2) - edges(1)), 'LineWidth', 2);
        xlim([0, 1]);
        ylim([0, max(pdf_values) * (edges(2) - edges(1)) * 1.1]);

        % 设置标题和图例
        title(['$n = ' num2str(n) ', \ r = ' num2str(r) '$'], 'Interpreter', 'latex');
        if i == length(n_values)
            xlabel('x');
        end
        if j == 1
            ylabel('Probability');
        end
%         legend({'Empirical Histogram', ...
%                 ['$\alpha = ' num2str(alpha) ', \ \beta = ' num2str(beta) '$']}, ...
%                 'Interpreter', 'latex', 'Location', 'best');
        grid on;
        hold off
    end
end

sgtitle("Empirical Frequency and Distribution Estimation for $||Q(:,i)||_2$", "Interpreter", "latex", "FontSize", 16);
% Create a separate axis for the legend
legend_handle = axes('Position', [0.1, 0.05, 0.8, 0.05], 'Visible', 'off');
legend(legend_handle, [h1, h2], {'Empirical Frequency', 'Fitted Beta PDF $\beta(\frac{r}{2}, \frac{n-r}{2})$'},...
                    'Interpreter', 'latex', 'Orientation', 'horizontal', 'FontSize', 12, 'Location', 'south');
