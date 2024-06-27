clear;
close all;

n_values = [100, 500, 1000]; 
m = 25;
ratio_values = [0.1, 0.3, 0.5, 0.7, 0.9]; 

num_rows = length(n_values);
num_cols = length(ratio_values);

% 设置 figure 的大小
figure('Position', [80, 60, 1100, 700]); 

% Margin settings
left_margin = 0.05;
right_margin = 0.05;
top_margin = 0.12; % Increased top margin to make space for the title
bottom_margin = 0.14; % Increased bottom margin to make space for the legend
h_spacing = 0.05;
v_spacing = 0.08;

% Calculate subplot width and height
subplot_width = (1 - left_margin - right_margin - (num_cols - 1) * h_spacing) / num_cols;
subplot_height = (1 - top_margin - bottom_margin - (num_rows - 1) * v_spacing) / num_rows;

for k = 1:length(n_values)
    n = n_values(k);
    step = n/(m*2);
    r_values = step:step:n/2; 

    variance_front = zeros(length(r_values), length(ratio_values));
    variance_back = zeros(length(r_values), length(ratio_values));

    for i = 1:length(r_values)
        for j = 1:length(ratio_values)
            r = fix(r_values(i));
            ratio = ratio_values(j);
            
            n1 = round(n * (1 - ratio));
            n2 = n - n1;
            
            data = randn(n, r);
            [U, ~] = qr(data, 0);
            
            y = U * randn(r, 1);
            y_tilde = y;
            last_elements = y(end - n2 + 1:end);
            y_tilde(end - n2 + 1:end) = last_elements(randperm(n2));
            
            beta = U \ y_tilde;
            y_tilde_hat = U * beta;
            
            values = abs(y_tilde - y_tilde_hat);
            
            front_norms = values(1:n1);
            back_norms = values(n1+1:end);
            
            variance_front(i, j) = var(front_norms);
            variance_back(i, j) = var(back_norms);
        end
    end

    for j = 1:length(ratio_values)
        % Calculate position for the subplot
        subplot_left = left_margin + (j - 1) * (subplot_width + h_spacing);
        subplot_bottom = 1 - top_margin - k * (subplot_height + v_spacing) + v_spacing;

        % Create subplot
        subplot_handle = subplot('Position', [subplot_left, subplot_bottom, subplot_width, subplot_height]);
        
        hold on;
        h1 = plot(r_values, variance_front(:, j), 'o-', 'DisplayName', 'Unshuffled part', 'MarkerSize', 5);
        h2 = plot(r_values, variance_back(:, j), 'x-', 'DisplayName', 'Shuffled part');
        if k == length(n_values)
            xlabel('r', 'FontSize', 15);
        end
        if j == 1
            ylabel('Variance', 'FontSize', 12);
        end
        ylim([0, max(max(variance_back(:, j)), max(variance_front(:, j)))*1.2]);
        xlim([0, max(r_values)]);
        title(['$n_2/n$ = ' num2str(ratio_values(j)), ', n = ' num2str(n)], 'Interpreter', 'latex', 'FontSize', 15);
        grid on;
        hold off;
    end
end

sgtitle('$\mathrm{Var}[(\tilde{y}-\hat{\tilde{y}})_j]$ between Shuffled Part and Unshuffled Part', 'FontSize', 18, 'Interpreter', 'latex');

% Create a separate axis for the legend
legend_handle = axes('Position', [0.15, 0.12, 0.8, 0.05], 'Visible', 'off');
legend([h1, h2], {'Unshuffled part ($j \leq n_1$))', 'Shuffled part ($j \geq n_1+1$)'},...
    'Interpreter', 'latex', 'Orientation', 'horizontal', 'FontSize', 16, 'Location', 'south');
