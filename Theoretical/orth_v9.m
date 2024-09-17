clear;
clc;
close all;

% M_values = [20, 50, 200];
M_values = [100, 500, 1000];
m = 25;
ratio_values = [0.1, 0.3, 0.5, 0.7, 0.9]; 
T = 10;  % Number of trials

num_rows = length(M_values);
num_cols = length(ratio_values);

figure('Position', [80, 60, 1400, 700]); 

% Margin settings
left_margin = 0.05;
right_margin = 0.05;
top_margin = 0.12; % Increased top margin to make space for the title
bottom_margin = 0.14; % Increased bottom margin to make space for the legend
h_spacing = 0.05;
v_spacing = 0.09;

% Calculate subplot width and height
subplot_width = (1 - left_margin - right_margin - (num_cols - 1) * h_spacing) / num_cols;
subplot_height = (1 - top_margin - bottom_margin - (num_rows - 1) * v_spacing) / num_rows;

sigma_xi_squared  = @(M, M2, r) (1+3*(1+sqrt(2)).*sqrt((M-r)./(r*M2*M))).* r.*(M-r).*M2./(M^2*(M-1)*(M+2)) + (M2/M)^2./M ;
sigma_eta_squared = @(M, M2, r) ((2*M-M2)/(M*2))* 2/M .* ( (r/M-1).^2 + 2.*r.*(M-r)/(M^2.*(M+2)) + (M2-1).*r.*(M-r)/(M*(M-1)*(M+2)) );

for k = 1:length(M_values)
    M = M_values(k);
    step = M/m;
    r_values = [step:step:M];
%     r_values = 1:7:7*step;

    variance_front = zeros(length(r_values), length(ratio_values), T);
    variance_back = zeros(length(r_values), length(ratio_values), T);
    variance_front_estimated = zeros(length(r_values), length(ratio_values));
    variance_back_estimated = zeros(length(r_values), length(ratio_values));
    
    for t = 1:T
        for i = 1:length(r_values)
            for j = 1:length(ratio_values)
                r = fix(r_values(i));
                ratio = ratio_values(j);

                M2 = round(M * ratio);
                M1 = M - M2;
                
                fprintf("M1=%d, M2=%d, r=%d\n\t",M1, M2, r)
                data = randn(M, r);
                [U, ~] = qr(data, 0);
                
                beta_init = randn(r, 1)/ sqrt(r);
                y = U * beta_init;
                fprintf("norm(y)=%.4f\n\t",norm(y))
                y_tilde = y;
                last_elements = y(M1+1:end);
                deranged_indices = derangement(M2);
                y_tilde(M1 + 1:end) = last_elements(deranged_indices);
                fprintf("norm(y_tilde)=%.4f\n\t",norm(y_tilde))
                
                beta = U \ y_tilde;
                y_tilde_hat = U * beta;
                fprintf("norm(y_tilde_hat)=%.4f\n\n",norm(y_tilde_hat))
                
                values = y_tilde - y_tilde_hat;

                front_norms = values(1:M1);
                back_norms = values(M1+1:end);
                variance_front(i, j, t) = var(front_norms);
                variance_back(i, j, t) = var(back_norms);
                
                variance_front_estimated(i, j) = sigma_xi_squared(M, M2, r);
                variance_back_estimated(i, j) = sigma_eta_squared(M, M2, r);
            end
        end
    end

    for j = 1:length(ratio_values)
        % Calculate position for the subplot
        subplot_left = left_margin + (j - 1) * (subplot_width + h_spacing);
        subplot_bottom = 1 - top_margin - k * (subplot_height + v_spacing) + v_spacing;

        % Create subplot
        subplot_handle = subplot('Position', [subplot_left, subplot_bottom, subplot_width, subplot_height]);
        
        hold on;
        % Calculate mean and std over T trials
        front_mean = mean(variance_front(:, j, :), 3);
        front_std = std(variance_front(:, j, :), 0, 3);
        back_mean = mean(variance_back(:, j, :), 3);
        back_std = std(variance_back(:, j, :), 0, 3);

        % Plot mean with shaded std for front
        fill([r_values, fliplr(r_values)], ...
            [front_mean - 2 * front_std; flipud(front_mean + 2 * front_std)], ...
            'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Unshuffled part ¡À2*std');
        h1 = plot(r_values, front_mean, 'o-', 'DisplayName', 'Unshuffled part', 'MarkerSize', 5);

        % Plot mean with shaded std for back
        fill([r_values, fliplr(r_values)], ...
            [back_mean - 2 * back_std; flipud(back_mean + 2 * back_std)], ...
            'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Shuffled part ¡À2*std');
        h2 = plot(r_values, back_mean, 'x-', 'DisplayName', 'Shuffled part');
        
        h3 = plot(r_values, variance_front_estimated(:, j), 'linewidth',2);
        h4 = plot(r_values, variance_back_estimated(:, j), 'linewidth',2);

        if k == length(M_values)
            xlabel('r', 'FontSize', 15);
        end
        if j == 1
            ylabel('Variance', 'FontSize', 12);
        end
        ylim([0, max(max(back_mean + 3 * back_std), max(front_mean + 3 * front_std)) * 1.4]);
        ax = gca;      
        ax.YAxis.Exponent = 0; 
        ytickformat('%.3f')
        xlim([0, max(r_values)]);
        title(['$M_2/M$ = ' num2str(ratio_values(j)), ', $M$ = ' num2str(M)], 'Interpreter', 'latex', 'FontSize', 15);
        grid on;
        hold off;
    end
end

sgtitle('$\mathrm{Var}[(\mathbf{\tilde{y}}-\mathbf{\hat{\tilde{y}}})_j]$ between Shuffled Part and Unshuffled Part', 'FontSize', 18, 'Interpreter', 'latex');

% Create a separate axis for the legend
legend([h1, h2, h3, h4], ...
    {'Unshuffled part ($j \leq M_1$))', 'Shuffled part ($j \geq M_1+1$)', 'Modeled unshuffled part ($j \leq M_1$))', 'Modeled shuffled part ($j \geq M_1+1$)'},...
    'Interpreter', 'latex', 'Orientation', 'horizontal', 'FontSize', 16, 'Location', 'south', 'Position', [0.28, 0.015, 0.4, 0.05]);

function d = derangement(n)
    % Generates a derangement of 1:n
    while true
        p = randperm(n);
        if all(p ~= 1:n)
            d = p;
            return;
        end
    end
end
