clear
close all

M_values = [100, 500, 1000]; 
rM_ratios = [1/0.9, 2, 10, 25];
shuffled_ratio = 0.4;
num_samples = 8000;  

num_rows = length(rM_ratios);
num_cols = length(M_values);

figure('Position', [80, 60, 900, 600]); 

% Margin settings
left_margin = 0.05;
right_margin = 0.05;
top_margin = 0.12; % Increased top margin to make space for the title
bottom_margin = 0.20; % Increased bottom margin to make space for the legend
h_spacing = 0.06;
v_spacing = 0.06;

% Calculate subplot width and height
subplot_width = (1 - left_margin - right_margin - (num_cols - 1) * h_spacing) / num_cols;
subplot_height = (1 - top_margin - bottom_margin - (num_rows - 1) * v_spacing) / num_rows;

for i = 1:length(rM_ratios)
    for j = 1:length(M_values)
        M = M_values(j);
        r = round(M / rM_ratios(end - i + 1));
        M2 = round(M * shuffled_ratio); 

        % Store all sampled f values
        all_f = zeros(num_samples, 1);
        
        mu_sum = M2 * r / M;
        var_sum = 2 * M2 * (M - M2) * r * (M - r) / (M^2 * (M + 2) * (M - 1));
        std_sum = sqrt(var_sum);
        
        for k = 1:num_samples
            % Generate a random matrix and perform QR decomposition
            A = randn(M, r) / sqrt(M);  % Ensure A is M x r
            [Q, ~] = qr(A, 0);
            Q = Q(end-M2+1:end, :);

            % Calculate F matrix and related quantities
            F = Q * Q';
            f = trace(F);

            all_f(k) = f;
        end

        % Calculate position for the subplot
        subplot_left = left_margin + (j - 1) * (subplot_width + h_spacing);
        subplot_bottom = 1 - top_margin - i * (subplot_height + v_spacing) + v_spacing;

        % Create subplot
        subplot_handle = subplot('Position', [subplot_left, subplot_bottom, subplot_width, subplot_height]);
        
        hold on

        % Set histogram boundaries and display
        edges = linspace(mu_sum - 5.5 * std_sum, mu_sum + 4.5 * std_sum, 120 + 1);
        h1 = histogram(all_f, 'BinEdges', edges, 'Normalization', 'probability');

        % Add vertical lines for mean and ¡À3¦Ò positions
        xline(mu_sum, '--r', 'LineWidth', 1.5, 'Label', '\mu', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top', 'FontSize', 12);
        xline(mu_sum - 3 * std_sum, '--b', 'LineWidth', 1.5, 'Label', '\mu-3\sigma', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top', 'FontSize', 12);
        xline(mu_sum + 3 * std_sum, '--b', 'LineWidth', 1.5, 'Label', '\mu+3\sigma', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top', 'FontSize', 12);
        
        % Estimate new Beta distribution parameters
        [alpha_prime, beta_prime] = beta_parameters(mu_sum/M2, var_sum/M2^2);

        if alpha_prime > 0 && beta_prime > 0  % Ensure valid parameters
            x = linspace(0, 1, 100*M2);
            beta_pdf = betapdf(x, alpha_prime, beta_prime);
            beta_proba = beta_pdf * (edges(2) - edges(1))/M2;
            h2 = plot(x*M2, beta_proba, 'r-', 'LineWidth', 1.8);
            
            % Add normal distribution curve
            normal_x = linspace(min(edges), max(edges), 1000);
            normal_pdf = normpdf(normal_x, mu_sum, std_sum);
            normal_proba = normal_pdf*(edges(2) - edges(1));
            h3 = plot(normal_x, normal_proba, 'y--', 'LineWidth', 1.8);
            
        else
            warning('Invalid alpha_prime or beta_prime for n = %d, r = %d', M, r);
        end
        
        % Set title
        title(['$M = ' num2str(M) ', \ r = ' num2str(r) '$'], 'Interpreter', 'latex');
        xlim([min(edges), max(edges)]);
        ylim([0, max(h1.Values) * 1.2]);
        if i == length(rM_ratios)
            xlabel('$\mathrm{tr}(\mathbf{U}^{(2)}{\mathbf{U}^{(2)}}^\top)$','Interpreter', 'latex', 'FontSize', 15);
        end
        if j == 1
            ylabel('Probability','Interpreter', 'latex', 'FontSize', 15);
        end
        grid on;
        hold off
    end
end

% Create a separate axis for the legend
legend_handle = axes('Position', [0.15, 0.05, 0.8, 0.05], 'Visible', 'off');
legend(legend_handle, [h1, h2, h3],...
    {'Empirical Frequency','Fitted Beta PDF $\beta(\frac{(1-\mu)*\mu^2}{\sigma^2}, \frac{(1-\mu)*\mu^2*(1-\mu)}{\sigma^2*\mu})$','Fitted Gaussian PDF $\mathcal{N}(\mu, \sigma^2)$'},...
    'Interpreter', 'latex', 'Orientation', 'horizontal', 'FontSize', 15, 'Location', 'south');

sgtitle("Approximation for $\mathrm{tr}(\mathbf{U}^{(2)}{\mathbf{U}^{(2)}}^\top)$ with Beta and Gaussian Distribution Given $\mu$ and $\sigma^2$", 'Interpreter', 'latex', 'FontSize', 18);

function [alpha, beta] = beta_parameters(mu, variance)
    % Check for valid mean and variance values
    if mu <= 0 || mu >= 1
        error('Mean must be between 0 and 1 (exclusive).');
    end
    
    if variance <= 0
        error('Variance must be positive.');
    end

    % Compute alpha and beta
    alpha = ((1 - mu) * mu^2 / variance) - mu;
    beta = alpha * (1 - mu) / mu;
end
