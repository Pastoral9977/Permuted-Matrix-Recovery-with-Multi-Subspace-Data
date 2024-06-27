
% ������������Q����ģ����ƽ����Q*Q'�ĶԽ�Ԫ

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

        % �洢���в����õ��� f ֵ
        all_f = zeros(num_samples * n2, 1);

        for k = 1:num_samples
            % ����������󲢽���QR�ֽ�
            A = randn(n, r) / sqrt(n);  % ȷ�� A �Ĵ�С�� n x r
            [Q, ~] = qr(A, 0);
            Q = Q(end-n2+1:end, :);

            % ���� F ����������
            F = Q * Q';
            f = diag(F);

            % ����ǰ�����õ��� f ֵ�洢����
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

        % ����ֱ��ͼ�ı߽����ʾ
        edges = linspace(0, 1, 150 + 1);
        h1 = histogram(all_f, 'BinEdges', edges, 'Normalization', 'probability');
        
        % ���� Beta �ֲ��Ĳ���
        alpha = r/2;  % ��״���� alpha
        beta = (n - r)/2;   % ��״���� beta

        % ���� x ���ֵ���� 0 �� 1��
        x = linspace(0, 1, 1000);

        % ���� Beta �ֲ��ĸ����ܶȺ���ֵ
        pdf_values = betapdf(x, alpha, beta);

        % ���� Beta �ֲ��ĸ����ܶȺ���ͼ��
        h2 = plot(x, pdf_values * (edges(2) - edges(1)), 'LineWidth', 2);
        xlim([0, 1]);
        ylim([0, max(pdf_values) * (edges(2) - edges(1)) * 1.1]);

        % ���ñ����ͼ��
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
