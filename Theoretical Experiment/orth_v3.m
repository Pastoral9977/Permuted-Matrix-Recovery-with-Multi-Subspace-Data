close all
clear

n = 50;
shuffed_ratio = 0.5;
r = 25;
n1 = fix(n*(1-shuffed_ratio));
n2 = n-n1;

data = randn(n, r);
[U, ~] = qr(data,0);

y = U*randn(r,1);
y_tilde = y;
last_elements = y(end-n2+1:end);
y_tilde(end-n2+1:end) = last_elements(randperm(n2));

beta = U\y_tilde;
y_tilde_hat = U*beta;

figure;
values = abs(y_tilde-y_tilde_hat);
plot(values, 'Marker','.', 'MarkerSize', 10)
x_line = [n1, n1];  
y_line = [min(values), max(values)*1.3]; 
line(x_line, y_line, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1.5); 
grid on

% 添加title、xlabel、ylabel和legend
title('Residuals between $\tilde{y}$ and $\hat{\tilde{y}}$', 'Interpreter', 'latex')
xlabel('Index', 'Interpreter', 'latex')
ylabel('Absolute Residual', 'Interpreter', 'latex')
legend({'$|\tilde{y} - \hat{\tilde{y}}|$', 'Shuffled Boundary'}, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize',  15)

% 在图中显示 r, n2 和 n
dim = [.15 .6 .3 .3];
sigma1 = r/n*sqrt(n2/n);
sigma2 = r/n;
str = sprintf('$r = %d$, $n_2 = %d$, $n = %d$\n$\\sigma_1 = \\frac{r}{n}\\sqrt{\\frac{n_2}{n}} = %.4f$\n$\\sigma_2 = \\frac{r}{n} = %.4f$', r, n2, n, sigma1, sigma2);
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'Interpreter', 'latex', 'BackgroundColor', 'white');

% 计算没有乱序部分和乱序部分的标准差
std_no_shuffle = std(values(1:n1));
std_shuffle = std(values(n1+1:end));

% 比较标准差比值和 sqrt(n2/n)
std_ratio = std_shuffle / std_no_shuffle;
diff_ratio = std_ratio / sqrt(n2/n);

% % 显示标准差比较结果
% dim2 = [.15 .4 .3 .2];
% str2 = sprintf('No shuffle std = %.4f\nShuffle std = %.4f\nstd ratio = %.4f\n$\\sqrt{\\frac{n_2}{n}}$ = %.4f\ndiff ratio = %.4f', ...
%                std_no_shuffle, std_shuffle, std_ratio, sqrt(n2/n), diff_ratio);
% annotation('textbox', dim2, 'String', str2, 'FitBoxToText', 'on', 'Interpreter', 'latex', 'BackgroundColor', 'white');

