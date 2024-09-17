close all
clear
clc

M_values = [100, 500, 1000];
M2_ratios = [0.1, 0.3, 0.5, 0.7, 0.9];
step_ratio = 25;

sigma_1 = @(M, M2, r) sqrt(2.*r.*(M-r)/(M^2*(M+2)));
sigma_2 = @(M, M2, r) sqrt(r.*(M-r)/(M*(M-1)*(M+2)));
kappa = @(m) norminv(1-1/m, 0, 1);
tau = @(m) 1 / (normpdf(kappa(m), 0, 1)*m);
% sigma_xi = @(M, M2, r) sqrt( 2*M2./r.*sigma_2(M, M2, r).^2.*(r/M+3*sigma_1(M, M2, r)+3*sigma_2(M, M2, r)) );
% sigma_xi  = @(M, M2, r) sqrt( UB(M, M2, r) );
% sigma_eta = @(M, M2, r) sqrt( (2*M-M2)/M * 1/M .* ( (r/M-1).^2 + 2.*r.*(M-r)/(M^2.*(M+2)) + (M2-1).*r.*(M-r)/(M*(M-1)*(M+2)) ) );

sigma_xi = @(M, M2, r) sqrt((1+3*(1+sqrt(2)).*sqrt((M-r)./(r*M2*M))).* r.*(M-r).*M2./(M^2*(M-1)*(M+2)) + (M2/M)^2./M) ;
sigma_eta = @(M, M2, r) sqrt(((2*M-M2)/(M*2))* 2/M .* ( (r/M-1).^2 + 2.*r.*(M-r)/(M^2.*(M+2)) + (M2-1).*r.*(M-r)/(M*(M-1)*(M+2))) );

% Initialize figure
figure;

% Outer loop over each M value
for j = 1:length(M_values)
    M = M_values(j);
    step = M/step_ratio;
    r = step:step:M;
    
    % Inner loop over each M2 value
    for i = 1:length(M2_ratios)
        M2 = M2_ratios(i) * M;
        
        % Calculate sigma_xi and sigma_eta
        s_xi = sigma_xi(M, M2, r);
        s_eta = sigma_eta(M, M2, r);
        
        % Calculate phi and res
        phi = (s_eta * kappa(M2) - s_xi * kappa(M - M2)) ./ sqrt(s_xi.^2 * tau(M - M2)^2 + s_eta.^2 * tau(M2)^2);
        res = 2*normcdf(phi)-1;
        
        % Plot res vs r in the corresponding subplot
        subplot(3, 1, j);
        plot(r, res, 'DisplayName', ['M_2 = ', num2str(M2)], 'linewidth', 2);
        hold on;
    end
    
    % Final touches for res plot
    subplot(3, 1, j);
    title(['$', '\Pr\left(j \geq M_1+1\right)$ vs $r$ for $M = ', num2str(M), '$'], 'Interpreter', 'latex');
    xlabel('r');
    ylabel('Probability');
    legend('show', 'Location', 'best');
    hold off;
end