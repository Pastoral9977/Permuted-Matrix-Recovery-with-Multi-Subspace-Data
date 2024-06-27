color1 = '#0F7002';
color2 = '#513ED8';
color5 = '#0000AE';
color4 = '#800000';
color6 = '#be002f';
color3 = '#ca6924';

sl = length(shuffled_ratio_list);
s_id = 1: sl;
fg = figure;
errorbar( mean_MUL_1, std_MUL_1,'Color',color1,'LineWidth',1.3);hold on;
errorbar( mean_MUL_2, std_MUL_2,'Color',color2,'LineWidth',1.3);hold on;
errorbar( mean_MUL_3, std_MUL_3,'Color',color3,'LineWidth',1.3);hold on;
errorbar( mean_LSR, std_LSR,'Color',color6,'LineWidth',1.3);hold on;
errorbar(mean_tilde, std_tilde,'k--', 'LineWidth',1.5);hold on;
xlim([1,10]);xticklabels({ '', '0.2', '', '0.4', '','0.6', '', '0.8', '', '1.0'});
% lgd = legend('$MUL_{L=1}$','$MUL_{L=2}$','$MUL_{L=3}$', 'upca', '$\tilde{X}$',  ...
%     'Interpreter', 'Latex', 'Location', 'NorthEastOutside');
% lgd.NumColumns = 1;
% legend('boxoff')
set(gca,'Ydir','Normal', 'FontSize', 18, 'Fontname','times new Roman');
xlabel('\alpha','FontSize', 18);
h = ylabel('$\textbf{\textit{RE\quad\quad}}$','Interpreter', 'Latex');
set(h, 'Rotation', pi/2, 'Position', get(h, 'Position')-[0.8,0,0] );
% title('Educational')
fg.Units = 'centimeters';
fg.Position = [5 5 12 9];

ax = gca;
ax.Units = 'centimeters';
ax.Position = [3 1.8 7. 6];  % 根据需要调整这个值
