color1 = '#0F7002';  % 绿色
color2 = '#513ED8';  % 深紫色
color3 = '#ca6924';  % 棕橙色
color4 = '#be002f';  % 深红色
color5 = '#0000AE';  % 深蓝色
color6 = '#800000';  % 栗色
color7 = '#FF5733';  % 橙红色
color8 = '#33FF57';  % 浅绿色
color9 = '#5733FF';  % 紫蓝色
color10 = '#FFD700';  % 金色
color11 = '#00FFFF';  % 青色
color12 = '#FF00FF';  % 洋红色
color13 = '#800080';  % 紫色
color14 = '#808080';  % 灰色
color15 = '#FFA500';  % 橙色
color16 = '#8B4513';  % 棕色
close all
sl = length(shuffled_ratio_list);
s_id = 1: sl;
fg = figure;
errorbar( mean_MUL_1, std_MUL_1,'Color',color10,'LineWidth',1.3);hold on;
errorbar( mean_MUL_2, std_MUL_2,'Color',color2,'LineWidth',1.3);hold on;
errorbar( mean_MUL_3, std_MUL_3,'Color',color12,'LineWidth',1.3);hold on;
errorbar( mean_LSR, std_LSR,'Color',color4,'LineWidth',1.3, 'LineStyle', '--','Marker', 'h', 'MarkerSize', 15);hold on;
errorbar( mean_rpca, std_rpca,'Color',color1,'LineWidth',1.3, 'LineStyle', ':');hold on;
errorbar( mean_rkpca, std_rkpca,'Color',color7,'LineWidth',1.3, 'LineStyle', ':');hold on;
errorbar( mean_ssc, std_ssc,'Color',color6,'LineWidth',1.3, 'LineStyle', ':');hold on;
errorbar(mean_tilde, std_tilde,'k-.', 'LineWidth',1.5);hold on;

xlim([1,10]);xticklabels({ '', '0.2', '', '0.4', '','0.6', '', '0.8', '', '1.0'});
lgd = legend('$\textbf{PMSDR}_{L=1}$', '$\textbf{PMSDR}_{L=2}$', '$\textbf{PMSDR}_{L=3}$', 'UPCA', 'RPCA', 'RKPCA', 'SSC', '$\tilde{X}$', ...
    'Interpreter', 'latex', 'Location', 'NorthEastOutside');
lgd.ItemTokenSize = [28, 18];
lgd.NumColumns = 1;
legend('boxoff')
set(gca,'Ydir','Normal', 'FontSize', 18, 'Fontname','times new Roman');
xlabel('shuffled\_ratio','FontSize', 18);

ylabel('$\textbf{\textit{RE}}$','Interpreter', 'Latex')
title('Medical')

title('Educational')
ylim([0.022, 0.2])

fg.Units = 'centimeters';
fg.Position = [10 5 18 11];

ax = gca;
ax.Units = 'centimeters';
ax.Position = [2.5 1.8 9.5 8.2]; 
