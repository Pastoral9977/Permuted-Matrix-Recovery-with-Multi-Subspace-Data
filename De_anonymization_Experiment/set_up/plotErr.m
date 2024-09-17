color1 = '#0F7002';  % ��ɫ
color2 = '#513ED8';  % ����ɫ
color3 = '#ca6924';  % �س�ɫ
color4 = '#be002f';  % ���ɫ
color5 = '#0000AE';  % ����ɫ
color6 = '#800000';  % ��ɫ
color7 = '#FF5733';  % �Ⱥ�ɫ
color8 = '#33FF57';  % ǳ��ɫ
color9 = '#5733FF';  % ����ɫ
color10 = '#FFD700';  % ��ɫ
color11 = '#00FFFF';  % ��ɫ
color12 = '#FF00FF';  % ���ɫ
color13 = '#800080';  % ��ɫ
color14 = '#808080';  % ��ɫ
color15 = '#FFA500';  % ��ɫ
color16 = '#8B4513';  % ��ɫ
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
