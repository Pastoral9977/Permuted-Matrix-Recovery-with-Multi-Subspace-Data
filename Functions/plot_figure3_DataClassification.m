
f = figure;

subplot(1,2,1)
imagesc(outlier_ratio_set, shuffled_ratio_set, abs(MISSRATE_mean), [0,0.01]);
c = flipud(gray);colormap(c);
colorbar;
set(gca,'Ydir','Normal', 'FontSize', 12, 'Fontname','times new Roman');
yticks(outlier_ratio_set);
xticks(shuffled_ratio_set);
yticklabels(outlier_ratio_set);
xticklabels(shuffled_ratio_set);
ylabel('outlier-ratio');
xlabel('shuffled-ratio');
title('rank = 4 && noise-free ~ MISSRATE')


subplot(1,2,2)
imagesc(outlier_ratio_set, shuffled_ratio_set, abs(MISSRATE_1_mean), [0,0.3]);
c = flipud(gray);colormap(c);
colorbar;
set(gca,'Ydir','Normal', 'FontSize', 12, 'Fontname','times new Roman');
yticks(outlier_ratio_set);
xticks(shuffled_ratio_set);
yticklabels(outlier_ratio_set);
xticklabels(shuffled_ratio_set);
ylabel('outlier-ratio');
xlabel('shuffled-ratio');
title('rank = 4 && noisy ~ MISSRATE-1')

