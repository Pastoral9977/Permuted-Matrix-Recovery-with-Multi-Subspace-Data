f = figure;

subplot(1,2,1)
imagesc(outlier_ratio_set, shuffled_ratio_set, abs(err_ratio_mean));
c = flipud(gray);
colormap(c);
colorbar;
set(gca,'Ydir','Normal', 'FontSize', 12, 'Fontname','times new Roman');
xticks(outlier_ratio_set);
yticks(shuffled_ratio_set);
xticklabels(outlier_ratio_set);
yticklabels(shuffled_ratio_set);
xlabel('outlier-ratio');
ylabel('shuffled-ratio');
title('rank = 4 && noise-free ~ fro-err-ratio')



subplot(1,2,2)
imagesc(outlier_ratio_set, shuffled_ratio_set, abs(err_ratio_1_mean));
c = flipud(gray);colormap(c);
colorbar;
set(gca,'Ydir','Normal', 'FontSize', 12, 'Fontname','times new Roman');
xticks(outlier_ratio_set);
yticks(shuffled_ratio_set);
xticklabels(outlier_ratio_set);
yticklabels(shuffled_ratio_set);
xlabel('outlier-ratio');
ylabel('shuffled-ratio');
title('rank = 4 && noisy ~ fro-err-ratio-1')


