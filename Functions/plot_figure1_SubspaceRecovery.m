close all
f1 = figure;

subplot(3,2,1);
imagesc(shuffled_ratio_set, outlier_ratio_set, abs(LargestDeg_mean_max));
c1 = flipud(gray);colormap(c1);
colorbar;
set(gca,'Ydir','Normal', 'FontSize', 12, 'Fontname','times new Roman');
yticks(outlier_ratio_set);
xticks(shuffled_ratio_set);
yticklabels(outlier_ratio_set);
xticklabels(shuffled_ratio_set);
ylabel('outlier-ratio');
xlabel('shuffled-ratio');
title('rank = 4 && GroupMax && Noise-free ~ LargestDeg')



subplot(3,2,2);
imagesc(shuffled_ratio_set, outlier_ratio_set, abs(LargestDeg1_mean_max));
c2 = flipud(gray);colormap(c2);
colorbar;
set(gca,'Ydir','Normal', 'FontSize', 12, 'Fontname','times new Roman');
yticks(outlier_ratio_set);
xticks(shuffled_ratio_set);
yticklabels(outlier_ratio_set);
xticklabels(shuffled_ratio_set);
ylabel('outlier-ratio');
xlabel('shuffled-ratio');
title('rank = 4 && GroupMax && SNoisy ~ LargestDeg')


subplot(3,2,3);
imagesc(shuffled_ratio_set, outlier_ratio_set, abs(MissRate_mean),[0,0.01]);
c3 = flipud(gray);colormap(c3);
colorbar;
set(gca,'Ydir','Normal', 'FontSize', 12, 'Fontname','times new Roman');
yticks(outlier_ratio_set);
xticks(shuffled_ratio_set);
yticklabels(outlier_ratio_set);
xticklabels(shuffled_ratio_set);
ylabel('outlier-ratio');
xlabel('shuffled-ratio');
title('rank = 4 && Noise-free ~missrate ')


subplot(3,2,4);
imagesc(shuffled_ratio_set, outlier_ratio_set, abs(MissRate1_mean),[0,0.01]);
c4 = flipud(gray);colormap(c4);
colorbar;
set(gca,'Ydir','Normal', 'FontSize', 12, 'Fontname','times new Roman');
yticks(outlier_ratio_set);
xticks(shuffled_ratio_set);
yticklabels(outlier_ratio_set);
xticklabels(shuffled_ratio_set);
ylabel('outlier-ratio');
xlabel('shuffled-ratio');
title('rank = 4 && Noisy ~missrate ')


subplot(3,2,5);
imagesc(shuffled_ratio_set, outlier_ratio_set, abs(Err_sel_mean), [0,1]);
c5 = flipud(gray);colormap(c5);
colorbar;
set(gca,'Ydir','Normal', 'FontSize', 12, 'Fontname','times new Roman');
yticks(outlier_ratio_set);
xticks(shuffled_ratio_set);
yticklabels(outlier_ratio_set);
xticklabels(shuffled_ratio_set);
ylabel('outlier-ratio');
xlabel('shuffled-ratio');
title('rank = 4 && Noisy ~ecc-sel ')

subplot(3,2,6)
imagesc(shuffled_ratio_set, outlier_ratio_set, abs(Err_sel_mean_1),[0,1]);
c6 = flipud(gray);colormap(c6);
colorbar;
set(gca,'Ydir','Normal', 'FontSize', 12, 'Fontname','times new Roman');
yticks(outlier_ratio_set);
xticks(shuffled_ratio_set);
yticklabels(outlier_ratio_set);
xticklabels(shuffled_ratio_set);
ylabel('outlier-ratio');
xlabel('shuffled-ratio');
title('rank = 4 && Noisy ~err-sel ')





