
figure;

subplot(1,3,1)
imagesc(rank_set, shuffled_ratio_set, abs(err_IN_mean),[0,0.2]);
c = flipud(gray);colormap(c);
colorbar;
set(gca,'Ydir','Normal', 'FontSize', 25, 'Fontname','times new Roman');
yticks(0.15:0.09:5);
xticks(2:6:50);
yticklabels([5,10,15,20,25,30]);
xticklabels(shuffled_ratio_set);
ylabel('rank');
xlabel('shuffled-ratio');
title('Err-in')


subplot(1,3,2)
imagesc(rank_set, shuffled_ratio_set, abs(missrate_OUT_mean),[0,0.2]);
c = flipud(gray);colormap(c);
colorbar;
set(gca,'Ydir','Normal', 'FontSize', 25, 'Fontname','times new Roman');
yticks(0.15:0.09:5);
xticks(2:6:50);
yticklabels([5,10,15,20,25,30]);
xticklabels(shuffled_ratio_set);
ylabel('rank');
xlabel('shuffled-ratio');
title('missrate-out')


subplot(1,3,3)
imagesc(rank_set, shuffled_ratio_set, abs(err_ratio_dir_mean),[0,0.2]);
c = flipud(gray);colormap(c);
colorbar;
set(gca,'Ydir','Normal', 'FontSize', 25, 'Fontname','times new Roman');
yticks(0.15:0.09:5);
xticks(2:6:50);
yticklabels([5,10,15,20,25,30]);
xticklabels(shuffled_ratio_set);
ylabel('rank');
xlabel('shuffled-ratio');
title('outliers recovery')


