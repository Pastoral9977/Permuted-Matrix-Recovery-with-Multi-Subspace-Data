for j = [12, 23, 46, 54]
    mkdir('C:\Users\Pastoral\Desktop\texfile\Gathering\images\fig_face\random', [num2str(sp),'_',num2str(j)])
    path = ['C:\Users\Pastoral\Desktop\texfile\Gathering\images\fig_face\random\',num2str(sp),'_',num2str(j)];
    f_id = 0;
for i = 1:num_groups

img_title = 'raw';
a = image_face(X_raw(:, sum(nn(1:i))+j), f_id, hh, ww, img_title);
eval(['cd ' path]);
str = [num2str(i),'a.png'];
saveas(gcf,str);
cd 'C:\Users\Pastoral\Desktop\Pursuit\Face Experiment'
f_id = f_id + 1;

img_title = 'outlier';
b = image_face(X_tilde(:, sum(nn(1:i))+j), f_id, hh, ww, img_title);
eval(['cd ' path]);
str = [num2str(i),'b.png'];
saveas(gcf,str);
cd 'C:\Users\Pastoral\Desktop\Pursuit\Face Experiment'
f_id = f_id + 1;

img_title = 'Recovered';
c = image_face(X_solved(:,sum(nn(1:i))+j), f_id, hh, ww, img_title);
eval(['cd ' path]);
str = [num2str(i),'c.png'];
saveas(gcf,str);
cd 'C:\Users\Pastoral\Desktop\Pursuit\Face Experiment'
f_id = f_id + 1;

end
end