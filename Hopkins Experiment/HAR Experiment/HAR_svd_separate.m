close all
clear 
addpath(genpath('C:\Users\Pastoral\Desktop\Copy_SSC_ADMM_v1.1'));
addpath(genpath('C:\Users\Pastoral\Desktop\Unlabeled_PCA_Code'));
addpath(genpath('C:\Users\Pastoral\Desktop\Pursuit\Functions'));
load('HARtrain.mat');

w = max(subjecttrain,[],'all');
idv = cell(1,w);
figure(1);
for i = 1:w
    idv{1,i} = ytrain( find(subjecttrain == i) );
    [~,D,~] = svd(Xtrain(:,idv{1,i}));
    d = diag(D);
    plot(d);hold on;
end
xlim([0,15]);
title('per person');
hold off

h = max(ytrain,[],'all');
act = cell(1,h);
figure(2);
for i = 1:h
    act{1,i} = subjecttrain( find(ytrain == i) );
    [~,D,~] = svd(Xtrain(:,act{1,i}));
    d = diag(D);
    plot(d);hold on;
end
xlim([0,20]);
title('per action');
hold off
