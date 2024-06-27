close all
clear 
addpath(genpath('C:\Users\Pastoral\Desktop\Copy_SSC_ADMM_v1.1'));
addpath(genpath('C:\Users\Pastoral\Desktop\Unlabeled_PCA_Code'));
addpath(genpath('C:\Users\Pastoral\Desktop\Pursuit\Functions'));
addpath(genpath('F:\Github download\LibADMM-master'));
addpath(genpath('F:\Matlab toolbox\lrr'));

load('HARtrain.mat');

[n,N] = size(Xtrain);
rank = 12;  num_groups = max(ytrain);
X_tilde = Xtrain;

Y = X_tilde;
t = ytrain';

% LRR 
lambda = 0.15;
[Z,E] = lrr(Y,lambda);
[idx] = clu_ncut(Z,num_groups);
[acc] = compacc(idx,t);




