clear;
load('HARtrain.mat');
X = Xtrain(:,ytrain == 1);
[U,D,~] = svd(X);
t = 0;r = 0;d = diag(D);
Y = U(:,1:10)*U(:,1:10)'*X;
[~,D1,~] = svd(Y);d1 = diag(D1);
while t/sum(d)<0.75
    r = r+1;
    t = t+d(r);
end
plot(d);
hold on
plot(d1);
