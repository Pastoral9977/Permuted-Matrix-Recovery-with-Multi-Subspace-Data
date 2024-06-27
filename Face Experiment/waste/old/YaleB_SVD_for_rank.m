clear;
load('YaleB_cell.mat');
h = 192; w = 168; 
hh = 48; ww = 42;
X = DSP(YaleB_cell{1,5},h,w,hh,ww);
[~,D,~] = svd(X);
t = 0;r = 0;d = diag(D);
while t/sum(d)<0.75
    r = r+1;
    t = t+d(r);
end
plot(d);
image_face(X(:, 23), 10, hh, ww, 'inlier');




