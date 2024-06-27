function [Inliers_id, Outliers_id, Z] = divide(X_tilde,hh,ww)
%Parameters for SSC
affine = false ;alpha = 100000;
X_std = X_tilde./vecnorm(X_tilde,1);
X0_std = X_std; 
[X11_std,X12_std,X21_std,X22_std] = halve(X_std,hh,ww);

Z0 = admmLasso_mat_func(X0_std, affine, alpha); 
Z1 = admmLasso_mat_func(X11_std, affine, alpha); 
Z2 = admmLasso_mat_func(X12_std, affine, alpha); 
Z3 = admmLasso_mat_func(X21_std, affine, alpha); 
Z4 = admmLasso_mat_func(X22_std, affine, alpha); 

value0 = vecnorm(Z0,1);
value1 = vecnorm(Z1,1);
value2 = vecnorm(Z2,1);
value3 = vecnorm(Z3,1);
value4 = vecnorm(Z4,1);

T0 = kmeans(value0',2);
T1 = kmeans(value1',2);
T2 = kmeans(value2',2);
T3 = kmeans(value3',2);
T4 = kmeans(value4',2);

I0 = find(T0 == 1); J0 = find(T0 == 2);
I1 = find(T1 == 1); J1 = find(T1 == 2);
I2 = find(T2 == 1); J2 = find(T2 == 2);
I3 = find(T3 == 1); J3 = find(T3 == 2);
I4 = find(T4 == 1); J4 = find(T4 == 2);

I0_mean = mean(value0(:,I0));J0_mean = mean(value0(:,J0));
I1_mean = mean(value1(:,I1));J1_mean = mean(value1(:,J1));
I2_mean = mean(value2(:,I2));J2_mean = mean(value2(:,J2));
I3_mean = mean(value3(:,I3));J3_mean = mean(value3(:,J3));
I4_mean = mean(value4(:,I4));J4_mean = mean(value4(:,J4));

y0 = max(I0_mean,J0_mean)/min(I0_mean,J0_mean);
y1 = max(I1_mean,J1_mean)/min(I1_mean,J1_mean);
y2 = max(I2_mean,J2_mean)/min(I2_mean,J2_mean);
y3 = max(I3_mean,J3_mean)/min(I3_mean,J3_mean);
y4 = max(I4_mean,J4_mean)/min(I4_mean,J4_mean);

if y0 == max([y0,y1,y2,y3,y4],[],'all')
    I_mean = I0_mean;J_mean = J0_mean;
    I = I0; J = J0; Z = Z0;
elseif y1 == max([y0,y1,y2,y3,y4],[],'all')
    I_mean = I1_mean;J_mean = J1_mean;
    I = I1; J = J1; Z = Z1;
elseif y2 == max([y0,y1,y2,y3,y4],[],'all')
    I_mean = I2_mean;J_mean = J2_mean;
    I = I2; J = J2; Z = Z2;
elseif y3 == max([y0,y1,y2,y3,y4],[],'all')
    I_mean = I3_mean;J_mean = J3_mean;
    I = I3; J = J3; Z = Z3;
elseif y4 == max([y0,y1,y2,y3,y4],[],'all')
    I_mean = I4_mean;J_mean = J4_mean;
    I = I4; J = J4; Z = Z4;
end

if (I_mean < J_mean)
    Inliers_id = J;Outliers_id = I;
else
    Inliers_id = I;Outliers_id = J;
end
end
    








