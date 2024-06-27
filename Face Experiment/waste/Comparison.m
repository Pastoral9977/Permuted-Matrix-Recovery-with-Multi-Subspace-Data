
%Parameters for SSC
affine = false ;alpha = 100000;
X_std = X_tilde;
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

if (I0_mean < J0_mean)
    Inliers0_id = J0;Outliers0_id = I0;
else
    Inliers0_id = I0;Outliers0_id = J0;
end

if (I1_mean < J1_mean)
    Inliers1_id = J1;Outliers1_id = I1;
else
    Inliers1_id = I1;Outliers1_id = J1;
end

if (I2_mean < J2_mean)
    Inliers2_id = J2;Outliers2_id = I2;
else
    Inliers2_id = I2;Outliers2_id = J2;
end

if (I3_mean < J3_mean)
    Inliers3_id = J3;Outliers3_id = I3;
else
    Inliers3_id = I3;Outliers3_id = J3;
end

if (I4_mean < J4_mean)
    Inliers4_id = J4;Outliers4_id = I4;
else
    Inliers4_id = I4;Outliers4_id = J4;
end

err0_in = length(setdiff(Inliers0_id,inliers_ID))/length(Inliers0_id);
err0_out = length(setdiff(Outliers0_id,outliers_ID))/length(Outliers0_id);

err1_in = length(setdiff(Inliers1_id,inliers_ID))/length(Inliers1_id);
err1_out = length(setdiff(Outliers1_id,outliers_ID))/length(Outliers1_id);

err2_in = length(setdiff(Inliers2_id,inliers_ID))/length(Inliers2_id);
err2_out = length(setdiff(Outliers2_id,outliers_ID))/length(Outliers2_id);

err3_in = length(setdiff(Inliers3_id,inliers_ID))/length(Inliers3_id);
err3_out = length(setdiff(Outliers3_id,outliers_ID))/length(Outliers3_id);

err4_in = length(setdiff(Inliers4_id,inliers_ID))/length(Inliers4_id);
err4_out = length(setdiff(Outliers4_id,outliers_ID))/length(Outliers4_id);



    











