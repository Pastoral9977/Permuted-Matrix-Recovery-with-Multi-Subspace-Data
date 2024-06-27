%--------------------------------------------------------------------------
% This is the main function to run the SSC algorithm for the motion
% segmentation problem on the Hopkins 155 dataset.
%
% cd to the main folder containing the Hopkins 155 sequences
% add the path to the folder "SSC_motion_face" containing these m-files
%
% avgmissrate1: the n-th element contains the average clustering error for 
% sequences with n motions (using 2F-dimensional data)
% avgmissrate2: the n-th element contains the average clustering error for 
% sequences with n motions (using 4n-dimensional data)
% medmissrate1: the n-th element contains the median clustering error for 
% sequences with n motions (using 2F-dimensional data)
% medmissrate2: the n-th element contains the median clustering error for 
% sequences with n motions (using 4n-dimensional data)
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------

clc, clear all, close all

cd 'C:\Users\Pastoral\Desktop\Hopkins155';
addpath 'C:\Users\Pastoral\Desktop\Copy_SSC_ADMM_v1.1';
addpath(genpath('C:\Users\Pastoral\Desktop\Pursuit\Functions'));
addpath(genpath('C:\Users\Pastoral\Desktop\vmc-master\vmc-master'));


outlier_ratio = 0.5; shuffled_ratio = 0.4;
r = 0; affine = true; outlier = true; rho = 0.7;alpha = 150;
maxNumGroup = 5;
for i = 1:maxNumGroup
    num(i) = 0;
end

d = dir;
for i = 1:length(d)
    if ( (d(i).isdir == 1) && ~strcmp(d(i).name,'.') && ~strcmp(d(i).name,'..') )
        filepath = d(i).name;
        eval(['cd ' filepath]);
        
        f = dir;
        foundValidData = false;
        for j = 1:length(f)
            if ( ~isempty(strfind(f(j).name,'_truth.mat')) )
                ind = j;
                foundValidData = true;
                break
            end
        end
        eval(['load ' f(ind).name]);
        ...;
        disp(['Dataset : ' f(ind).name]);
        cd ..
        
        if (foundValidData)
            n = max(s);
            N = size(x,2);
            F = size(x,3);
            D = 2*F;
            X_raw = reshape(permute(x(1:2,:,:),[1 3 2]),D,N);

            num_shuffled = fix(N*outlier_ratio);
            outliers_ID = randperm(N,num_shuffled);
            inliers_ID = setdiff(1:N,outliers_ID);
            
            xx = shuffle_h(x, shuffled_ratio, outliers_ID);
            X_tilde = reshape(permute(xx(1:2,:,:),[1 3 2]),D,N);
                  
            X_std = X_tilde./vecnorm(X_tilde,2);

% Outliers_detection
affine2 = true ;alpha2 = 10^9; detect_flag = 1;
[Inliers_id, Outliers_id] = Outlier_Detect(X_std, affine2, detect_flag, alpha2);


err_in = length(setdiff(Inliers_id,inliers_ID))/length(Inliers_id);
err_out = length(setdiff(Outliers_id,outliers_ID))/length(Outliers_id);

str = ['err_in = ' num2str(err_in) ', err_out = ' num2str(err_out)];
disp(str);
            

[missrate_in,~,~,ggrps] = SSC(X_tilde(:,Inliers_id),r,affine,alpha,outlier,rho,s(Inliers_id));
str = ['missrate_in = ' num2str(missrate_in)];
disp(str);

           
label_hat = [];
for k = 1:n
    idx = find(ggrps == k);
    idxx{k} = Inliers_id(idx);
    X_sel{1,k} = X_tilde(:,idxx{k}); 
    X_sel_std{1,k} = X_std(:,idxx{k});
    % The following procedure is designed for the final assessment.
    label_hat = [label_hat; idxx{k}, k*ones(length(idxx{k}),1)];
end
Q = sortrows(label_hat);
Label_sel = Q(:,2);

Basis_estimated = zeros(D, 4, n);
Basis_estimated_std = zeros(D, 4, n);
for k = 1:n
        [U,~,~] = svd(X_sel{1,k});
        [UU,~,~] = svd(X_sel_std{1,k});
       Basis_estimated(:,:,k) = U(:,1:4) ;
       Basis_estimated_std(:,:,k) = UU(:,1:4) ;
end
           
if (shuffled_ratio >= 0.6)
    [X_solved,~, Label_unsel] = Solve_all(X_tilde, Outliers_id, Basis_estimated);
else
    [X_solved,~, Label_unsel] = Solve_partial_new(X_tilde, Outliers_id, ...
        Basis_estimated, Basis_estimated_std);
end          
            
IDX = zeros(N,1);
for p = 1:size(X_tilde,2)
    if (ismember(p,Inliers_id) == 1)
        IDX(p,1) = Label_sel(Inliers_id == p);
    elseif (ismember(p,Outliers_id) == 1)
        IDX(p,1) = Label_unsel(Outliers_id == p);
    end
end   

missrate_out = sum(IDX(Outliers_id) ~= s(Outliers_id))/length(s(Outliers_id));
str = ['missrate_out = ' num2str(missrate_out)];
disp(str);
[err_ratio] = Evaluate(X_solved, X_raw);
[err_ratio_in] = Evaluate(X_solved(:,inliers_ID), X_raw(:,inliers_ID));
[err_ratio_out] = Evaluate(X_solved(:,outliers_ID),X_raw(:,outliers_ID));
            

            num(n) = num(n) + 1;
            missrate_in_Tot{n}(num(n)) = missrate_in;
            missrate_out_Tot{n}(num(n)) = missrate_out;
            err_in_Tot{n}(num(n)) = err_in;
            err_out_Tot{n}(num(n)) = err_out;
            err_ratio_in_Tot{n}(num(n)) = err_ratio_in;
            err_ratio_out_Tot{n}(num(n)) = err_ratio_out;
            
            
            eval(['cd ' filepath]);
            save SSC_MS.mat missrate_in missrate_out err_in err_out err_ratio_in err_ratio_out alpha
            cd ..
        end   
    end
end

L = [2 3];
for i = 1:length(L)
    j = L(i);
    avgmissrate_in(j) = mean(missrate_in_Tot{j});
    medmissrate_in(j) = median(missrate_in_Tot{j});
    avgmissrate_out(j) = mean(missrate_out_Tot{j});
    medmissrate_out(j) = median(missrate_out_Tot{j});
    avg_err_in(j) = mean(err_in_Tot{j});
    med_err_in(j) = median(err_in_Tot{j});
    avg_err_out(j) = mean(err_out_Tot{j});
    med_err_out(j) = median(err_out_Tot{j});
    avg_err_ratio_in(j) = mean(err_ratio_in_Tot{j});
    med_err_ratio_in(j) = median(err_ratio_in_Tot{j});
    avg_err_ratio_out(j) = mean(err_ratio_out_Tot{j});
    med_err_ratio_out(j) = median(err_ratio_out_Tot{j});
end
save SSC_MS.mat missrate_in_Tot avgmissrate_in medmissrate_in ...
    missrate_out_Tot avgmissrate_out medmissrate_out ...
   err_in_Tot avg_err_in med_err_in  ...
   err_out_Tot avg_err_out med_err_out ...
   err_ratio_in_Tot avg_err_ratio_in med_err_ratio_in  ...
   err_ratio_out_Tot avg_err_ratio_out med_err_ratio_out ...
   alpha
   