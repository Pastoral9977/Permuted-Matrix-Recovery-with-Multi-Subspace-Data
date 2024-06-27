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
%% Initialization
clear
close all
addpath(genpath('..\Copy_SSC_ADMM_v1.1'));
addpath(genpath('..\Functions'));
cd 'C:\Users\Pastoral\Desktop\Hopkins155';

outlier_ratio = 0.4; 
shuffled_ratio = 0.4;
maxNumGroup = 5;
rr = 4;
involve_gt_bases = true;
num = zeros(1,maxNumGroup);
missrate_in_Tot = cell(1, maxNumGroup);
missrate_out_Tot = cell(1, maxNumGroup);
err_ratio_out_Tot = cell(1, maxNumGroup);
missrate_out_Tot_gtbases = cell(1, maxNumGroup);
err_ratio_out_Tot_gtbases = cell(1, maxNumGroup);

d = dir;
dataset_count = 0;
for i = 1:length(d)
    
    %% check the data dir
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
        dataset_count = dataset_count + 1;
        eval(['load ' f(ind).name]);
        fprintf([repmat('=', 1, 80), '\n']);
        fprintf(['Dataset(' num2str(dataset_count) ') : ' f(ind).name '\n']);
        fprintf([repmat('=', 1, 80), '\n']);
        cd ..
        
        if (foundValidData)
            rng(42);
           %% data Information
            n = length(unique(s)); % s-->gt labels
            num(n) = num(n) + 1;
            N = size(x,2); 
            F = size(x,3);
            D = 2*F;
            X_gt = reshape(permute(x(1:2,:,:),[1 3 2]),D,N);
            
            num_shuffled = fix(N*outlier_ratio);
            outliers_ID = sort(randperm(N,num_shuffled));
            inliers_ID = setdiff(1:N,outliers_ID);
            
            fprintf('--Data Information:--\n')
            fprintf(...
                    ['\tnum_groups = '                  num2str(n)                              ...
                    '\n\tnumber of total samples = '    num2str(size(X_gt,2))                   ...
                    '\n\tdimension of each sample = '   num2str(size(X_gt,1))                   ...
                    '\n\tnumber of outliers = '         num2str(num_shuffled)                   ...
                    '\n\tshuffled_ratio = '             num2str(shuffled_ratio) '\n\n']         ...
                    )     
            
            X_tilde = X_gt;
            for q = 1:num_shuffled
                j = outliers_ID(q);
                X_tilde(:,j) = shuffle(X_gt(:,j),shuffled_ratio);
            end

% %             Outliers_detection
%             affine2 = true ;alpha2 = 1e10; detect_flag = 1;validation = -1;
%             [Inliers_id, Outliers_id] = Outlier_Detect_v2(cnormalize(X_tilde),detect_flag,alpha2,validation,affine2);
            %% Outlier Detection
            fprintf('--Outlier Detection:--\n')
            alpha = 5;
            lambda = 0.95;
            thres = 1e-3;
            [Inliers_id, Outliers_id] = Outlier_Detect_v3(X_tilde, alpha, lambda, thres);
            remark(Inliers_id, Outliers_id, inliers_ID, outliers_ID);
            
            %% SSC Subspace Clustering
            fprintf('\n--Subspace Clustering:--\n')
            r = 0; 
            affine = true; 
            outlier = true; 
            rho = 0.7; 
            alpha = 500;
            fprintf('\talpha = %d, affine = %d, outlier = %d, rho = %f, r = %f \n',alpha, affine, outlier, rho, r)
            [missrate_in, ggrps] = SSC(X_tilde(:,Inliers_id),r,affine,alpha,outlier,rho,s(Inliers_id));
            fprintf(['\tmissrate_in = ' num2str(missrate_in) '\n\n']);
           
            %% Reconstructed Bases Estimation
            label_hat = zeros(length(Inliers_id), 2);
            X_selected_inlier_groups = cell(1, n);
            index = 1;
            for k = 1:n
                kth_group_ids = Inliers_id(ggrps == k);
                X_selected_inlier_groups{k} = X_tilde(:,kth_group_ids); 
                % The following procedure is designed for the final assessment.
                num_kth_group_ids = length(kth_group_ids);
                label_hat(index:index+num_kth_group_ids-1, :) = [reshape(kth_group_ids, [num_kth_group_ids, 1]), k*ones(num_kth_group_ids, 1)];
                index = index + num_kth_group_ids;
            end
            Q = sortrows(label_hat);
            inlier_class_solved = Q(:,2);
            
            U_esti = zeros(D, rr, n);
            for k = 1:n
                [U,~,~] = svd(X_selected_inlier_groups{k});
                U_esti(:,:,k) = U(:,1:rr) ;
            end
             
            %% Matrix Recovery with gt bases
            fprintf('--Matrix Recovery with gt bases:--\n')
            if involve_gt_bases
                U_gt = ground_truth_bases_construction(X_gt, s, rr);
                [X_solved_gtbases, outlier_class_solved_gtbases] = Solve_partial_LSR(X_tilde, outliers_ID, U_gt);
                allpoints_class_solved_gtbases = zeros(N,1);
                for p = 1:size(X_tilde,2)
                    if (ismember(p,outliers_ID) == 1)
                        allpoints_class_solved_gtbases(p) = outlier_class_solved_gtbases(outliers_ID == p);
                    end
                end
                missrate_out_gtbases = sum(allpoints_class_solved_gtbases(outliers_ID) ~= s(outliers_ID)) / length(outliers_ID);
                fprintf('\tmissrate_out with ground truth bases = %0.4f\n', missrate_out_gtbases);
                [err_ratio_out_gtbases] = Evaluate(X_solved_gtbases(:,outliers_ID),X_gt(:,outliers_ID));
                fprintf('\terr_ratio_out = %.4f\n\n', err_ratio_out_gtbases);
                missrate_out_Tot_gtbases{n}(num(n)) = missrate_out_gtbases;
                err_ratio_out_Tot_gtbases{n}(num(n)) = err_ratio_out_gtbases;
            end
             
            %% Matrix Recovery with reconstructed bases
            fprintf('--Matrix Recovery with Reconstructed Bases:--\n')
            [X_solved, outlier_class_solved] = Solve_partial_LSR(X_tilde, Outliers_id, U_esti);     
            allpoints_class_solved = zeros(N,1);
            for p = 1:size(X_tilde,2)
                if (ismember(p,Inliers_id) == 1)
                    allpoints_class_solved(p,1) = inlier_class_solved(Inliers_id == p);
                elseif (ismember(p,Outliers_id) == 1)
                    allpoints_class_solved(p,1) = outlier_class_solved(Outliers_id == p);
                end
            end   
            
            undetected_outlier_num = length(setdiff(outliers_ID, Outliers_id));
            detected_true_Outliers_id = setdiff(Outliers_id, inliers_ID);
            missclassified_outlier_num = sum(allpoints_class_solved(detected_true_Outliers_id) ~= s(detected_true_Outliers_id));
            missrate_out = (missclassified_outlier_num) / (length(outliers_ID)-undetected_outlier_num);
            undetected_outlier_num_ratio = undetected_outlier_num / length(outliers_ID);
            fprintf('\tmissrate_out  = %0.4f, where undetected_outlier_num_ratio = %.4f\n', missrate_out, undetected_outlier_num_ratio);
            
            [err_ratio_out] = Evaluate(X_solved(:,outliers_ID),X_gt(:,outliers_ID));
            fprintf('\terr_ratio_out = %.4f\n\n', err_ratio_out);

            missrate_in_Tot{n}(num(n)) = missrate_in;
            missrate_out_Tot{n}(num(n)) = missrate_out;
            err_ratio_out_Tot{n}(num(n)) = err_ratio_out;
            
            eval(['cd ' filepath]);
            save SSC_MS.mat missrate_in missrate_out err_ratio_out alpha
            cd ..
        end   
    end
end

L = [1 2 3 4 5];
avgmissrate_in = zeros(1,length(L));
medmissrate_in = zeros(1,length(L));
avgmissrate_out = zeros(1,length(L));
medmissrate_out = zeros(1,length(L));
avg_err_ratio_out = zeros(1,length(L));
med_err_ratio_out = zeros(1,length(L));
for i = 1:length(L)
    j = L(i);
    avgmissrate_in(j)   = mean(missrate_in_Tot{j});
    medmissrate_in(j)   = median(missrate_in_Tot{j});
    avgmissrate_out(j)  = mean(missrate_out_Tot{j});
    medmissrate_out(j)  = median(missrate_out_Tot{j});
    avg_err_ratio_out(j) = mean(err_ratio_out_Tot{j});
    med_err_ratio_out(j) = median(err_ratio_out_Tot{j});
end

if involve_gt_bases
    avgmissrate_out_gtbases = zeros(1,length(L));
    medmissrate_out_gtbases = zeros(1,length(L));
    avg_err_ratio_out_gtbases = zeros(1,length(L));
    med_err_ratio_out_gtbases = zeros(1,length(L));
    for i = 1:length(L)
        j = L(i);
        avgmissrate_out_gtbases(j)  = mean(missrate_out_Tot_gtbases{j});
        medmissrate_out_gtbases(j)  = median(missrate_out_Tot_gtbases{j});
        avg_err_ratio_out_gtbases(j) = mean(err_ratio_out_Tot_gtbases{j});
        med_err_ratio_out_gtbases(j) = median(err_ratio_out_Tot_gtbases{j});
    end  
end

allmissrate_in     = [missrate_in_Tot{2}, missrate_in_Tot{3}, missrate_in_Tot{5}];
allmissrate_out     = [missrate_out_Tot{2}, missrate_out_Tot{3}, missrate_out_Tot{5}];
allerr_ratio_out    = [err_ratio_out_Tot{2}, err_ratio_out_Tot{3}, err_ratio_out_Tot{5}];
med_all_missrate_in    = median(allmissrate_in);
med_all_missrate_out    = median(allmissrate_out);
med_all_err_ratio_out   = median(allerr_ratio_out);
avg_all_missrate_in    = mean(allmissrate_in);
avg_all_missrate_out    = mean(allmissrate_out);
avg_all_err_ratio_out   = mean(allerr_ratio_out);
save SSC_MS.mat         missrate_in_Tot         avgmissrate_in          medmissrate_in          ...
     missrate_out_Tot   avgmissrate_out         medmissrate_out         ...
     alpha              avg_err_ratio_out       med_err_ratio_out       err_ratio_out_Tot
 
if involve_gt_bases
    allmissrate_out_gtbases     = [missrate_out_Tot_gtbases{2}, missrate_out_Tot_gtbases{3}, missrate_out_Tot_gtbases{5}];
    allerr_ratio_out_gtbases    = [err_ratio_out_Tot_gtbases{2}, err_ratio_out_Tot_gtbases{3}, err_ratio_out_Tot_gtbases{5}];
    med_all_missrate_out_gtbases    = median(allmissrate_out_gtbases);
    med_all_err_ratio_out_gtbases   = median(allerr_ratio_out_gtbases);
    avg_all_missrate_out_gtbases    = mean(allmissrate_out_gtbases);
    avg_all_err_ratio_out_gtbases   = mean(allerr_ratio_out_gtbases);
    save SSC_MS.mat         missrate_in_Tot         avgmissrate_in          medmissrate_in          ...
         missrate_out_Tot   avgmissrate_out         medmissrate_out         missrate_out_Tot_gtbases...
         alpha              avg_err_ratio_out       med_err_ratio_out       err_ratio_out_Tot       ...
         err_ratio_out_Tot_gtbases                  avg_err_ratio_out_gtbases                       ...
         med_err_ratio_out_gtbases                  medmissrate_out_gtbases                         ...
         avgmissrate_out_gtbases
end

 
fprintf([repmat('=', 1, 60), '\n']);
fprintf('Summary with SSC alpha = %d \n', alpha);
fprintf([repmat('=', 1, 60), '\n']);

disp(['Median  missrate_in   for 2,3,5 group  =   ' num2str([medmissrate_in(2), medmissrate_in(3), medmissrate_in(5)], '%.3f, %.3f, %.3f')]);
disp(['Average missrate_in   for 2,3,5 group  =   ' num2str([avgmissrate_in(2), avgmissrate_in(3), avgmissrate_in(5)], '%.3f, %.3f, %.3f')]);

disp(['Median  missrate_out  for 2,3,5 group  =   ' num2str([medmissrate_out(2), medmissrate_out(3), medmissrate_out(5)], '%.3f, %.3f, %.3f')]);
disp(['Average missrate_out  for 2,3,5 group  =   ' num2str([avgmissrate_out(2), avgmissrate_out(3), avgmissrate_out(5)], '%.3f, %.3f, %.3f')]);

disp(['Median  err_ratio_out for 2,3,5 group  =   ' num2str([med_err_ratio_out(2), med_err_ratio_out(3), med_err_ratio_out(5)], '%.3f, %.3f, %.3f')]);
disp(['Average err_ratio_out for 2,3,5 group  =   ' num2str([avg_err_ratio_out(2), avg_err_ratio_out(3), avg_err_ratio_out(5)], '%.3f, %.3f, %.3f')]);

if involve_gt_bases
    disp(['Median  missrate_out_gtbases  for 2,3,5 group  =   ' num2str([medmissrate_out_gtbases(2), medmissrate_out_gtbases(3), medmissrate_out_gtbases(5)], '%.3f, %.3f, %.3f')]);
    disp(['Average missrate_out_gtbases  for 2,3,5 group  =   ' num2str([avgmissrate_out_gtbases(2), avgmissrate_out_gtbases(3), avgmissrate_out_gtbases(5)], '%.3f, %.3f, %.3f')]);

    disp(['Median  err_ratio_out_gtbases for 2,3,5 group  =   ' num2str([med_err_ratio_out_gtbases(2), med_err_ratio_out_gtbases(3), med_err_ratio_out_gtbases(5)], '%.3f, %.3f, %.3f')]);
    disp(['Average err_ratio_out_gtbases for 2,3,5 group  =   ' num2str([avg_err_ratio_out_gtbases(2), avg_err_ratio_out_gtbases(3), avg_err_ratio_out_gtbases(5)], '%.3f, %.3f, %.3f')]);
end
 