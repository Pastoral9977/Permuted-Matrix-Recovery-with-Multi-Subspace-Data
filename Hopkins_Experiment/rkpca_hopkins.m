%%%%%%%%%%%%%%%%%%%%
% Hopkins for RKPCA 
%%%%%%%%%%%%%%%%%%%%
%% Initialization
clear
close all
addpath(genpath(pwd));
addpath(genpath('..\Functions'));
cd 'Hopkins155';
 

outlier_ratio = 0.4; 
shuffled_ratio = 0.4;
maxNumGroup = 5;
rr = 4;
seed = 42;
involve_gt_bases = true;
num = zeros(1,maxNumGroup);
missrate_in_Tot = cell(1, maxNumGroup);

missrate_out_Tot = cell(1, maxNumGroup);
missrate_out_Tot_gtbases = cell(1, maxNumGroup);

err_ratio_out_Tot_PMSDR = cell(1, maxNumGroup);
err_ratio_out_Tot_gtbases_PMSDR = cell(1, maxNumGroup);
err_ratio_out_Tot_RKPCA = cell(1, maxNumGroup);
err_ratio_out_Tot_RKPCA_PMSDR = cell(1, maxNumGroup);
err_ratio_out_Tot_RKPCA_PMSDR_gt = cell(1, maxNumGroup);

d = dir;


Alpha_set = [10];
for trial = 1:length(Alpha_set)
    c = 0.1;
    filename = sprintf('Hopkins_Report_trial_%d.txt', trial);
    diary(filename);
    dataset_count = 0;
    verbose = false;
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
            rng(seed);
           %% data Information
            n = length(unique(s)); % s-->gt labels
            num(n) = num(n) + 1;
            N = size(x,2); 
            F = size(x,3);
            D = 2*F;
            X_gt = reshape(permute(x(1:2,:,:),[1 3 2]),D,N);
            Vnorms = vecnorm(X_gt);
            
            num_shuffled = fix(N*outlier_ratio);
            outliers_ID = sort(randperm(N,num_shuffled));
            inliers_ID = setdiff(1:N,outliers_ID);
            
            fprintf('--Data Information:--\n')            
            X_tilde = X_gt;
            for q = 1:num_shuffled
                j = outliers_ID(q);
                X_tilde(:,j) = shuffle(X_gt(:,j),shuffled_ratio);
            end

            %% Outlier Detection
            fprintf('--Outlier Detection:--\n')
            alpha = 5;
            lambda = 0.95;
            thres = 1e-3;
            [Inliers_id, Outliers_id] = outlier_detection(X_tilde, alpha, lambda, thres);
%             remark(Inliers_id, Outliers_id, inliers_ID, outliers_ID);
            
            %% RKPCA Subspace Clustering
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
                [X_solved_gtbases, outlier_class_solved_gtbases] = Solve_partial_LSR(X_tilde, outliers_ID, U_gt, false);
                allpoints_class_solved_gtbases = zeros(N,1);
                for p = 1:size(X_tilde,2)
                    if (ismember(p,outliers_ID) == 1)
                        allpoints_class_solved_gtbases(p) = outlier_class_solved_gtbases(outliers_ID == p);
                    end
                end
                missrate_out_gtbases = sum(allpoints_class_solved_gtbases(outliers_ID) ~= s(outliers_ID)) / length(outliers_ID);
                fprintf('\tmissrate_out with ground truth bases = %0.4f\n', missrate_out_gtbases);
                [err_ratio_out_gtbases] = Evaluate(X_solved_gtbases,X_gt);
                fprintf('\terr_ratio_out_gtbases = %.4f (not only outliers part)\n\n', err_ratio_out_gtbases);
                missrate_out_Tot_gtbases{n}(num(n)) = missrate_out_gtbases;
                err_ratio_out_Tot_gtbases_PMSDR{n}(num(n)) = err_ratio_out_gtbases;
            end
             
            %% Matrix Recovery with Reconstructed Bases
            fprintf('--Matrix Recovery with Reconstructed Bases:--\n')
            [X_solved, outlier_class_solved] = Solve_partial_LSR(X_tilde, Outliers_id, U_esti, false);     
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
            misclassified_outlier_num = sum(allpoints_class_solved(detected_true_Outliers_id) ~= s(detected_true_Outliers_id));
            missrate_out = (misclassified_outlier_num) / (length(outliers_ID)-undetected_outlier_num);
            undetected_outlier_num_ratio = undetected_outlier_num / length(outliers_ID);
            fprintf('\tmissrate_out  = %0.4f, where undetected_outlier_num_ratio = %.4f\n', missrate_out, undetected_outlier_num_ratio);
            
            [err_ratio_out_PMSDR] = Evaluate(X_solved,X_gt);
            fprintf('\terr_ratio_out_PMSDR = %.4f\n\n', err_ratio_out_PMSDR);
            
            missrate_in_Tot{n}(num(n)) = missrate_in;
            missrate_out_Tot{n}(num(n)) = missrate_out;
            err_ratio_out_Tot_PMSDR{n}(num(n)) = err_ratio_out_PMSDR;
            
            %% Matrix Recovery with RKPCA
            fprintf('--Matrix Recovery with RKPCA:--\n')
            lambdas = [0.15:0.05:0.6];
            ker.type='rbf';
            ker.par=[0];
            options.p=1;
            options.c=c;
            options.verbose=verbose;
            err_ratio_out_RKPCA = inf;
            display_progress()
            for ii = 1:length(lambdas)
                par = lambdas(ii);
                [X_RKPCA,E_t,f]=RKPCA_PLMAdSS(X_tilde./Vnorms,ker,par,options);
                [err] = Evaluate(X_RKPCA.*Vnorms, X_gt);
                if err < err_ratio_out_RKPCA
                    lambda = par;
                    err_ratio_out_RKPCA = err;
                end
                update_progress(ii, length(lambdas))
            end
            fprintf('\terr_ratio_out_RKPCA = %.4f\n\n', err_ratio_out_RKPCA);
            fprintf('\tchosen lambda = %.1f\n', lambda);
            err_ratio_out_Tot_RKPCA{n}(num(n)) = err_ratio_out_RKPCA;
            
            %% Matrix Recovery with RKPCA PMSDR
            fprintf('--Matrix Recovery with RKPCA PMSDR:--\n')
            X_PMSDR_RKPCA = zeros(size(X_tilde));
            XX_tilde = X_tilde./Vnorms;
            for id = 1:n
                if sum(allpoints_class_solved==id) == 0
                    continue
                end
                XX = XX_tilde(:, allpoints_class_solved==id);
                ker.type='rbf';
                ker.par=[0];
                options.p=1;
                options.c=c;
                options.verbose=verbose;
                [XX_hat,E_t,f]=RKPCA_PLMAdSS(XX,ker,lambda,options);
                X_PMSDR_RKPCA(:, allpoints_class_solved==id) = XX_hat;
            end
            [err_ratio_out_RKPCA_PMSDR] = Evaluate(X_PMSDR_RKPCA.*Vnorms,X_gt);
            fprintf('\terr_ratio_out_RKPCA_PMSDR = %.4f\n\n', err_ratio_out_RKPCA_PMSDR);
            err_ratio_out_Tot_RKPCA_PMSDR{n}(num(n)) = err_ratio_out_RKPCA_PMSDR;
            
            %% Matrix Recovery with RKPCA PMSDR gt
            fprintf('--Matrix Recovery with RKPCA PMSDR GT:--\n')
            X_PMSDR_RKPCA_gt = zeros(size(X_tilde));
            XX_tilde = X_tilde./Vnorms;
            for id = 1:n
                if sum(s==id) == 0
                    continue
                end
                XX = XX_tilde(:, s==id);
                ker.type='rbf';
                ker.par=[0];
                options.p=1;
                options.c=c;
                options.verbose=verbose;
                [XX_hat,E_t,f]=RKPCA_PLMAdSS(XX,ker,lambda,options);
                X_PMSDR_RKPCA_gt(:, s==id) = XX_hat;
            end
            [err_ratio_out_RKPCA_PMSDR_gt] = Evaluate(X_PMSDR_RKPCA_gt.*Vnorms,X_gt);
            fprintf('\terr_ratio_out_RKPCA_PMSDR_gt = %.4f\n\n', err_ratio_out_RKPCA_PMSDR_gt);
            err_ratio_out_Tot_RKPCA_PMSDR_gt{n}(num(n)) = err_ratio_out_RKPCA_PMSDR_gt;
            
       
            eval(['cd ' filepath]);
            save RKPCA_MS.mat missrate_in missrate_out err_ratio_out_gtbases err_ratio_out_PMSDR err_ratio_out_RKPCA err_ratio_out_RKPCA_PMSDR alpha
            cd ..
        end   
    end
end


L = [1 2 3 4 5];
avgmissrate_in = zeros(1,length(L));
medmissrate_in = zeros(1,length(L));
avgmissrate_out = zeros(1,length(L));
medmissrate_out = zeros(1,length(L));
avg_err_ratio_out_PMSDR = zeros(1,length(L));
med_err_ratio_out_PMSDR = zeros(1,length(L));
avg_err_ratio_out_RKPCA = zeros(1,length(L));
med_err_ratio_out_RKPCA = zeros(1,length(L));
avg_err_ratio_out_RKPCA_PMSDR = zeros(1,length(L));
med_err_ratio_out_RKPCA_PMSDR = zeros(1,length(L));
avg_err_ratio_out_RKPCA_PMSDR_gt = zeros(1,length(L));
med_err_ratio_out_RKPCA_PMSDR_gt = zeros(1,length(L));

for i = 1:length(L)
    j = L(i);
    avgmissrate_in(j)   = mean(missrate_in_Tot{j});
    medmissrate_in(j)   = median(missrate_in_Tot{j});
    avgmissrate_out(j)  = mean(missrate_out_Tot{j});
    medmissrate_out(j)  = median(missrate_out_Tot{j});
    avg_err_ratio_out_PMSDR(j) = mean(err_ratio_out_Tot_PMSDR{j});
    med_err_ratio_out_PMSDR(j) = median(err_ratio_out_Tot_PMSDR{j});
    avg_err_ratio_out_RKPCA(j) = mean(err_ratio_out_Tot_RKPCA{j});
    med_err_ratio_out_RKPCA(j) = median(err_ratio_out_Tot_RKPCA{j});
    avg_err_ratio_out_RKPCA_PMSDR(j) = mean(err_ratio_out_Tot_RKPCA_PMSDR{j});
    med_err_ratio_out_RKPCA_PMSDR(j) = median(err_ratio_out_Tot_RKPCA_PMSDR{j});
    avg_err_ratio_out_RKPCA_PMSDR_gt(j) = mean(err_ratio_out_Tot_RKPCA_PMSDR_gt{j});
    med_err_ratio_out_RKPCA_PMSDR_gt(j) = median(err_ratio_out_Tot_RKPCA_PMSDR_gt{j});
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
        avg_err_ratio_out_gtbases(j) = mean(err_ratio_out_Tot_gtbases_PMSDR{j});
        med_err_ratio_out_gtbases(j) = median(err_ratio_out_Tot_gtbases_PMSDR{j});
    end  
end

allmissrate_in     = [missrate_in_Tot{2}, missrate_in_Tot{3}, missrate_in_Tot{5}];
allmissrate_out     = [missrate_out_Tot{2}, missrate_out_Tot{3}, missrate_out_Tot{5}];
allerr_ratio_out_PMSDR    = [err_ratio_out_Tot_PMSDR{2}, err_ratio_out_Tot_PMSDR{3}, err_ratio_out_Tot_PMSDR{5}];
allerr_ratio_out_RKPCA    = [err_ratio_out_Tot_RKPCA{2}, err_ratio_out_Tot_RKPCA{3}, err_ratio_out_Tot_RKPCA{5}];
allerr_ratio_out_RKPCA_PMSDR    = [err_ratio_out_Tot_RKPCA_PMSDR{2}, err_ratio_out_Tot_RKPCA_PMSDR{3}, err_ratio_out_Tot_RKPCA_PMSDR{5}];
med_all_missrate_in    = median(allmissrate_in);
med_all_missrate_out    = median(allmissrate_out);
med_all_err_ratio_out_PMSDR   = median(allerr_ratio_out_PMSDR);
med_all_err_ratio_out_RKPCA   = median(allerr_ratio_out_RKPCA);
med_all_err_ratio_out_RKPCA_PMSDR   = median(allerr_ratio_out_RKPCA_PMSDR);
avg_all_missrate_in    = mean(allmissrate_in);
avg_all_missrate_out    = mean(allmissrate_out);
avg_all_err_ratio_out_PMSDR   = mean(allerr_ratio_out_PMSDR);
avg_all_err_ratio_out_RKPCA   = mean(allerr_ratio_out_RKPCA);
avg_all_err_ratio_out_RKPCA_PMSDR   = mean(allerr_ratio_out_RKPCA_PMSDR);

fprintf([repmat('=', 1, 60), '\n']);
fprintf('Summary with SSC alpha = %d\n', alpha);
fprintf([repmat('=', 1, 60), '\n']);

disp(['Median  missrate_in   for 2,3,5 group  =   ' num2str([medmissrate_in(2), medmissrate_in(3), medmissrate_in(5)], '%.3f, %.3f, %.3f')]);
disp(['Average missrate_in   for 2,3,5 group  =   ' num2str([avgmissrate_in(2), avgmissrate_in(3), avgmissrate_in(5)], '%.3f, %.3f, %.3f')]);

disp(['Median  missrate_out  for 2,3,5 group  =   ' num2str([medmissrate_out(2), medmissrate_out(3), medmissrate_out(5)], '%.3f, %.3f, %.3f')]);
disp(['Average missrate_out  for 2,3,5 group  =   ' num2str([avgmissrate_out(2), avgmissrate_out(3), avgmissrate_out(5)], '%.3f, %.3f, %.3f')]);


disp(['Median  missrate_out_gtbases  for 2,3,5 group  =   ' num2str([medmissrate_out_gtbases(2), medmissrate_out_gtbases(3), medmissrate_out_gtbases(5)], '%.3f, %.3f, %.3f')]);
disp(['Average missrate_out_gtbases  for 2,3,5 group  =   ' num2str([avgmissrate_out_gtbases(2), avgmissrate_out_gtbases(3), avgmissrate_out_gtbases(5)], '%.3f, %.3f, %.3f')]);

disp(repmat('--',1,20))

disp(['Median  err_ratio_out_PMSDR_gtbases for 2,3,5 group  =   ' num2str([med_err_ratio_out_gtbases(2), med_err_ratio_out_gtbases(3), med_err_ratio_out_gtbases(5)], '%.3f, %.3f, %.3f')]);
disp(['Average err_ratio_out_PMSDR_gtbases for 2,3,5 group  =   ' num2str([avg_err_ratio_out_gtbases(2), avg_err_ratio_out_gtbases(3), avg_err_ratio_out_gtbases(5)], '%.3f, %.3f, %.3f')]);

disp(['Median  err_ratio_out_PMSDR for 2,3,5 group  =   ' num2str([med_err_ratio_out_PMSDR(2), med_err_ratio_out_PMSDR(3), med_err_ratio_out_PMSDR(5)], '%.3f, %.3f, %.3f')]);
disp(['Average err_ratio_out_PMSDR for 2,3,5 group  =   ' num2str([avg_err_ratio_out_PMSDR(2), avg_err_ratio_out_PMSDR(3), avg_err_ratio_out_PMSDR(5)], '%.3f, %.3f, %.3f')]);

disp(repmat('--',1,20))

disp(['Median  err_ratio_out_RKPCA for 2,3,5 group  =   ' num2str([med_err_ratio_out_RKPCA(2), med_err_ratio_out_RKPCA(3), med_err_ratio_out_RKPCA(5)], '%.3f, %.3f, %.3f')]);
disp(['Average err_ratio_out_RKPCA for 2,3,5 group  =   ' num2str([avg_err_ratio_out_RKPCA(2), avg_err_ratio_out_RKPCA(3), avg_err_ratio_out_RKPCA(5)], '%.3f, %.3f, %.3f')]);

disp(['Median  err_ratio_out_RKPCA_PMSDR for 2,3,5 group  =   ' num2str([med_err_ratio_out_RKPCA_PMSDR(2), med_err_ratio_out_RKPCA_PMSDR(3), med_err_ratio_out_RKPCA_PMSDR(5)], '%.3f, %.3f, %.3f')]);
disp(['Average err_ratio_out_RKPCA_PMSDR for 2,3,5 group  =   ' num2str([avg_err_ratio_out_RKPCA_PMSDR(2), avg_err_ratio_out_RKPCA_PMSDR(3), avg_err_ratio_out_RKPCA_PMSDR(5)], '%.3f, %.3f, %.3f')]);

disp(['Median  err_ratio_out_RKPCA_PMSDR_gt for 2,3,5 group  =   ' ...
    num2str([med_err_ratio_out_RKPCA_PMSDR_gt(2), med_err_ratio_out_RKPCA_PMSDR_gt(3), med_err_ratio_out_RKPCA_PMSDR_gt(5)], '%.3f, %.3f, %.3f')]);
disp(['Average err_ratio_out_RKPCA_PMSDR_gt for 2,3,5 group  =   ' ...
    num2str([avg_err_ratio_out_RKPCA_PMSDR_gt(2), avg_err_ratio_out_RKPCA_PMSDR_gt(3), avg_err_ratio_out_RKPCA_PMSDR_gt(5)], '%.3f, %.3f, %.3f')]);

disp(repmat('--',1,20))

end
