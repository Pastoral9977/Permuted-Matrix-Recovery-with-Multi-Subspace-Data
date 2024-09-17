%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Copyright @ Liangqi Xie, 2012
%--------------------------------------------------------------------------
%% Initialization
clear
close all
addpath(genpath('..\Functions'));
addpath(genpath(pwd));
diary('MRUC_Hopkins.txt')
cd 'Hopkins155';

outlier_ratio = 0.4; 
shuffled_ratio = 0.4;
maxNumGroup = 5;
rr1 = 4;
rr2 = 2;
num = zeros(1,maxNumGroup);
err_ratio_Tot = cell(1, maxNumGroup);
err_ratio_PMSDR_Tot = cell(1, maxNumGroup);
err_ratio_PMSDR_gt_Tot = cell(1, maxNumGroup);
parts = 2;

d = dir;
dataset_count = 0;
SEED = randi(100000);
% SEED = 6429;
for i = 1:length(d)
    seed = SEED*i;
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
        fprintf(['Dataset(' num2str(dataset_count) ') : ' f(ind).name ', parts:%d, seed:%d   [i:%d]\n'],parts,seed,i);
        fprintf([repmat('=', 1, 80), '\n']);
        cd ..
        
        if (foundValidData)
           %% data Information
            rng(seed);
            s(s==0)=1;
            num_groups = length(unique(s));                                 % s-->gt labels
            num(num_groups) = num(num_groups) + 1;
            N = size(x,2); 
            F = size(x,3);
            D = 2*F;
            X_raw = reshape(permute(x(1:2,:,:),[1 3 2]),D,N);
%             X_raw = X_raw./vecnorm(X_raw);
            X_gt = zeros(size(X_raw));
            idx = 1;
            ss=[];
            for jj=1:num_groups
                X_gt(:,idx:idx+sum(s==jj)-1)=X_raw(:,s==jj);
                ss=[ss;ones(sum(s==jj),1)*jj];
                idx=idx+sum(s==jj);
            end
            s=ss;
            clear ss X_raw idx
            
            fprintf('--Data Information:--\n')
            fprintf(...
                    ['\tnum_groups = '                  num2str(num_groups)                              ...
                    '\n\tnumber of total samples = '    num2str(size(X_gt,2))                   ...
                    '\n\tdimension of each sample = '   num2str(size(X_gt,1))                   ...
                    '\n\toutlier_ratio = '              num2str(outlier_ratio)                   ...
                    '\n\tshuffled_ratio = '             num2str(shuffled_ratio) '\n\n']         ...
                    )     
            
            
            %% data preparation
            [DATA_SUB,DATA_TOTAL,perm_labels,X_tilde,outliers_ID,seq_gt] = prepare_mruc_hopkins_data(X_gt,parts,s,shuffled_ratio,outlier_ratio,seed);
            inliers_ID = setdiff(1:N,outliers_ID);
            [DATA_SUB_solved, seq] = PMSDR_MRUC_data_handling(DATA_TOTAL, num_groups, rr1, s(inliers_ID), s(outliers_ID));
            
            % ensure seq_gt is right
            X_test = [];
            for jj = 1:num_groups
                X_test = [X_test,DATA_SUB.B1s{jj}];
                for kk = 1:parts
                    X_test = [X_test,DATA_SUB.Bns{jj}{kk}];
                end
            end
            if norm(X_test-X_tilde(:,seq_gt))~=0
                error('Something Wrong with seq_gt!');
            end
            
            % ensure seq is right
            X_new_tilde = X_tilde(:,[inliers_ID,outliers_ID]);
            X_test = [];
            for jj = 1:num_groups
                X_test = [X_test,DATA_SUB_solved.B1_SUB_solved{jj}];
                for kk = 1:parts
                    X_test = [X_test,DATA_SUB_solved.Bn_SUB_solved{jj}{kk}];
                end
            end
            if norm(X_test-X_new_tilde(:,seq))~=0
                error('Something Wrong with seq!');
            end
            %% MRUC Recovery
            eps_init=0.05;
            eps_min=0.000001;
            verbose=true;
            max_out_iter=200;
            max_in_iter=1000;
            omega=0.5;
            beta=100;
            gd = 1;
            X_new_gt = X_gt(:,[inliers_ID,outliers_ID]);
            lr = 0;
            for init_rank=rr2
                fprintf('initial_rank = %d\n',rr2)
                tol = 10;
                %% total
                fprintf('-----------------MRUC-----------------\n')
                use_B1=DATA_TOTAL.B1_TOTAL;
                use_Bn=DATA_TOTAL.Bn_TOTAL;
                use_A1=DATA_TOTAL.A1_TOTAL;
                use_An=DATA_TOTAL.An_TOTAL;
                use_A1_row_ind=DATA_TOTAL.A1_row_ind_TOTAL;
                use_An_row_ind=DATA_TOTAL.An_row_ind_TOTAL;
                use_test_ind=DATA_TOTAL.test_ind_TOTAL;
                use_test_data=DATA_TOTAL.test_label_TOTAL;
                [best_error,best_P,best_B,final_cert,history,result]=CD_complete_hopkins(use_B1,use_Bn,use_A1,use_An,use_A1_row_ind,use_An_row_ind,init_rank, ... 
                    eps_init,verbose,max_out_iter,max_in_iter,use_test_ind,use_test_data,tol,omega,beta,eps_min,perm_labels,lr,gd,X_new_gt);
                Err = Evaluate(best_B,X_new_gt);
                fprintf('---------------------------------------------------------\n')
                fprintf('\tMRUC:          \tError = %.4f\n',Err)
                fprintf('---------------------------------------------------------\n')
                
                %% sub
                X_hat = [];
                fprintf('-----------------PMSDR MRUC-----------------\n')
                for jj = 1:num_groups
                    use_B1 = DATA_SUB_solved.B1_SUB_solved{jj};
                    use_Bn = DATA_SUB_solved.Bn_SUB_solved{jj};
                    use_A1 = DATA_SUB_solved.A1_SUB_solved{jj};
                    use_An = DATA_SUB_solved.An_SUB_solved{jj};
                    use_A1_row_ind = DATA_SUB_solved.A1_row_ind_SUB_solved{jj};
                    use_An_row_ind = DATA_SUB_solved.An_row_ind_SUB_solved{jj};
                    use_test_ind=DATA_SUB_solved.test_inds_solved{jj};
                    use_test_data=DATA_SUB_solved.test_labels_solved{jj};
                    
                    [best_error,best_P,best_B,final_cert,history,result]=CD_complete_hopkins(use_B1,use_Bn,use_A1,use_An,use_A1_row_ind,use_An_row_ind,init_rank,... 
                    eps_init,verbose,max_out_iter,max_in_iter,use_test_ind,use_test_data,tol,omega,beta,eps_min,perm_labels,lr,gd);
                    X_hat=[X_hat,best_B];
                    fprintf('\n')
                end
                Err_PMSDR = Evaluate(X_hat,X_new_gt(:,seq));
                fprintf('---------------------------------------------------------\n')
                fprintf('\tPMSDR MRUC:    \tError = %.4f\n',Err_PMSDR)
                fprintf('---------------------------------------------------------\n')
                
                %% gt
                X_hat = [];
                fprintf('-----------------PMSDR MRUC gt-----------------\n')
                for jj = 1:num_groups
                    use_B1 = DATA_SUB.B1s{jj};
                    use_Bn = DATA_SUB.Bns{jj};
                    use_A1 = DATA_SUB.A1s{jj};
                    use_An = DATA_SUB.Ans{jj};
                    use_A1_row_ind = DATA_SUB.A1_row_inds{jj};
                    use_An_row_ind = DATA_SUB.An_row_inds{jj};
                    use_test_ind=DATA_SUB.test_inds{jj};
                    use_test_data=DATA_SUB.test_labels{jj};
                    [best_error,best_P,best_B,final_cert,history,result]=CD_complete_hopkins(use_B1,use_Bn,use_A1,use_An,use_A1_row_ind,use_An_row_ind,init_rank,... 
                    eps_init,verbose,max_out_iter,max_in_iter,use_test_ind,use_test_data,tol,omega,beta,eps_min,perm_labels,lr,gd);
                    X_hat=[X_hat,best_B];
                    fprintf('\n')
                end
                Err_PMSDR_gt = Evaluate(X_hat,X_gt(:,seq_gt));
                fprintf('---------------------------------------------------------\n')
                fprintf('\tPMSDR MRUC gt: \tError = %.4f\n',Err_PMSDR_gt)
                fprintf('---------------------------------------------------------\n')
            end
            
            err_ratio_Tot{num_groups}(num(num_groups)) = Err;
            err_ratio_PMSDR_Tot{num_groups}(num(num_groups)) = Err_PMSDR;
            err_ratio_PMSDR_gt_Tot{num_groups}(num(num_groups)) = Err_PMSDR_gt;
            
            eval(['cd ' filepath]);
            cd ..
        end   
    end
end

L = [1 2 3 4 5];
avg_err_ratio = zeros(1,length(L));
med_err_ratio = zeros(1,length(L));

avg_err_ratio_PMSDR = zeros(1,length(L));
med_err_ratio_PMSDR = zeros(1,length(L));

avg_err_ratio_PMSDR_gt = zeros(1,length(L));
med_err_ratio_PMSDR_gt = zeros(1,length(L));

for i = 1:length(L)
    j = L(i);
    avg_err_ratio(j) = mean(err_ratio_Tot{j});
    med_err_ratio(j) = median(err_ratio_Tot{j});
    
    avg_err_ratio_PMSDR(j) = mean(err_ratio_PMSDR_Tot{j});
    med_err_ratio_PMSDR(j) = median(err_ratio_PMSDR_Tot{j});
    
    avg_err_ratio_PMSDR_gt(j) = mean(err_ratio_PMSDR_gt_Tot{j});
    med_err_ratio_PMSDR_gt(j) = median(err_ratio_PMSDR_gt_Tot{j});
end

save('Results_Hopkins_MRUC')

fprintf([repmat('=', 1, 60), '\n']);
fprintf('Summary for MRUC Recovery \n');
fprintf([repmat('=', 1, 60), '\n']);

disp(['Median  err_ratio for 2,3,5 group  =   ' num2str([med_err_ratio(2), med_err_ratio(3), med_err_ratio(5)], '%.3f, %.3f, %.3f')]);
disp(['Average err_ratio for 2,3,5 group  =   ' num2str([avg_err_ratio(2), avg_err_ratio(3), avg_err_ratio(5)], '%.3f, %.3f, %.3f')]);
disp(['Median  err_ratio_PMSDR for 2,3,5 group  =   ' num2str([med_err_ratio_PMSDR(2), med_err_ratio_PMSDR(3), med_err_ratio_PMSDR(5)], '%.3f, %.3f, %.3f')]);
disp(['Average err_ratio_PMSDR for 2,3,5 group  =   ' num2str([avg_err_ratio_PMSDR(2), avg_err_ratio_PMSDR(3), avg_err_ratio_PMSDR(5)], '%.3f, %.3f, %.3f')]);
disp(['Median  err_ratio_PMSDR_gt for 2,3,5 group  =   ' num2str([med_err_ratio_PMSDR_gt(2), med_err_ratio_PMSDR_gt(3), med_err_ratio_PMSDR_gt(5)], '%.3f, %.3f, %.3f')]);
disp(['Average err_ratio_PMSDR_gt for 2,3,5 group  =   ' num2str([avg_err_ratio_PMSDR_gt(2), avg_err_ratio_PMSDR_gt(3), avg_err_ratio_PMSDR_gt(5)], '%.3f, %.3f, %.3f')]);



 