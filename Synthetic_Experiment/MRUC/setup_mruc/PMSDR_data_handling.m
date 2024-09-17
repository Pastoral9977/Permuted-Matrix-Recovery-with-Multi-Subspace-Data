function [B1_SUB_solved, A1_SUB_solved,A1_row_ind_SUB_solved,An_SUB_solved,An_row_ind_SUB_solved] ...
                                    = PMSDR_data_handling(B1_TOTAL, A1_TOTAL, An_TOTAL, A1_row_ind_TOTAL, An_row_ind_TOTAL, num_groups, rr)
                                
        X = B1_TOTAL;
        [D, parts] = size(An_TOTAL);
        s = [];
        basic_count = floor(size(B1_TOTAL, 2) / num_groups);
        remainder = mod(size(B1_TOTAL, 2), num_groups);
        for ii = 1:num_groups
            count = basic_count + (ii <= remainder);
            s = [s, ii * ones(1, count)];
        end
        
        r = 0;
        affine = false;
        outlier = false;
        alpha = 800;
        rho = 1;
        broadcast = false;

        [missrate_in,ggrps] = SSC(X,r,affine,alpha,outlier,rho,s,broadcast);
        fprintf('\tmissrate_in (Subspace Clustering) : %.4f\n', missrate_in)
        Bases = zeros(D, rr, num_groups);
        B1_SUB_solved = cell(1, num_groups);
        A1_SUB_solved = cell(1, num_groups);
        A1_row_ind_SUB_solved = cell(1, num_groups);
        An_SUB_solved = cell(1, num_groups);
        An_row_ind_SUB_solved = cell(1, num_groups);
        
        % inliers part
        for jj = 1:num_groups
            B1_SUB_solved{jj} = B1_TOTAL(:, ggrps==jj);
            [U, ~, ~] = svd(X(:,ggrps==jj));
            Bases(:,:,jj) = U(:, 1:rr);
            A1_SUB_solved{jj} = cell(D, 1);
            A1_row_ind_SUB_solved{jj} = cell(D, 1);
            An_SUB_solved{jj} = cell(D, parts);
            An_row_ind_SUB_solved{jj} = cell(D, parts);
            for ii = 1:D
                A1_SUB_solved{jj}{ii} = A1_TOTAL{ii}(ggrps==jj);
%                 A1_row_ind_SUB_solved{jj}{ii} = A1_row_ind_TOTAL{ii}(ggrps==jj);
                A1_row_ind_SUB_solved{jj}{ii} = 1:sum(ggrps==jj);
            end
        end
        
        % outliers part
        num_each_part = length(An_TOTAL{1,1});
        N = parts*num_each_part;
        Y = zeros(D, N);
        for ii = 1:D
            a = [];
            for kk = 1:parts
                a = [a, An_TOTAL{ii, kk}];
            end
            Y(ii,:) = a;
        end
        outliers_class_solved = zeros(1, N);
        for t = 1:N
            y = Y(:, t);
            p = outlier_classification(y, Bases);
            outliers_class_solved(t) = p;
            part_label = ceil(t/num_each_part);
            part_ind = mod(t,num_each_part);
            if part_ind == 0
                part_ind = num_each_part;
            end
            for ii = 1:D
                An_SUB_solved{p}{ii, part_label}(end+1) = An_TOTAL{ii, part_label}(part_ind);
                An_row_ind_SUB_solved{p}{ii, part_label}(end+1)= An_row_ind_TOTAL{ii, part_label}(part_ind);
            end
        end
        for kk = 1:parts
            for jj = 1:num_groups
                for ii = 1:D
                    An_row_ind_SUB_solved{jj}{ii, kk} = 1:length(An_row_ind_SUB_solved{jj}{ii, kk});
                end
            end
        end
        
        % outlier_classification error
        gt_outliers_class = [];
        basic_count = floor(num_each_part / num_groups);
        remainder = mod(num_each_part, num_groups);
        cc = [];
        for jj = 1:num_groups
            count = basic_count + (jj <= remainder);
            cc = [cc, jj * ones(1, count)];
        end
        for kk = 1:parts
            gt_outliers_class = [gt_outliers_class, cc];
        end
        assert(length(gt_outliers_class) == N, 'gt_outliers_class length must be %d.', N);
        missrate_out = sum(gt_outliers_class ~= outliers_class_solved)/N;
        fprintf('\tmissrate_out (Outlier Classification): %.4f\n', missrate_out)
end

function outlier_class_solved = outlier_classification(y,U_solved)
    [~, ~, num_groups] = size(U_solved);
    L = zeros(1,num_groups);
    for jj = 1:num_groups
        [~, res] = LSR_v3(U_solved(:,:,jj), y, 0.5);
        L(jj) = res.subspace_cos_dist;
    end
    [~,outlier_class_solved] = min(L);
end