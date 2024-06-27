%% Assessment
function Assessment()
    fprintf("Assessment...\n")

    global Involve_reconstructed_bases Involve_gt_bases DATA LABEL
    
    % For ground truth bases
    if Involve_gt_bases
        LABEL.allpoints_class_solved_gtbases = zeros(DATA.N, 1);
        for i = 1:DATA.N
            if (ismember(i, LABEL.outliers_ID) == 1)
                LABEL.allpoints_class_solved_gtbases(i) = LABEL.outlier_class_solved_gtbases(LABEL.outliers_ID == i);
            end
        end
        
        % Calculate missrate for gt bases
        missrate_out_gt = sum(LABEL.allpoints_class_solved_gtbases(LABEL.outliers_ID) ~= LABEL.allpoints_class_gt(LABEL.outliers_ID)) / length(LABEL.outliers_ID);
        fprintf('\tmissrate_out with ground truth bases = %0.4f\n', missrate_out_gt);
    end

    % For reconstructed bases
    LABEL.allpoints_class_solved_reconstructedbases = zeros(DATA.N, 1);
    if Involve_reconstructed_bases
        for i = 1:DATA.N
            if ismember(i, LABEL.Outliers_id)
                LABEL.allpoints_class_solved_reconstructedbases(i) = LABEL.outlier_class_solved_reconstructedbases(LABEL.Outliers_id == i);
            end
        end
        
        % Calculate missrate for reconstructed bases
        undetected_outlier_num = length(setdiff(LABEL.outliers_ID, LABEL.Outliers_id));
        detected_true_Outliers_id = setdiff(LABEL.Outliers_id, LABEL.inliers_ID);
        missclassified_outlier_num = sum(LABEL.allpoints_class_solved_reconstructedbases(detected_true_Outliers_id) ~= LABEL.allpoints_class_gt(detected_true_Outliers_id));
        missrate_out_reconstructed = (missclassified_outlier_num+undetected_outlier_num) / length(LABEL.outliers_ID);
        undetected_outlier_num_ratio = undetected_outlier_num/length(LABEL.outliers_ID);
        fprintf('\tmissrate_out with reconstructed bases = %0.4f, where undetected_outlier_num_ratio = %.4f\n\n', missrate_out_reconstructed, undetected_outlier_num_ratio);
    end
    
end