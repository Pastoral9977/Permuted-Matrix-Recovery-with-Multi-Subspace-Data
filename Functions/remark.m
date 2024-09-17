function [test] = remark(Inliers_id, Outliers_id, inliers_ID, outliers_ID)
    TP = length(intersect(Inliers_id, inliers_ID)); 
    FP = length(intersect(Inliers_id, outliers_ID)); 
    FN = length(intersect(Outliers_id, inliers_ID)); 
    TN = length(intersect(Outliers_id, outliers_ID)); 
    
    test = struct;
    
    test.accuracy = (TP+TN)/(TP+TN+FP+FN);
    test.precision = TP/(TP+FP);
    test.recall = TP/(TP+FN);
    test.specificity = TN/(TN+FP);
    test.F1 = 2 * test.precision * test.recall / (test.precision + test.recall);
    
    fprintf('\t\t number of detected inliers/ true inliers = %d/%d ;\n',   length(Inliers_id),  length(inliers_ID) )
    fprintf('\t\t number of detected inliers/ true outliers = %d/%d ;\n',   length(Outliers_id),  length(outliers_ID) )
    fprintf('\t\t number of undetected true inliers/ true inliers = %d/%d;\n', length(setdiff(inliers_ID,Inliers_id)), length(inliers_ID))
    fprintf('\t\t number of undetected true outliers/ true outliers = %d/%d;\n', length(setdiff(outliers_ID,Outliers_id)), length(outliers_ID))

end