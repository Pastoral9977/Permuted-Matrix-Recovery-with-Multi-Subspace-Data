function U_gt = ground_truth_bases_construction(data, allpoints_class_gt, rank)
    
    if length(allpoints_class_gt) ~= size(data, 2)
        error("The input data can not match the ground truth label!")
    end
    
    classes = unique(allpoints_class_gt);
    num_groups = length(classes);
    [V, ~] = size(data);
    U_gt = ones(V, rank, num_groups);

    for idx = 1:num_groups
        class_label = classes(idx);
        Mat = data(:, allpoints_class_gt == class_label);
        [U, ~, ~] = svd(Mat);
        U_gt(:, :, idx) = U(:, 1:rank);
    end

end