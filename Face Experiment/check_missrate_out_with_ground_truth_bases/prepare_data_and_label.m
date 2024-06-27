%% Prepare data and label
function prepare_data_and_label()
    global LABEL INIT
    run Data_Loading.m
    run Initialization.m
    run Data_Construction.m

    allpoints_class_gt = zeros(N, 1);
    for i = 1:INIT.num_groups
        allpoints_class_gt(sum(nn(1:i))+1:sum(nn(1:i+1))) = i;
    end
    LABEL.allpoints_class_gt = allpoints_class_gt;
    
end