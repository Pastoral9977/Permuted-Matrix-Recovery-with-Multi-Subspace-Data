function [X_gt, U_gt, Subspace_Label_gt] = Generate_data(D,V,num_groups,rrank)

    M_gt = zeros(D,V,num_groups);
    X = zeros(D,num_groups*V);
    U_gt = zeros(D,rrank,num_groups);

            
    for i = 1:num_groups
        [Basis,~,~] = svd(randn(D,rrank));
        U_gt(:,:,i) = Basis(:,1:rrank);
        B_initial = U_gt(:,:,i) * randn(rrank,V);
        M_gt(:,:,i) = B_initial./vecnorm(B_initial);
        X(:,(i-1)*V+1:i*V) = M_gt(:,:,i);
    end

    Subspace_Label_gt = zeros(size(X,2),1);
    for i = 1:num_groups
        for j = V*(i-1)+1:V*i
            Subspace_Label_gt(j) = i;
        end
    end

    X_gt = reshape(M_gt, [D, V*num_groups]);
       % Ground_truth data and labels will be involved in assessment of data recovery
       % procedure.

end
    
    
    
    
        
        
        
        
        
        