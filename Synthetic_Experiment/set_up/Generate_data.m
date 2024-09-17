function [X_gt, U_gt, Subspace_Label_gt] = Generate_data(D,V,num_groups,rrank)

    U_gt = zeros(D,rrank,num_groups);
    XX = [];
    for ii = 1:num_groups
        [U_gt(:,:,ii),~] = qr(randn(D,rrank), 0);
        B_initial = U_gt(:,:,ii) * randn(rrank,V);
        XX = [XX B_initial./vecnorm(B_initial)];
    end
    
    X_gt = XX;
    Subspace_Label_gt = zeros(num_groups*V,1);
    for ii = 1:num_groups
        for jj = V*(ii-1)+1:V*ii
            Subspace_Label_gt(jj) = ii;
        end
    end
    
end
    
    
    
    
        
        
        
        
        
        