function [init_image, sub_matrices] = initialize_image(perm_miss_image, np, mp)
    % Assume raw image to be of size 48*42
    
    [n,m] = size(perm_miss_image);
    sub_matrices={np*mp};
    sub_known_ind={np*mp};
    factor_nn=n/np;
    factor_mm=m/mp;
    init_image=zeros(n, m);
    lambda=1.5;

    for ii=1:np
        for jj=1:mp
            sub_matrices{ii+(jj-1)*np}=perm_miss_image(((ii-1)*factor_nn+1):(ii*factor_nn),((jj-1)*factor_mm+1):(jj*factor_mm));
            sub_known_ind{ii+(jj-1)*np}=(sub_matrices{ii+(jj-1)*np}~=0);
            nowB=sub_matrices{ii+(jj-1)*np}; 
            B=nowB;
            B_old=zeros(n/np,m/mp);
            while norm(B_old-B,'fro')>0.0001  
                B_old=B;
                B(sub_known_ind{ii+(jj-1)*np})=nowB(sub_known_ind{ii+(jj-1)*np});
                B=ProxNC(B,lambda);  
            end
            B(sub_known_ind{ii+(jj-1)*np})=nowB(sub_known_ind{ii+(jj-1)*np});
            init_image(((ii-1)*factor_nn+1):(ii*factor_nn),((jj-1)*factor_mm+1):(jj*factor_mm))=B;
        end
    end
    [~,init_image]=best_approx(init_image,1);

end

function mat=ProxNC(B,lambda)
    [U,S,V] = svd(B);
    sigv=diag(S);
    n=length(sigv(sigv>lambda));
    mat=U(:,1:n)*S(1:n,1:n)*V(:,1:n)';
end