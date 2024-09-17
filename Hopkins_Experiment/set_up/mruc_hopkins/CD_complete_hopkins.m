function [best_test_err,best_P,best_B,final_cert,history,result]=CD_complete_hopkins(B1,Bn,A1,An,A1_row_ind,An_row_ind,init_rank, ... 
    eps_init,verbose,max_out_iter,max_in_iter,test_ind,test_label,tol,omega,beta,eps_min,perm_label,default_lr,gradient_step,rawimage)
    if nargin < 21
        print_err = false;
    else
        print_err = true;
    end
    [Bn,test_ind,test_label] = clean(Bn,test_ind,test_label);
%     Bn=process_Bn(Bn);
    final_cert=1;
    fix=false;
    n=size(B1,1);
    one_vector=ones(n,1); 
    parts=length(Bn);
    Bn_gradient = cell(parts,1);
    softP = cell(parts,1);
    softP_old = cell(parts,1);
    PP = cell(parts,1);
    lr = zeros(parts,1);
    B1_num_columns=size(B1,2);
    Bn_num_columns=zeros(parts,1);
    offset=B1_num_columns;
    Bn_ind=cell(parts,1);
    now_rank=init_rank;
    B_raw = B1;
    for pp=1:parts
        Bn_num_columns(pp)=size(Bn{pp},2);
        Bn_gradient{pp} = zeros(n,Bn_num_columns(pp));
        lr(pp)=0;
        softP{pp}=zeros(n,n);  
        Bn_ind{pp}=[offset+1,offset+Bn_num_columns(pp)];
        offset=offset+Bn_num_columns(pp);
        B_raw = [B_raw, Bn{pp}];
    end
    if print_err
        err = Evaluate(B_raw,rawimage);
        fprintf('Evaluate(B_raw,rawimage) = %.4f\n--\n',err)
    end
    best_test_err=inf(parts+1,1);
    history=zeros(6,0);
    cert=zeros(parts,1);
    diff=zeros(parts,1);
    result=[99999999999,999999999999,9999999999]';
    counts=0;
    eps_counts=0;
    best_B1=B1;
    best_Bn = Bn;
    eps=eps_init;
    B1=best_B1;
    Bn=best_Bn;  
    lambda=999;
    for out=1:max_out_iter
        counts=counts+1;
        eps_counts=eps_counts+1;
        for ii=1:n
           B1(ii,A1_row_ind{ii})=A1{ii};
        end
        for gg=1:gradient_step
            B=B1;
            for pp=1:parts 
                Bn{pp}=Bn{pp}-lr(pp)*Bn_gradient{pp};
                B=[B,Bn{pp}];
            end
            if sum(lr)>=0.02
                [B,last_sig]=ProxNC(B,now_rank);
                lambda=last_sig(now_rank);
            end
            B1=B(:,1:B1_num_columns);
            for pp=1:parts
                indx=Bn_ind{pp};
                nowB=B(:,indx(1):indx(2));
                Bn{pp}=nowB;
            end
        end
        rvs=0;
        for pp=1:parts
            indx=Bn_ind{pp};
            nowB=B(:,indx(1):indx(2));
%                 Bn{pp}=nowB;
            C=zeros(n);
            for ii=1:n
                right_diff_sum=zeros(1,n);
                for jj=1:n
                    right_diff_sum(jj)=sum((An{jj,pp}-nowB(ii,An_row_ind{jj,pp})).^2);
                end
                C(ii,:)=right_diff_sum;
            end
            % inner loop
            if ~fix
                f=zeros(n,1);
                g=zeros(n,1);
                softP_old{pp}=softP{pp};
                softPP_old=softP{pp};
                for kk=1:max_in_iter
                    f=Softmin(C-f*one_vector'-one_vector*g',eps)+f;
                    g=Softmin(C-f*one_vector'-one_vector*g',eps,2)+g;
                    if mod(kk,10)==0
                        softPP=exp((f*one_vector'+one_vector*g'-C)/eps);
                        if (norm(softPP-softPP_old,'fro')<1e-4) && (norm(sum(softPP,2)-1)<0.01)
                            break
                        end
                        softPP_old=softPP;
                    end
                end
                softP{pp}=softPP;
                [~, argmax] = max(softP{pp},[],2);
                PP{pp} = argmax;
            end
            Bgradient = zeros(n,Bn_num_columns(pp));
            for ii=1:n
                for jj=1:n            
                    Bgradient(ii,An_row_ind{jj,pp})=Bgradient(ii,An_row_ind{jj,pp})+softPP(ii,jj)*(nowB(ii,An_row_ind{jj,pp})-An{jj,pp});
                end
            end
            Bn_gradient{pp}=Bgradient;
            rvs=rvs+sum(sum(C.*softPP));
            cert(pp)=mean(abs(max(softPP,[],2)-1));
            diff(pp)=mean(abs(softPP-softP_old{pp}),'all');
            if default_lr~=0
                lr(pp)=default_lr;
            else
                lr(pp) = min(max((0.5*(1-cert(pp))^omega+0.01)*(1-min(beta*diff(pp),0.99)),0.01),0.15);
            end
        end
        nnorm=sum(svd(B));
        rvs_total=rvs+lambda*nnorm;
        test_err=zeros(parts+1,1);
        test_length=zeros(parts+1,1);
        test_err(1)=norm(B1(test_ind{1})-test_label{1}')^2; 
        test_length(1)=length(test_ind{1});
        perm_err=0;
        for testii=2:parts+1
            nowB=Bn{testii-1};
            test_err(testii)=norm(nowB(test_ind{testii})-test_label{testii}').^2;
            test_length(testii)=length(test_ind{testii});
            perm_err=perm_err+sum(abs(PP{testii-1}-perm_label{testii-1}'));
%             perm_err=perm_err+norm(PP{testii-1}-perm_label{testii-1},'fro')^2/2;
        end
        total_test_err=sqrt(sum(test_err)/sum(test_length));
        history(:,end+1)=[rvs_total;total_test_err;mean(lr);lambda;eps;perm_err];

        part_test_err=sqrt(test_err./test_length);
%         fprintf('out:%d, part_test_err_sum:%.4f\n',out,sum(part_test_err))
        
        if  (total_test_err<=(0.9999*result(3))) || any(part_test_err<=(0.9999*best_test_err))
            if (total_test_err<=(0.9999*result(3)))
                result(1)=rvs;
                result(2)=rvs_total;
                result(3)=total_test_err;
            end
            if part_test_err(1)<=(0.9999*best_test_err(1))
                best_B1=B1;
            end
            for ppp=1:parts
                if part_test_err(ppp+1)<=(0.9999*best_test_err(ppp+1))
                    best_Bn{ppp}=Bn{ppp};
                end
            end
            best_B=best_B1;
            for ppp=1:parts 
                best_B=[best_B,best_Bn{ppp}];
            end
            best_test_err(part_test_err<=(0.9999*best_test_err))=part_test_err(part_test_err<=(0.9999*best_test_err));
            best_P=PP;
            best_softP=softP;
            final_cert=mean(cert);
            counts=0;
            eps_counts=0;
        end

        if verbose && mod(out,50)==0
            fprintf('%dth iteration: Outer loop:(%.5f,%.5f,%.5f,%.5f,%.1f),\n\tcert:%.5f,diff:%.5f,lr:%.5f,eps:%.5f,lambda:%.5f,rank:%d\n',...
                      out,rvs,total_test_err,nnorm,rvs_total,perm_err,mean(cert),mean(diff),mean(lr),eps,lambda,rank(B));
            if print_err
                if mod(out,300) == 0
                    err = Evaluate(best_B,rawimage);
                    fprintf(' Evaluate(best_B,rawimage); = %.4f\n---\n\n',err)
                end
            end
        end
        if eps<eps_min
            fprintf('Earlier Termination~\n---\n')
            break
        end
        if eps_counts>=tol
%             break
            eps_counts=0;
            eps=eps/2;
            default_lr=default_lr*0.99;
        end
    end
end


function [mat,sigv]=ProxNC(B,rank)
    [U,S,V] = svd(B);
    sigv=diag(S);
    n=rank;
    mat=U(:,1:n)*S(1:n,1:n)*V(:,1:n)';
end

function value=Softmin(M,eps,axis)
    if nargin<=2
        axis=1;
    end
    if axis==1
        rmin=min(M,[],2);
        value=rmin-eps*log(sum(exp((rmin-M)/eps),2));
    else
        rmin=min(M);
        value=rmin-eps*log(sum(exp((rmin-M)/eps)));
        value=value';
    end
end

function perm = permutation_from_matrix(P)
    n = size(P, 1);
    perm = zeros(1, n);
    for i = 1:n
        perm(i) = find(P(i, :) == 1);
    end
end

function [clean_Bn,clean_ind,clean_label] = clean(Bn,test_ind,test_label)
    parts = length(Bn);
    clean_Bn = {};
    clean_ind = cell(parts+1,1);
    clean_label = cell(parts+1,1);
    clean_ind{1} = test_ind{1};
    clean_label{1} = test_label{1};
    idx = 1;
    c = 0;
    for ii = 1:parts
        if ~isempty(Bn{ii})
            clean_Bn{idx} = Bn{ii};
            clean_ind{idx+1} = test_ind{ii+1};
            clean_label{idx+1} = test_label{ii+1};
            idx = idx+1;
        else
            c = c+1;
        end
    end
    if c>0
        fprintf('\tclean works!\n')
    end
end