function [best_recover,softP,now_direction]=MRUC_FaceRecovery(max_in_iter,eps_init,max_out_iter,verbose,use_rank,init_image,sub_observes,lr1,lr2,np,mp,perm_index)
    recover_image=init_image;
    eps=eps_init;
    eps_counts=0;
    fix_eps=false;
    best_obj=9999;
    best_cert=1;
    best_recover=recover_image;
    length_perm=length(perm_index);
    sub_recovers={};
    lambda=1.2;
    [n,m] = size(init_image);
    now_direction=zeros(n,m); % same size with large elementwise matrix
    factor_ind_nn = n/np;
    factor_ind_mm = m/mp;
    for iter=1:max_out_iter
        eps_counts=eps_counts+1;
        recover_image_old=recover_image;

        for ii=1:length_perm
           ind1=mod(perm_index(ii)-1,np)+1;
           ind2=ceil(perm_index(ii)/np);
           sub_recovers{perm_index(ii)}=recover_image(((ind1-1)*factor_ind_nn+1):(ind1*factor_ind_nn),((ind2-1)*factor_ind_mm+1):(ind2*factor_ind_mm));
        end

        now_norm=sum(svd(recover_image));
        cost=compute_cost(sub_recovers,sub_observes,perm_index);  % size: (length_perm,length_perm)
        softP=Sinkhorn(cost,eps,max_in_iter);  % size: (length_perm,length_perm)
        now_obj=sum(cost.*softP,'all')+lambda*now_norm;

        for ii=1:(np*mp)
           ind1=mod(ii-1,np)+1;
           ind2=ceil(ii/np);
           now_recover=recover_image(((ind1-1)*factor_ind_nn+1):(ind1*factor_ind_nn),((ind2-1)*factor_ind_mm+1):(ind2*factor_ind_mm));
           if ismember(ii,perm_index)
               direction=zeros(factor_ind_nn, factor_ind_mm); % same size with each block, not large blockwise matrix
               for jj=1:length_perm
                   direction=direction+softP(find(perm_index==ii),jj)*sub_observes{perm_index(jj)};   % may match with Bgradient in CD_complete.m
               end
               obs_ind=(direction~=0);
               now_recover(obs_ind)=(1-lr1)*now_recover(obs_ind)+lr1*direction(obs_ind);
               now_direction(((ind1-1)*factor_ind_nn+1):(ind1*factor_ind_nn),((ind2-1)*factor_ind_mm+1):(ind2*factor_ind_mm))=direction;
           else
               now_observes=sub_observes{ii};
               now_obs_ind=(now_observes~=0);
               %now_recover(now_obs_ind)=now_observes(now_obs_ind);
               now_recover(now_obs_ind)=(1-lr2)*now_recover(now_obs_ind)+lr2*now_observes(now_obs_ind);
               now_direction(((ind1-1)*factor_ind_nn+1):(ind1*factor_ind_nn),((ind2-1)*factor_ind_mm+1):(ind2*factor_ind_mm))=now_observes;
           end 
           recover_image(((ind1-1)*factor_ind_nn+1):(ind1*factor_ind_nn),((ind2-1)*factor_ind_mm+1):(ind2*factor_ind_mm))=now_recover;
        end
        [lambda,recover_image]=best_approx(recover_image,use_rank);

        if now_obj<(0.9999*best_obj)
            best_obj=now_obj;
            eps_counts=0;
            best_recover=recover_image;
        end
        
        cert=mean(abs(max(softP,[],2)-1));
        if cert<0.001
            [~,recover_perm]=max(softP,[],2);
            while norm(recover_image_old-recover_image,'fro')>0.001
                recover_image_old=recover_image;
                [~,recover_image]=best_approx(recover_image,15);
                for ii=1:(np*mp)
                   ind1=mod(ii-1,np)+1;
                   ind2=ceil(ii/np);
                   now_recover=recover_image(((ind1-1)*factor_ind_nn+1):(ind1*factor_ind_nn),((ind2-1)*factor_ind_mm+1):(ind2*factor_ind_mm));
                   if ismember(ii,perm_index)
                       now_observes=sub_observes{perm_index(recover_perm(find(perm_index==ii)))};
                   else
                       now_observes=sub_observes{ii};

                   end 
                   now_obs_ind=(now_observes~=0);
                   now_recover(now_obs_ind)=now_observes(now_obs_ind);
                   now_direction(((ind1-1)*factor_ind_nn+1):(ind1*factor_ind_nn),((ind2-1)*factor_ind_mm+1):(ind2*factor_ind_mm))=now_observes;
                   recover_image(((ind1-1)*factor_ind_nn+1):(ind1*factor_ind_nn),((ind2-1)*factor_ind_mm+1):(ind2*factor_ind_mm))=now_recover;
                end
            end
            best_recover=recover_image;
            break
        end

        if cert<(0.9999*best_cert)
            best_cert=cert;
            eps_counts=0;
            best_recover=recover_image;
        end

        if (eps_counts>100) && (~fix_eps) && (eps>0.0001)
            eps=eps/2;
            eps_counts=0;
            recover_image=best_recover;
        end

        
        if verbose
            fprintf('%dth iteration,norm:%.5f,obj:%.5f,cert:%.5f,eps:%.5f \n',iter,now_norm,now_obj,eps);
        end
    end
end

function C=compute_cost(sub_recovers,sub_observes,perm_index)
    perm_n=length(perm_index);
    C=zeros(perm_n);
    for ii=1:perm_n
        now_recover=sub_recovers{perm_index(ii)};
        for jj=1:perm_n
            now_obs=sub_observes{perm_index(jj)};
            obs_ind=(now_obs>0);
            C(ii,jj)=norm(now_recover(obs_ind)-now_obs(obs_ind),'fro')^2;
        end
    end 
end

function softP=Sinkhorn(C,eps,max_iter)
    [n,m]=size(C);
    f=zeros(n,1);
    g=zeros(m,1);
    one_vector_f=ones(m,1);
    one_vector_g=ones(n,1);
    for kk=1:max_iter
        f=Softmin(C-f*one_vector_f'-one_vector_g*g',eps)+f;
        g=Softmin(C-f*one_vector_f'-one_vector_g*g',eps,2)+g;
        if mod(kk,10)==0
            softP=exp((f*one_vector_f'+one_vector_g*g'-C)/eps);
            if norm(sum(softP,2)-1)<0.001
                break
            end
        end
    end
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

