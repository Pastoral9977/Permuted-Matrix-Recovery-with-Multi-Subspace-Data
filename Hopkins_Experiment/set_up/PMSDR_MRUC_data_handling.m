function [DATA_SUB_solved, seq]=PMSDR_MRUC_data_handling(DATA_TOTAL,num_groups,rr,s_in,s_out)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seq map indices in [B1_TOTAL,Bn_TOTAL{1},...,Bn_TOTAL{parts}] (which is identical to X_new_tilde = X_tilde(:,[id_in,id_out]))
% into 
% [B1_SUB_solved{1},Bn_SUB_solved{1}{1},...,Bn_SUB_solved{1}{parts},...,Bn_SUB_solved{num_groups}{1},...,Bn_SUB_solved{num_groups}{parts}]
% (which is identical to X_new_tilde(:,seq))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5
    print_missrate_out=false;
else
    print_missrate_out=true;
end

B1_TOTAL=DATA_TOTAL.B1_TOTAL;
Bn_TOTAL=DATA_TOTAL.Bn_TOTAL;
A1_TOTAL=DATA_TOTAL.A1_TOTAL;
An_TOTAL=DATA_TOTAL.An_TOTAL;
An_row_ind_TOTAL=DATA_TOTAL.An_row_ind_TOTAL;
test_ind_TOTAL=DATA_TOTAL.test_ind_TOTAL;
test_label_TOTAL=DATA_TOTAL.test_label_TOTAL;

X_in=B1_TOTAL;
[D,parts]=size(An_TOTAL);

X_gt = X_in;
for ii = 2:parts+1
    X_gt = [X_gt,reshape(test_label_TOTAL{ii},[D,length(test_label_TOTAL{ii})/D])];
end

test_inds_solved=cell(num_groups,1);
test_labels_solved=cell(num_groups,1);
for jj = 1:num_groups
    test_labels_solved{jj}=cell(parts+1,1);
    test_inds_solved{jj}=cell(parts+1,1);
end


%% subspace clustering
fprintf('--subspace clustering--\n')
r=0; 
affine=true; 
outlier=true; 
alpha=500;
rho=0.7; 
broadcast=false;
[missrate_in,ggrps]=SSC(X_in,r,affine,alpha,outlier,rho,s_in);
fprintf('\tmissrate_in(Subspace Clustering):%.4f\n', missrate_in)

Bases=zeros(D,rr,num_groups);
B1_SUB_solved=cell(num_groups,1);
Bn_SUB_solved=cell(num_groups,1);
A1_SUB_solved=cell(num_groups,1);
A1_row_ind_SUB_solved=cell(num_groups,1);
An_SUB_solved=cell(num_groups,1);
An_row_ind_SUB_solved=cell(num_groups,1);

%% inlierspart
seq_in=cell(num_groups,1);
for jj=1:num_groups
    B1_SUB_solved{jj}=B1_TOTAL(:,ggrps==jj);
    seq_in{jj}=find(ggrps==jj);
    Bn_SUB_solved{jj}=cell(parts,1);
    X = X_in(:,ggrps==jj);
    test_labels_solved{jj}{1} = X(:);
    test_inds_solved{jj}{1} = 1:length(X(:));
    [U,~,~]=svd(X);
    Bases(:,:,jj)=U(:,1:rr);
    A1_SUB_solved{jj}=cell(D,1);
    A1_row_ind_SUB_solved{jj}=cell(D,1);
    An_SUB_solved{jj}=cell(D,parts);
    An_row_ind_SUB_solved{jj}=cell(D,parts);
    for ii=1:D
        A1_SUB_solved{jj}{ii}=A1_TOTAL{ii}(ggrps==jj);
        A1_row_ind_SUB_solved{jj}{ii}=1:sum(ggrps==jj);
    end
end

%% outlierspart
part_nums=[];
for jj=1:parts
   part_nums=[part_nums, length(An_TOTAL{1,jj})];
end

N_out=sum(part_nums);
Y=[];
for jj=1:parts
    Y=[Y, Bn_TOTAL{jj}];
end

seq_out=cell(num_groups,parts);
N_in=length(s_in);
outlier_classes=[];
X_out_in_parts = cell(num_groups,parts);
for t=1:N_out
    y=Y(:,t);
    subspace_label=outlier_classification(y,Bases);
    outlier_classes=[outlier_classes,subspace_label];
    [part_label,part_ind]=part_handle(t,part_nums);
%     fprintf('subspace:%d, part:%d\n',subspace_label,part_label)
    X_out_in_parts{subspace_label,part_label} = [X_out_in_parts{subspace_label,part_label},X_gt(:,N_in+t)];
    Bn_SUB_solved{subspace_label}{part_label}=[Bn_SUB_solved{subspace_label}{part_label},y];
    seq_out{subspace_label,part_label}=[seq_out{subspace_label,part_label},t+N_in];
    for ii=1:D
        An_SUB_solved{subspace_label}{ii,part_label}(end+1)=An_TOTAL{ii,part_label}(part_ind);
        An_row_ind_SUB_solved{subspace_label}{ii,part_label}(end+1)=An_row_ind_TOTAL{ii,part_label}(part_ind);
    end
end

if print_missrate_out
    if N_out~=length(s_out)
        error('Something Wrong with N_out')
    end
    missrate_out=sum(s_out~=outlier_classes')/N_out;
    fprintf('\tmissrate_out(Outlier Classification):%.4f\n',missrate_out)
end

for kk=1:parts
    for jj=1:num_groups
        for ii=1:D
            An_row_ind_SUB_solved{jj}{ii, kk}=1:length(An_row_ind_SUB_solved{jj}{ii, kk});
        end
    end
end

%% test part
for jj = 1:num_groups
    for kk = 2:parts+1
        test_labels_solved{jj}{kk}=X_out_in_parts{jj,kk-1}(:);
        test_inds_solved{jj}{kk}=1:length(test_labels_solved{jj}{kk}(:));
    end
end
 
%% sequence part
seq=[];
for jj=1:num_groups
    seq=[seq;seq_in{jj}];
    for kk=1:parts
        seq=[seq;seq_out{jj,kk}'];
    end
end

%% SUM UP
DATA_SUB_solved=struct();
DATA_SUB_solved.B1_SUB_solved=B1_SUB_solved;
DATA_SUB_solved.Bn_SUB_solved=Bn_SUB_solved;
DATA_SUB_solved.A1_SUB_solved=A1_SUB_solved;
DATA_SUB_solved.A1_row_ind_SUB_solved=A1_row_ind_SUB_solved;
DATA_SUB_solved.An_SUB_solved=An_SUB_solved;
DATA_SUB_solved.An_row_ind_SUB_solved=An_row_ind_SUB_solved;
DATA_SUB_solved.test_inds_solved=test_inds_solved;
DATA_SUB_solved.test_labels_solved=test_labels_solved;

end

function outlier_class_solved=outlier_classification(y,U_solved)
    [~, ~, num_groups]=size(U_solved);
    L=zeros(1,num_groups);
    for jj=1:num_groups
        [~, res]=LSR_v3(U_solved(:,:,jj), y, 0.5);
        L(jj)=res.subspace_cos_dist;
    end
    [~,outlier_class_solved]=min(L);
end

function [part_label,part_ind]=part_handle(t,part_nums)
part_label=0;
now=0;
while now < t
   part_label=part_label+1;
    now=now+part_nums(part_label);
end
part_ind=t-now+part_nums(part_label);
end