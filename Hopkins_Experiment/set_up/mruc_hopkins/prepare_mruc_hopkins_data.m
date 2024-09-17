function [DATA_SUB,DATA_TOTAL,perm_labels,X_tilde,outliers_ID,seq_gt] = prepare_mruc_hopkins_data(X_gt,parts,s,shuffled_ratio,outlier_ratio,seed)

rng(seed)
n = size(X_gt,1);
[X_tilde,outliers_ID_in_parts,perm_labels] = permute_data(X_gt, shuffled_ratio, outlier_ratio, s, parts, seed);
outliers_ID = [];
for kk = 1:parts
    outliers_ID = [outliers_ID, outliers_ID_in_parts{kk}];
end

inliers_ID = setdiff(1:size(X_gt,2),outliers_ID);
num_groups = length(unique(s));
A1s=cell(num_groups,1); 
A1_row_inds=cell(num_groups,1); 
Ans=cell(num_groups,1);                 
An_row_inds=cell(num_groups,1);
B1s=cell(num_groups,1);
Bns=cell(num_groups,1);
test_inds=cell(num_groups,1);
test_labels=cell(num_groups,1);
seq_in=cell(num_groups,1);
seq_out=cell(num_groups,parts);
%% SUB
for jj = 1:num_groups
    %%
    test_inds{jj} = cell(parts+1,1);
    test_labels{jj} = cell(parts+1,1);
    in_inds_jj = inliers_or_outliers_in_subspaceJ(jj,s,inliers_ID);
    seq_in{jj}=in_inds_jj;
    XX = X_gt(:,in_inds_jj);
    B1s{jj} = XX;
    test_inds{jj}{1} = 1:length(XX(:));
    test_labels{jj}{1} = XX(:);
    
    %%
    A1s{jj}=cell(n,1);
    A1_row_inds{jj}=cell(n,1);
    for ii=1:n
        A1s{jj}{ii}=X_gt(ii,in_inds_jj);
        A1_row_inds{jj}{ii}=1:length(in_inds_jj);
    end
    
    %%
    Ans{jj}=cell(n,parts);                 
    An_row_inds{jj}=cell(n,parts);
    Bns{jj}=cell(parts,1);
    for kk = 1:parts
        out_inds_jj_kk = inliers_or_outliers_in_subspaceJ(jj, s, outliers_ID_in_parts{kk});
        seq_out{jj,kk}=out_inds_jj_kk;
        pmat = X_tilde(:,out_inds_jj_kk);
        Bns{jj}{kk} = pmat;
%         Bns{jj}{kk} = randn(n,length(out_inds_jj_kk))*mean(abs(pmat(:)));
        for ii = 1:n
            Ans{jj}{ii,kk}=pmat(ii,:);
            An_row_inds{jj}{ii,kk}=1:length(out_inds_jj_kk);
        end
        XX = X_gt(:,out_inds_jj_kk);
        test_inds{jj}{kk+1} = 1:length(XX(:));
        test_labels{jj}{kk+1} = XX(:);
    end
end

%% seq_gt
seq_gt=[];
for jj=1:num_groups
    seq_gt=[seq_gt;seq_in{jj}];
    for kk=1:parts
        seq_gt=[seq_gt;seq_out{jj,kk}];
    end
end

%% TOTAL
B1_TOTAL = X_tilde(:, inliers_ID);

test_ind_TOTAL=cell(parts+1,1);
test_label_TOTAL=cell(parts+1,1);
test_ind_TOTAL{1} = 1:n*length(inliers_ID);
test_label_TOTAL{1} = B1_TOTAL(:);

A1_TOTAL = cell(n,1);
A1_row_ind_TOTAL = cell(n,1);
for ii = 1:n
    A1_TOTAL{ii}=B1_TOTAL(ii,:);
    A1_row_ind_TOTAL{ii}=1:length(inliers_ID);
end

Bn_TOTAL = cell(parts,1);
An_TOTAL = cell(n,parts);
An_row_ind_TOTAL = cell(n,parts);
for kk = 1:parts
    id_out_kk = outliers_ID_in_parts{kk};
    pmat = X_tilde(:,id_out_kk);
    tr_mat = X_gt(:,id_out_kk);
    Bn_TOTAL{kk} = pmat;
    for ii = 1:n
        An_TOTAL{ii,kk} = X_tilde(ii,id_out_kk);
        An_row_ind_TOTAL{ii,kk} = 1:length(id_out_kk);
    end
    test_label_TOTAL{kk+1}=tr_mat(:);
    test_ind_TOTAL{kk+1}=1:length(tr_mat(:));
end
    
%% SUM UP
DATA_SUB = struct();
DATA_TOTAL = struct();

DATA_SUB.B1s=B1s;
DATA_SUB.Bns=Bns;
DATA_SUB.A1s=A1s;
DATA_SUB.Ans=Ans;
DATA_SUB.A1_row_inds=A1_row_inds;
DATA_SUB.An_row_inds=An_row_inds;
DATA_SUB.test_inds=test_inds;
DATA_SUB.test_labels=test_labels;

DATA_TOTAL.B1_TOTAL=B1_TOTAL;
DATA_TOTAL.Bn_TOTAL=Bn_TOTAL;
DATA_TOTAL.A1_TOTAL=A1_TOTAL;
DATA_TOTAL.An_TOTAL=An_TOTAL;
DATA_TOTAL.A1_row_ind_TOTAL=A1_row_ind_TOTAL;
DATA_TOTAL.An_row_ind_TOTAL=An_row_ind_TOTAL;
DATA_TOTAL.test_ind_TOTAL=test_ind_TOTAL;
DATA_TOTAL.test_label_TOTAL=test_label_TOTAL;
end

function inds = inliers_or_outliers_in_subspaceJ(j, s, inliers_or_outliers_ID)
    ids = find(s==j);                                               % find(): find all non-zero entries' indices, not logical
    inds = ids(ismember(ids,inliers_or_outliers_ID));               % all inliers in group j
end