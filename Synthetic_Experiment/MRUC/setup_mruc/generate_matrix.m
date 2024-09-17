function [B1,A1,An,A1_row_ind,An_row_ind,oracle_image,permutation_matrix,test_ind,test_label,norm_constant,rawimage,outlier_perm_belongings]...
    =generate_matrix(D,ranks,parts,missing_percent,shuffle_columns,perm_num,width,snr,seed_for_rawimage_generation,seed_for_permutation_matrix)

rng(seed_for_rawimage_generation)

A=randn(D,ranks);
C=randn(ranks,width);
noimage=A*C;    
rawimage=noimage + randn(D,width)*snr;  
[n,m]=size(rawimage); 
shuffle_columns_per_part=shuffle_columns/parts;
rawimage=rawimage*2/max([n,m]);        
noimage=noimage*2/max([n,m]);          
norm_constant=mean(abs(rawimage(:)));  
permutation_label=cell(parts,1);       
permutation_matrix=cell(parts,1);
R1=rawimage(:,1:(end-shuffle_columns)); 
A1_row_ind=cell(n,1);                   
ind=1:n*(m-shuffle_columns);    
selected_missing_entries = randperm(n*(m-shuffle_columns), round(n*(m-shuffle_columns)*missing_percent/100));
R1(ind(selected_missing_entries))=0; 
oracle_image=R1;
A1=cell(n,1);
test_ind=cell(parts,1);
test_label=cell(parts,1);
for ii=(round(n*(m-shuffle_columns)*missing_percent/100)+1):(n*(m-shuffle_columns))
    row_ind=mod(ind(ii)-1,n)+1;          
    col_ind=fix((ind(ii)-1)/n)+1;        
    A1_row_ind{row_ind}(end+1)=col_ind;  
end
for ii=1:n
    A1{ii}=R1(ii,A1_row_ind{ii}); 
end
An=cell(n,parts);                 
An_row_ind=cell(n,parts);
test_ind{1}=[1:(m-shuffle_columns)*n]';    
nowdata=noimage(:,1:(end-shuffle_columns)); 
test_label{1}=nowdata(:);  

rng(seed_for_permutation_matrix)
for ii=1:parts
    perm_matrix=zeros(n,n); 
    perm_label=[randperm(perm_num),(perm_num+1):n]; 
    pmat = rawimage(perm_label,(end-shuffle_columns+(ii-1)*shuffle_columns_per_part+1):(end-shuffle_columns+ii*shuffle_columns_per_part));  
    nowdata=noimage(:,(end-shuffle_columns+(ii-1)*shuffle_columns_per_part+1):(end-shuffle_columns+ii*shuffle_columns_per_part)); 
    test_ind{ii+1}=[1:(shuffle_columns_per_part)*n]';  
    test_label{ii+1}=nowdata(:);
    ind=1:n*shuffle_columns_per_part; 
    selected_missing_entries = randperm(n*shuffle_columns_per_part, round(n*shuffle_columns_per_part*missing_percent/100));
    pmat(ind(selected_missing_entries))=0;  
    for jj=(round(n*shuffle_columns_per_part*missing_percent/100)+1):(n*shuffle_columns_per_part) 
        row_ind=mod(ind(jj)-1,n)+1;
        col_ind=fix((ind(jj)-1)/n)+1;
        An_row_ind{row_ind,ii}(end+1)=col_ind;
    end 
    for ss=1:n
        perm_matrix(ss,perm_label(ss))=1;
    end
    perm_matrix=perm_matrix';   
    for jj=1:n
        An{jj,ii}=pmat(jj,An_row_ind{jj,ii});
    end
    permutation_label{ii}=perm_label;
    permutation_matrix{ii}=perm_matrix;
    oracle_image=[oracle_image,perm_matrix*pmat];
end
B1 = R1;
outlier_perm_belongings = [];
for kk = 1:parts
    outlier_perm_belongings = [outlier_perm_belongings, kk*ones(1, shuffle_columns_per_part)];
end
end