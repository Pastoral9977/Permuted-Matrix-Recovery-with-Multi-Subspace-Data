function Bn=generate_Bn(parts,n,shuffle_columns_per_parts,norm_constant,seed)

rng(seed)

Bn = cell(parts,1);
for pp=1:parts
    Bn{pp} = randn(n,shuffle_columns_per_parts(pp))*norm_constant;
end

end