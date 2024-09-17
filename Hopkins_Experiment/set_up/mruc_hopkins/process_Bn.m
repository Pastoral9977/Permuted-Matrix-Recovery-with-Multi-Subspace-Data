function new_Bn = process_Bn(Bn)
parts = length(Bn);
new_Bn = Bn;
for ii = 1:parts
    nowB=Bn{ii};
    norm_constant = mean(abs(nowB(:)));
    A = randn(size(nowB));
    A = A./vecnorm(A);
    new_Bn{ii} = A*norm_constant;
end 
end