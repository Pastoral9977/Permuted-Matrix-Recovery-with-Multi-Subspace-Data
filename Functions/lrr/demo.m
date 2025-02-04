function [] = demo()
%This routine demonstrates an example of using LRR to do subspace segmentation. We cosntruct 5 independent subspaces, each of which has a rank of 10, 
%sample 200 points of dimension 100 from each subspae, and randomly choose some points to corrupt.
[X,cids] = generate_data();
ls = [0.0001 0.0005 0.001 0.002 0.004 0.008 0.01 0.02 0.04 0.08 0.1]; %parameter lambda
rs = [];
accs=[];
for i=1:length(ls)
    Z = solve_lrr(X,X,ls(i));
    L = abs(Z)+abs(Z');
    disp('Perfoming NCut ...');
    idx = clu_ncut(L,5);
    acc = compacc(idx,cids);
    disp(['lambda=' num2str(ls(i)) ',seg acc=' num2str(acc)]);
    rs = [rs,rank(Z,1e-3*norm(Z,2))];
    accs = [accs,acc];
end

close all;
figure;
subplot(1,2,1);
plot(ls,accs);
xlabel('parameter \lambda');
ylabel('segmentation accuracy');
subplot(1,2,2);
plot(ls,rs);
xlabel('parameter \lambda');
ylabel('rank(Z)');
function [X,cids] = generate_data()
n = 200;
d = 10;
D = 100;
[U,~,~] = svd(rand(D));
cids = [];
U1 = U(:,1:d);
X = U1*rand(d,n);
cids = [cids,ones(1,n)];

for i=2:5
    R = orth(rand(D));
    U1 = R*U1;
    X = [X,U1*rand(d,n)];
    cids = [cids,i*ones(1,n)];
end
nX = size(X,2);
norm_x = sqrt(sum(X.^2,1));
norm_x = repmat(norm_x,D,1);
gn = norm_x.*randn(D,nX);
inds = rand(1,nX)<=0.3;
X(:,inds) = X(:,inds) + 0.3*gn(:,inds);

function [idx] = clu_ncut(L,K)
% this routine groups the data X into K subspaces by NCut
% inputs:
%       L -- an N*N affinity matrix, N is the number of data points
%       K -- the number of subpaces (i.e., clusters)
L = abs(L)+abs(L');

D = diag(sum(L,2).^(-1./2));
L = eye(size(L,1)) - D*L*D;
[U,~,~] = svd(L);

V = U(:,end-K+1:end);
idx = kmeans(V,K,'emptyaction','singleton','replicates',10,'display','off');
idx = idx';

function [acc] = compacc(Segmentation,RefSegmentation)
ngroups = length(unique(RefSegmentation));
if(size(RefSegmentation,2)==1)
    RefSegmentation=RefSegmentation';
end
if(size(Segmentation,2)==1)
    Segmentation=Segmentation';
end
Permutations = perms(1:ngroups);
miss = zeros(size(Permutations,1),size(Segmentation,1));
for k=1:size(Segmentation,1)
    for j=1:size(Permutations,1)
        miss(j,k) = sum(abs(Segmentation(k,:)-Permutations(j,RefSegmentation))>0.1);
    end
end

[miss,~] = min(miss,[],1);

acc = 1 - miss/length(Segmentation);